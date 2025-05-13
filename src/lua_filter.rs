use anyhow::Result;
use mlua::{Lua, Function, UserDataFields, UserDataMethods};
use rust_htslib::bam::{record::{Aux, Cigar, Record}, pileup::Alignment};
use perbase_lib::read_filter::ReadFilter;

pub struct LuaReadFilter<'a> {
    pub(crate) lua: &'a Lua,
    pub(crate) filter_func: Function,
}

impl<'a> LuaReadFilter<'a> {
    // Create a new LuaReadFilter instance with the given expression
    pub fn new(expression: &str, lua: &'a Lua) -> Result<Self> {
        let filter_func = lua.load(expression).into_function()?;
        lua.register_userdata_type::<Record>(|reg| {
            reg.add_field_method_get("mapping_quality", |_, this| Ok(this.mapq()));
            reg.add_field_method_get("flags", |_, this| Ok(this.flags()));
            reg.add_field_method_get("tid", |_, this| Ok(this.tid()));
            reg.add_field_method_get("start", |_, this| Ok(this.pos()));
            reg.add_field_method_get("stop", |_, this| Ok(this.cigar().end_pos()));
            reg.add_field_method_get("length", |_, this| Ok(this.seq_len()));
            reg.add_field_method_get("insert_size", |_, this| Ok(this.insert_size()));
            reg.add_field_method_get("qname", |_, this| {
                let q = this.qname();
                Ok(std::str::from_utf8(q).unwrap_or("").to_string())
            });
            reg.add_field_method_get("sequence", |_, this| {
                let seq = this.seq();
                Ok(std::str::from_utf8(&seq.as_bytes())
                    .unwrap_or("")
                    .to_string())
            });
            reg.add_field_method_get("strand", |_, this| {
                let flags = this.flags();
                // Forward strand logic
                let is_forward = ((flags & 1) == 0 && (flags & 16) == 0) || 
                                ((flags & 1) != 0 && ((flags & 99) == 99 || (flags & 147) == 147));
                // Reverse strand logic
                let is_reverse = ((flags & 1) == 0 && (flags & 16) != 0) || 
                                ((flags & 1) != 0 && ((flags & 83) == 83 || (flags & 163) == 163));
                
                Ok(if is_forward { 1 } else if is_reverse { -1 } else { 0 })
            });
            reg.add_function("qpos", |_, this: mlua::AnyUserData| {
                let r: Result<usize, mlua::Error> = this.named_user_value("qpos");
                r
            });
            reg.add_field_function_get("bq", |_, this: mlua::AnyUserData| {
                let qpos: usize = match this.named_user_value::<usize>("qpos") {
                    Ok(qpos) => qpos,
                    Err(_) => {
                        return Ok(-1);
                    }
                };
                this.borrow_scoped::<Record, i32>(|r| match qpos {
                    usize::MAX => -1,
                    _ => r.qual()[qpos] as i32,
                })
            });
            reg.add_field_function_get("distance_from_5prime", |_, this| {
                let qpos: usize = match this.named_user_value("qpos") {
                    Ok(qpos) => qpos,
                    Err(_) => {
                        return Ok(-1);
                    }
                };
                this.borrow_scoped::<Record, i32>(|r| {
                    if r.is_reverse() {
                        r.seq_len() as i32 - qpos as i32
                    } else {
                        qpos as i32
                    }
                })
            });
            reg.add_field_function_get("distance_from_3prime", |_, this| {
                let qpos: usize = match this.named_user_value("qpos") {
                    Ok(qpos) => qpos,
                    Err(_) => {
                        return Ok(usize::MAX);
                    }
                };
                this.borrow_scoped::<Record, usize>(|r| {
                    if r.is_reverse() {
                        qpos
                    } else {
                        r.seq_len() - qpos
                    }
                })
            });

            reg.add_method("n_proportion_3_prime", |_, this, n_bases: usize| {
                let seq = this.seq();
                let mut count = 0;
                let reverse = this.is_reverse();
                for i in 0..n_bases {
                    let base =
                        seq[if reverse { i } else { seq.len() - 1 - i }].to_ascii_uppercase();
                    if base == b'N' {
                        count += 1;
                    }
                }
                Ok(count as f64 / n_bases as f64)
            });

            reg.add_method("n_proportion_5_prime", |_, this, n_bases: usize| {
                let seq = this.seq();
                let mut count = 0;
                let reverse = this.is_reverse();
                for i in 0..n_bases {
                    let base =
                        seq[if reverse { seq.len() - 1 - i } else { i }].to_ascii_uppercase();
                    if base == b'N' {
                        count += 1;
                    }
                }
                Ok(count as f64 / n_bases as f64)
            });

            reg.add_field_method_get("indel_count", |_, this| {
                let cigar = this.cigar();
                let mut count = 0;
                for op in cigar.iter() {
                    match op {
                        Cigar::Ins(_) | Cigar::Del(_) => {
                            count += 1;
                        }
                        _ => {}
                    }
                }
                Ok(count)
            });

            reg.add_field_method_get("soft_clips_3_prime", |_, this| {
                let cigar = this.cigar();
                if this.is_reverse() {
                    Ok(cigar.leading_softclips())
                } else {
                    Ok(cigar.trailing_softclips())
                }
            });
            reg.add_field_method_get("soft_clips_5_prime", |_, this| {
                let cigar = this.cigar();
                if this.is_reverse() {
                    Ok(cigar.trailing_softclips())
                } else {
                    Ok(cigar.leading_softclips())
                }
            });
            reg.add_field_method_get("average_base_quality", |_, this| {
                let qual = this.qual();
                let sum = qual.iter().map(|q| *q as u64).sum::<u64>();
                let count = qual.len();
                Ok(sum as f64 / count as f64)
            });

            reg.add_method("tag", |lua, this: &Record, tag: String| {
                let tag = tag.as_bytes();
                let aux = this.aux(tag).map_err(mlua::Error::external)?;
                let lua_val: mlua::Value = match aux {
                    Aux::Char(v) => mlua::Value::String(lua.create_string(&[v])?),
                    Aux::I8(v) => mlua::Value::Number(v as f64),
                    Aux::U8(v) => mlua::Value::Number(v as f64),
                    Aux::I16(v) => mlua::Value::Number(v as f64),
                    Aux::U16(v) => mlua::Value::Number(v as f64),
                    Aux::I32(v) => mlua::Value::Number(v as f64),
                    Aux::U32(v) => mlua::Value::Number(v as f64),
                    Aux::Float(v) => mlua::Value::Number(v as f64),
                    Aux::Double(v) => mlua::Value::Number(v as f64),
                    Aux::String(v) => mlua::Value::String(lua.create_string(&v)?),
                    Aux::ArrayFloat(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(f32::NAN) as f32);
                        }
                        mlua::Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayI32(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(i32::MIN) as i32);
                        }
                        mlua::Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayI8(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(i8::MIN) as i8);
                        }
                        mlua::Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayU8(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(u8::MIN) as u8);
                        }
                        mlua::Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayU16(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(u16::MIN) as u16);
                        }
                        mlua::Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayU32(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(u32::MIN) as u32);
                        }
                        mlua::Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayI16(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(i16::MIN) as i16);
                        }
                        mlua::Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::HexByteArray(v) => {
                        let lstr = String::from_utf8_lossy(v.as_bytes()).to_string();
                        mlua::Value::String(lua.create_string(&lstr)?)
                    }
                };
                Ok(Some(lua_val))
            })
        })?;
        Ok(Self { lua, filter_func })
    }
}

impl<'a> ReadFilter for LuaReadFilter<'a> {
    /// Filter reads based user expression.
    #[inline]
    fn filter_read(&self, read: &Record, alignment: Option<&Alignment>) -> bool {
        let r = self.lua.scope(|scope| {
            let globals = self.lua.globals();
            let ud = scope.create_any_userdata_ref(read)?;
            ud.set_named_user_value("qpos", alignment.unwrap().qpos().unwrap_or(usize::MAX))?;

            globals.set("read", ud).expect("error setting read");

            self.filter_func.call::<bool>(())
        });

        match r {
            Ok(r) => r,
            Err(e) => {
                eprintln!("Error evaluating expression: {}", e);
                false
            }
        }
    }
} 