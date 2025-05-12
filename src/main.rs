#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod cached_faidx;
mod processor;

use anyhow::Result;
use cached_faidx::CachedFaidx;
use clap::Parser;
use processor::{excluded, BasicProcessor};

use mlua::prelude::*;
use mlua::{Function, Value};
use perbase_lib::{
    par_granges::{self, RegionProcessor},
    position::pileup_position::PileupPosition,
    read_filter::ReadFilter,
};
use rust_htslib::bam::{
    self,
    pileup::Alignment,
    record::{Aux, Cigar, Record},
    Read,
};
use std::path::PathBuf;

struct LuaReadFilter<'a> {
    lua: &'a Lua,
    filter_func: Function,
}

impl<'a> LuaReadFilter<'a> {
    // Create a new LuaReadFilter instance with the given expression
    fn new(expression: &str, lua: &'a Lua) -> Result<Self> {
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
            reg.add_function("qpos", |_, this: mlua::AnyUserData| {
                let r: Result<usize, LuaError> = this.named_user_value("qpos");
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
                let aux = this.aux(tag).map_err(LuaError::external)?;
                let lua_val: Value = match aux {
                    Aux::Char(v) => Value::String(lua.create_string(&[v])?),
                    Aux::I8(v) => Value::Number(v as f64),
                    Aux::U8(v) => Value::Number(v as f64),
                    Aux::I16(v) => Value::Number(v as f64),
                    Aux::U16(v) => Value::Number(v as f64),
                    Aux::I32(v) => Value::Number(v as f64),
                    Aux::U32(v) => Value::Number(v as f64),
                    Aux::Float(v) => Value::Number(v as f64),
                    Aux::Double(v) => Value::Number(v as f64),
                    Aux::String(v) => Value::String(lua.create_string(&v)?),
                    Aux::ArrayFloat(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(f32::NAN) as f32);
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayI32(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(i32::MIN) as i32);
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayI8(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(i8::MIN) as i8);
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayU8(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(u8::MIN) as u8);
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayU16(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(u16::MIN) as u16);
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayU32(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(u32::MIN) as u32);
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::ArrayI16(v) => {
                        let mut arr = Vec::new();
                        for i in 0..v.len() {
                            arr.push(v.get(i).unwrap_or(i16::MIN) as i16);
                        }
                        Value::Table(lua.create_sequence_from(arr)?)
                    }
                    Aux::HexByteArray(v) => {
                        let lstr = String::from_utf8_lossy(v.as_bytes()).to_string();
                        Value::String(lua.create_string(&lstr)?)
                    }
                };
                Ok(Some(lua_val))
            })
        })?;
        Ok(Self { lua, filter_func })
    }
}

fn register_pile(lua: &Lua) -> mlua::Result<()> {
    lua.register_userdata_type::<PileupPosition>(|reg| {
        reg.add_field_method_get("depth", |_, this| Ok(this.depth));
        reg.add_field_method_get("a", |_, this| Ok(this.a));
        reg.add_field_method_get("c", |_, this| Ok(this.c));
        reg.add_field_method_get("g", |_, this| Ok(this.g));
        reg.add_field_method_get("t", |_, this| Ok(this.t));
        reg.add_field_method_get("n", |_, this| Ok(this.n));
        reg.add_field_method_get("fail", |_, this| Ok(this.fail));
        reg.add_field_method_get("ins", |_, this| Ok(this.ins));
        reg.add_field_method_get("del", |_, this| Ok(this.del));
        reg.add_field_method_get("ref_skip", |_, this| Ok(this.ref_skip));
        reg.add_field_method_get("pos", |_, this| Ok(this.pos));
    })
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

impl RegionProcessor for BasicProcessor {
    type P = PileupPosition;

    // This function receives an interval to examine.
    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P> {
        let mut reader = bam::IndexedReader::from_path(&self.bamfile).expect("Indexed reader");
        let mut fai = if let Some(fasta) = &self.fasta_path {
            reader.set_reference(fasta).expect("reference");
            Some(CachedFaidx::new(fasta).expect("error reading fasta"))
        } else {
            None
        };

        let header = reader.header().to_owned();
        let lua = Lua::new();

        let rf = LuaReadFilter::new(&self.expression, &lua).unwrap_or_else(|_| {
            panic!(
                "error creating lua read filter with expression {}",
                &self.expression
            )
        });

        let exclude_intervals = self.exclude_regions.as_ref().map(|regions_bed| {
            Self::bed_to_intervals(&header, regions_bed, true).expect("BED file")
        });

        let string_count = rf
            .lua
            .create_function(|_, (haystack, needle): (String, String)| {
                assert!(needle.len() == 1);
                let needle = needle.chars().next().unwrap();
                Ok(haystack.chars().filter(|c| *c == needle).count())
            })
            .expect("eror creating function");
        rf.lua
            .globals()
            .set("string_count", string_count)
            .expect("error setting string_count");

        // fetch the region
        reader.fetch((tid, start, stop)).expect("Fetched ROI");
        // Walk over pileups
        let mut p = reader.pileup();
        let chrom = unsafe { std::str::from_utf8_unchecked(header.target_names()[tid as usize]) };
        p.set_max_depth(self.max_depth);
        let mut result: Vec<PileupPosition> = p
            .flat_map(|p| {
                let pileup = p.expect("Extracted a pileup");
                // Verify that we are within the bounds of the chunk we are iterating on
                // Since pileup will pull reads that overhang edges.
                if pileup.pos() >= start
                    && pileup.pos() < stop
                    // and check if this position is excluded.
                    && !excluded(&exclude_intervals, &pileup)
                {
                    if self.mate_fix {
                        Some(PileupPosition::from_pileup_mate_aware(
                            pileup, &header, &rf, None,
                        ))
                    } else {
                        Some(PileupPosition::from_pileup(pileup, &header, &rf, None))
                    }
                } else {
                    None
                }
            })
            .collect();
        if let Some(fai) = &mut fai {
            result.iter_mut().for_each(|p| {
                let s = fai
                    .fetch_seq(chrom, p.pos as usize, (p.pos + 1) as usize)
                    .expect("error extracting reference base");
                p.ref_base = Some(s[0] as char);
            });
        }
        result
    }
}

#[derive(Parser, Default, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(help = "Path to the bamfile")]
    bam_path: PathBuf,
    #[clap(help = "Lua expression to evaluate")]
    expression: String,
    #[clap(short, long, default_value = "2", help = "Number of threads to use")]
    threads: usize,
    #[clap(
        short,
        long,
        default_value_t = 100000,
        help = "maximum depth in the pileup"
    )]
    max_depth: u32,
    #[clap(short, long, help = "optional path to the BED of include regions")]
    bedfile: Option<PathBuf>,
    #[clap(short, long, help = "optional path to the reference fasta file")]
    fasta: Option<PathBuf>,
    #[clap(short, long, help = "optional path to BED of exclude regions")]
    exclude: Option<PathBuf>,

    #[clap(
        long,
        help = "adjust depth to not double count overlapping mates",
        long_help = "note that for now this is much slower than the default"
    )]
    mate_fix: bool,

    #[clap(short, long, help = "optional expression required for the pileup")]
    pile_expression: Option<String>,
}

fn main() -> Result<()> {
    let opts = Args::parse();

    if !opts.expression.contains("return") {
        eprintln!("Expression '{}' must contain 'return'", opts.expression);
        std::process::exit(1);
    }

    let basic_processor = BasicProcessor {
        bamfile: PathBuf::from(&opts.bam_path),
        expression: String::from("") + opts.expression.as_str(),
        max_depth: opts.max_depth,
        exclude_regions: opts.exclude,
        mate_fix: opts.mate_fix,
        fasta_path: opts.fasta.clone(),
    };

    let par_granges_runner = par_granges::ParGranges::new(
        opts.bam_path,      // pass in bam
        opts.fasta,         // optional ref fasta
        opts.bedfile,       // bedfile to narrow regions
        None,               // optional bcf/vcf file to specify positions of interest
        true,               // Merge any overlapping regions in the BED file
        Some(opts.threads), // optional allowed number of threads, defaults to max
        None,               // optional chunksize modification
        None,               // optional modifier on the size of the channel for sending Positions
        basic_processor,
    );

    let pile_lua = Lua::new();
    register_pile(&pile_lua)?;

    let pile_expression: Option<Function> = if let Some(expression) = opts.pile_expression {
        Some(pile_lua.load(expression.as_str()).into_function()?)
    } else {
        None
    };

    // Run the processor
    let receiver = par_granges_runner.process()?;
    println!("#chrom\tpos0\tref_base\tdepth\ta\tc\tg\tt\tn");
    // Pull the in-order results from the receiver channel
    receiver
        .into_iter()
        .filter(|p| p.depth > 0)
        // filter on the pile expression
        .filter(|p| {
            if let Some(pile_expression) = &pile_expression {
                let r = pile_lua.scope(|scope| {
                    let globals = pile_lua.globals();
                    let ud = scope.create_any_userdata_ref(p)?;
                    globals.set("pile", ud).expect("error setting pile");

                    pile_expression.call::<bool>(())
                });
                match r {
                    Ok(r) => r,
                    Err(e) => {
                        eprintln!("Error evaluating expression: {}", e);
                        std::process::exit(1);
                    }
                }
            } else {
                true
            }
        })
        .for_each(|p: PileupPosition| {
            //p:PileupPosition { ref_seq: "chr2", pos: 196, ref_base: None, depth: 1, a: 1, c: 0, g: 0, t: 0, n: 0, ins: 0, del: 0, ref_skip: 0, fail: 1, near_max_depth: false }
            println!(
                "{chrom}\t{pos}\t{ref_base}\t{depth}\t{a}\t{c}\t{g}\t{t}\t{n}",
                chrom = p.ref_seq,
                pos = p.pos,
                depth = p.depth,
                ref_base = p.ref_base.unwrap_or('.'),
                a = p.a,
                c = p.c,
                g = p.g,
                t = p.t,
                n = p.n
            );
        });

    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;
    use mlua::Lua;
    use rust_htslib::bam;
    use rust_htslib::bam::pileup::Pileup;
    use rust_htslib::bam::record::Record;
    use rust_htslib::bam::{header::HeaderRecord, Header, HeaderView, IndexedReader, Read};
    use tempfile::NamedTempFile;

    #[test]
    fn test_read_bq() -> Result<()> {
        // Create a header with chr1
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", &1000000u32);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        // Create a test BAM record using SAM format
        let record = Record::from_sam(
            &header_view,
            b"test_read\t0\tchr1\t100\t30\t4M\t*\t0\t0\tACGT\t&&&&\tRG:Z:test", // Use standard ASCII for qual
        )
        .expect("Failed to create record from SAM");

        // Write record to a temporary BAM file
        let tmp = NamedTempFile::new()?;
        let path = tmp.path();
        {
            let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam)?;
            writer.write(&record)?;
        }
        bam::index::build(path, None, bam::index::Type::Bai, 1)?;

        // Create pileup
        let mut reader = IndexedReader::from_path(path)?;
        reader.fetch(("chr1", 100, 101))?; // Fetch the position of the record
        let mut pileups = reader.pileup();
        let pileup: Pileup = pileups.next().unwrap().expect("Failed to get pileup");

        // Get the first alignment from the pileup
        let alignment = pileup
            .alignments()
            .next()
            .expect("No alignment found in pileup");

        let lua = Lua::new();
        let rf = LuaReadFilter::new(
            "return read.bq > 0 and read.distance_from_5prime == 0 and read.distance_from_3prime > 0",
            &lua,
        )?; // Example expression

        // Test the bq functionality using the alignment from the pileup
        let result = rf.filter_read(&alignment.record(), Some(&alignment));
        assert!(result);

        Ok(())
    }

    #[test]
    fn test_pileup_position() -> mlua::Result<()> {
        let pileup_position = PileupPosition {
            depth: 10,
            a: 1,
            c: 2,
            g: 3,
            t: 4,
            n: 5,
            fail: 6,
            ins: 7,
            del: 8,
            ref_skip: 9,
            pos: 10,
            ..Default::default()
        };

        let lua = Lua::new();
        register_pile(&lua)?;
        let globals = lua.globals();
        for (expected, expression) in [
            (true, "pile.g > 3"),
            (true, "pile.a > 0"),
            (false, "pile.a > 10"),
            (false, "pile.ref_skip == 100"),
            (true, "pile.ref_skip == 9"),
        ] {
            eprintln!("Testing expression: {}", expression);
            lua.scope(|scope| {
                let p = scope
                    .create_any_userdata_ref(&pileup_position)
                    .expect("error creating user data");
                globals.set("pile", p)?;
                let f = lua
                    .load(&(String::from("return ") + expression))
                    .into_function()?;
                let result: bool = f.call(())?;
                Ok(result == expected)
            })
            .expect("error evaluating expression");
        }
        Ok(())
    }
}
