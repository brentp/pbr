#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod cached_faidx;
mod processor;

use anyhow::Result;
use cached_faidx::CachedFaidx;
use clap::Parser;
use processor::{excluded, BasicProcessor};

use mlua::prelude::*;
use mlua::{Function, UserData};
use perbase_lib::{
    par_granges::{self, RegionProcessor},
    position::pileup_position::PileupPosition,
    read_filter::ReadFilter,
};
use rust_htslib::bam::{self, pileup::Alignment, record::Record, Read};
use std::path::PathBuf;

struct LuaReadFilter<'a> {
    lua: &'a Lua,
    filter_func: Function<'a>,
}

impl<'a> LuaReadFilter<'a> {
    // Create a new LuaReadFilter instance with the given expression
    fn new(expression: &str, lua: &'a Lua) -> Result<Self> {
        let filter_func = lua.load(expression).into_function()?;
        lua.register_userdata_type::<MyRecord>(|reg| {
            reg.add_field_method_get("mapping_quality", |_, this| Ok(this.0.mapq()));
            reg.add_field_method_get("flags", |_, this| Ok(this.0.flags()));
            reg.add_field_method_get("tid", |_, this| Ok(this.0.tid()));
            reg.add_field_method_get("start", |_, this| Ok(this.0.pos()));
            reg.add_field_method_get("stop", |_, this| Ok(this.0.cigar().end_pos()));
            /*
            p.add_field_method_get("qpos", |_, this| Ok(this.1.qpos()));
            p.add_field_method_get("bq", |_, this| {
                if let Some(qpos) = this.1.qpos() {
                    let qual = this.0.qual()[qpos];
                    Ok(qual as i32)
                } else {
                    //Err(mlua::Error::RuntimeError("qpos is None".to_string()))
                    Ok(-1)
                }
            });
            */
        })?;

        Ok(Self { lua, filter_func })
    }
}

struct Pile<'a>(&'a PileupPosition);
impl<'a> UserData for Pile<'a> {
    fn add_fields<'lua, F: mlua::UserDataFields<'lua, Self>>(fields: &mut F) {
        fields.add_field_method_get("depth", |_, this| Ok(this.0.depth));
        fields.add_field_method_get("a", |_, this| Ok(this.0.a));
        fields.add_field_method_get("c", |_, this| Ok(this.0.c));
        fields.add_field_method_get("g", |_, this| Ok(this.0.g));
        fields.add_field_method_get("t", |_, this| Ok(this.0.t));
        fields.add_field_method_get("n", |_, this| Ok(this.0.n));
        fields.add_field_method_get("fail", |_, this| Ok(this.0.fail));
        fields.add_field_method_get("ins", |_, this| Ok(this.0.ins));
        fields.add_field_method_get("del", |_, this| Ok(this.0.del));
        fields.add_field_method_get("ref_skip", |_, this| Ok(this.0.ref_skip));
        fields.add_field_method_get("pos", |_, this| Ok(this.0.pos));
    }
}

struct MyRecord<'a>(&'a Record);
impl<'a> UserData for MyRecord<'a> {
    fn add_fields<'lua, F: mlua::UserDataFields<'lua, Self>>(fields: &mut F) {
        fields.add_field_method_get("mapping_quality", |_, this| Ok(this.0.mapq()));
        fields.add_field_method_get("flags", |_, this| Ok(this.0.flags()));
        fields.add_field_method_get("tid", |_, this| Ok(this.0.tid()));
        fields.add_field_method_get("start", |_, this| Ok(this.0.pos()));
        fields.add_field_method_get("stop", |_, this| Ok(this.0.cigar().end_pos()));
        /*
        fields.add_field_method_get("qpos", |_, this| Ok(this.1.qpos()));
        fields.add_field_method_get("distance_from_left_end", |_, this| {
            Ok(match this.1.qpos() {
                Some(qpos) => qpos as i32,
                None => -1,
            })
        });
        fields.add_field_method_get("distance_from_right_end", |_, this| {
            if this.1.qpos().is_none() {
                return Ok(-1);
            }
            let len = this.0.seq_len();
            Ok((len - this.1.qpos().unwrap_or(0)) as i32)
        });
        */

        // see:  https://github.com/rust-bio/rust-htslib/pull/393
        //fields.add_field_method_get("strand", |_, this| {
        //    Ok(String::from(this.0.strand().strand_symbol()))
        //});
        fields.add_field_method_get("length", |_, this| Ok(this.0.seq_len()));
        fields.add_field_method_get("insert_size", |_, this| Ok(this.0.insert_size()));
        fields.add_field_method_get("qname", |_, this| {
            let q = this.0.qname();
            Ok(std::str::from_utf8(q)
                .unwrap_or("BAD_READ_NAME")
                .to_string())
        });
        fields.add_field_method_get("base_qualities", |_, this| {
            let quals = this.0.qual().to_owned();
            Ok(quals)
        });
        /*
        fields.add_field_method_get("bq", |_, this| {
            if let Some(qpos) = this.1.qpos() {
                let qual = this.0.qual()[qpos];
                Ok(qual as i32)
            } else {
                //Err(mlua::Error::RuntimeError("qpos is None".to_string()))
                Ok(-1)
            }
        });
        */
        fields.add_field_method_get("sequence", |_, this| {
            let seq = this.0.seq();
            Ok(String::from(unsafe {
                std::str::from_utf8_unchecked(&seq.as_bytes())
            }))
        });
        /*
        fields.add_field_method_get("cigar", |_, this| {
            let cigar = this.0.cigar();
            Ok(cigar)
        });
        */
        fields.add_field_method_get("cigar_string", |_, this| {
            let cigar = this.0.cigar();
            Ok(cigar.to_string())
        });
    }
}

impl<'a> ReadFilter for LuaReadFilter<'a> {
    /// Filter reads based user expression.
    #[inline]
    fn filter_read(&self, read: &Record, _alignment: Option<&Alignment>) -> bool {
        let rec = MyRecord(read);
        let r = self.lua.scope(|scope| {
            let globals = self.lua.globals();
            let ud = scope.create_userdata_ref(&rec)?;

            globals.set("read", ud).expect("error setting read");

            self.filter_func.call::<_, bool>(())
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
                    let p = scope
                        .create_nonstatic_userdata(Pile(p))
                        .expect("error creating user data");
                    globals.set("pile", p).expect("error setting pile");

                    pile_expression.call::<_, bool>(())
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
                    .create_nonstatic_userdata(Pile(&pileup_position))
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
