#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod processor;

use processor::{excluded, BasicProcessor};

use anyhow::Result;
use clap::Parser;

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

        Ok(Self {
            lua: lua,
            filter_func,
        })
    }
}

struct MyRecord<'a>(&'a Record, &'a Alignment<'a>);

impl<'a> UserData for MyRecord<'a> {
    fn add_fields<'lua, F: mlua::UserDataFields<'lua, Self>>(fields: &mut F) {
        fields.add_field_method_get("mapping_quality", |_, this| Ok(this.0.mapq()));
        fields.add_field_method_get("flags", |_, this| Ok(this.0.flags()));
        fields.add_field_method_get("tid", |_, this| Ok(this.0.tid()));
        fields.add_field_method_get("pos", |_, this| Ok(this.0.pos()));
        fields.add_field_method_get("start", |_, this| Ok(this.0.pos()));
        fields.add_field_method_get("stop", |_, this| Ok(this.0.cigar().end_pos()));
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
        fields.add_field_method_get("bq", |_, this| {
            if let Some(qpos) = this.1.qpos() {
                let qual = this.0.qual()[qpos as usize];
                Ok(qual as i32)
            } else {
                //Err(mlua::Error::RuntimeError("qpos is None".to_string()))
                Ok(-1)
            }
        });
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
    fn filter_read(&self, read: &Record, alignment: Option<&Alignment>) -> bool {
        let r = MyRecord(read, alignment.expect("always have alignment here"));

        let r = self.lua.scope(|scope| {
            let globals = self.lua.globals();
            let r = scope
                .create_nonstatic_userdata(r)
                .expect("error creating user data");
            globals.set("read", r).expect("error setting read");

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
        let header = reader.header().to_owned();
        let lua = Lua::new();

        let rf = LuaReadFilter::new(&self.expression, &lua).expect(
            format!(
                "error creating lua read filter with expression {}",
                &self.expression
            )
            .as_str(),
        );

        let exclude_intervals = if let Some(regions_bed) = &self.exclude_regions {
            Some(Self::bed_to_intervals(&header, regions_bed, true).expect("BED file"))
        } else {
            None
        };

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
        p.set_max_depth(self.max_depth);
        let result: Vec<PileupPosition> = p
            .flat_map(|p| {
                let pileup = p.expect("Extracted a pileup");
                // Verify that we are within the bounds of the chunk we are iterating on
                // Since pileup will pull reads that overhang edges.
                if pileup.pos() >= start
                    && pileup.pos() < stop
                    // and check if this position is excluded.
                    && !excluded(&exclude_intervals, &pileup)
                {
                    Some(PileupPosition::from_pileup(pileup, &header, &rf, None))
                } else {
                    None
                }
            })
            .collect();
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
    };

    let par_granges_runner = par_granges::ParGranges::new(
        PathBuf::from(opts.bam_path), // pass in bam
        opts.fasta,                   // optional ref fasta
        opts.bedfile,                 // bedfile to narrow regions
        None,                         // optional bcf/vcf file to specify positions of interest
        true,                         // Merge any overlapping regions in the BED file
        Some(opts.threads),           // optional allowed number of threads, defaults to max
        None,                         // optional chunksize modification
        None, // optional modifier on the size of the channel for sending Positions
        basic_processor,
    );

    // Run the processor
    let receiver = par_granges_runner.process()?;
    println!("#chrom\tpos\tdepth\ta\tc\tg\tt\tn");
    // Pull the in-order results from the receiver channel
    receiver.into_iter().for_each(|p: PileupPosition| {
        // Note that the returned values are required to be `serde::Serialize`, so more fancy things
        // than just debug printing are doable.
        if p.depth > 0 {
            //p:PileupPosition { ref_seq: "chr2", pos: 196, ref_base: None, depth: 1, a: 1, c: 0, g: 0, t: 0, n: 0, ins: 0, del: 0, ref_skip: 0, fail: 1, near_max_depth: false }
            println!(
                "{chrom}\t{pos}\t{depth}\t{a}\t{c}\t{g}\t{t}\t{n}",
                chrom = p.ref_seq,
                pos = p.pos,
                depth = p.depth,
                a = p.a,
                c = p.c,
                g = p.g,
                t = p.t,
                n = p.n
            );
        }
    });

    Ok(())
}
