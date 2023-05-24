use anyhow::Result;
use mlua::prelude::*;
use mlua::UserData;
use perbase_lib::{
    par_granges::{self, RegionProcessor},
    position::pileup_position::PileupPosition,
    read_filter::ReadFilter,
};
use rust_htslib::bam::{self, pileup::Alignment, record::Record, Read};
use std::env;
use std::path::PathBuf;

struct LuaReadFilter {
    expression: String,
    lua: Lua,
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

// The actual implementation of `ReadFilter`
impl ReadFilter for LuaReadFilter {
    // Filter reads based SAM flags and mapping quality, true means pass
    #[inline]
    fn filter_read(&self, read: &Record, alignment: Option<&Alignment>) -> bool {
        let r = MyRecord(read, alignment.expect("always have alignment here"));

        self.lua
            .scope(|scope| {
                let globals = self.lua.globals();
                let r = scope
                    .create_nonstatic_userdata(r)
                    .expect("error creating user data");
                globals.set("read", r).expect("error setting read");

                let result: bool = self
                    .lua
                    // TODO: compile once, eval many
                    .load(self.expression.as_str())
                    .set_name("filter_read")
                    .expect("error setting chunk")
                    .eval()
                    .expect("error getting result");

                Ok(result)
            })
            .unwrap_or(false)
    }
}

struct BasicProcessor {
    // An indexed bamfile to query for the region we were passed
    bamfile: PathBuf,
    expression: String,
}

// Implementation of the `RegionProcessor` trait to process each region
impl RegionProcessor for BasicProcessor {
    type P = PileupPosition;

    // This function receives an interval to examine.
    fn process_region(&self, tid: u32, start: u32, stop: u32) -> Vec<Self::P> {
        let mut reader = bam::IndexedReader::from_path(&self.bamfile).expect("Indexed reader");
        let header = reader.header().to_owned();

        let rf = LuaReadFilter {
            lua: Lua::new(),
            expression: self.expression.to_string(),
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
        let result: Vec<PileupPosition> = reader
            .pileup()
            .flat_map(|p| {
                let pileup = p.expect("Extracted a pileup");
                // Verify that we are within the bounds of the chunk we are iterating on
                // Since pileup will pull reads that overhang edges.
                if pileup.pos() >= start && pileup.pos() < stop {
                    Some(PileupPosition::from_pileup(pileup, &header, &rf, None))
                } else {
                    None
                }
            })
            .collect();
        result
    }
}

fn main() -> Result<()> {
    let bam_path = env::args().nth(1).expect("bamfile");
    let expression = env::args().nth(2).expect("lua expression");

    if !expression.contains("return") {
        eprintln!("Expression '{}' must contain 'return'", expression);
        std::process::exit(1);
    }

    let basic_processor = BasicProcessor {
        bamfile: PathBuf::from(&bam_path),
        expression: String::from("") + expression.as_str(),
    };

    let par_granges_runner = par_granges::ParGranges::new(
        PathBuf::from(bam_path), // pass in bam
        None,                    // optional ref fasta
        None,                    // bedfile to narrow regions
        None,                    // optional bcf/vcf file to specify positions of interest
        true,                    // Merge any overlapping regions in the BED file
        None,                    // optional allowed number of threads, defaults to max
        None,                    // optional chunksize modification
        None, // optional modifier on the size of the channel for sending Positions
        basic_processor,
    );

    // Run the processor
    let receiver = par_granges_runner.process()?;
    // Pull the in-order results from the receiver channel
    receiver.into_iter().for_each(|p: PileupPosition| {
        // Note that the returned values are required to be `serde::Serialize`, so more fancy things
        // than just debug printing are doable.
        if p.depth > 0 {
            println!("p:{:?}", p);
        }
    });

    Ok(())
}
