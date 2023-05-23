use anyhow::Result;
use mlua::prelude::*;
use mlua::UserData;
use perbase_lib::{
    par_granges::{self, RegionProcessor},
    position::pileup_position::PileupPosition,
    read_filter::ReadFilter,
};
use rust_htslib::bam::{self, record::Record, Read};
use std::env;
use std::path::PathBuf;

struct LuaReadFilter {
    expression: String,
    lua: Lua,
}

struct MyRecord<'a>(&'a Record);

impl<'a> UserData for MyRecord<'a> {
    fn add_fields<'lua, F: mlua::UserDataFields<'lua, Self>>(fields: &mut F) {
        fields.add_field_method_get("mapping_quality", |_, this| Ok(this.0.mapq()));
        fields.add_field_method_get("flags", |_, this| Ok(this.0.flags()));
        fields.add_field_method_get("tid", |_, this| Ok(this.0.tid()));
        fields.add_field_method_get("pos", |_, this| Ok(this.0.pos()));
        fields.add_field_method_get("start", |_, this| Ok(this.0.pos()));
        fields.add_field_method_get("stop", |_, this| Ok(this.0.cigar().end_pos()));
        // see:  https://github.com/rust-bio/rust-htslib/pull/393
        //fields.add_field_method_get("strand", |_, this| {
        //    Ok(String::from(this.0.strand().strand_symbol()))
        //});
        fields.add_field_method_get("insert_size", |_, this| Ok(this.0.insert_size()));
        fields.add_field_method_get("qname", |_, this| Ok(this.0.qname().to_owned()));
        fields.add_field_method_get("base_qualities", |_, this| {
            let quals = this.0.qual().to_owned();
            Ok(quals)
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
    fn filter_read(&self, read: &Record) -> bool {
        let r = MyRecord(read);

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

    let basic_processor = BasicProcessor {
        bamfile: PathBuf::from(&bam_path),
        expression: String::from("return ") + expression.as_str(),
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
        println!("p:{:?}", p);
    });

    Ok(())
}
