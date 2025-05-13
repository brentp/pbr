use anyhow::{anyhow, Context, Result};
use bio::io::bed;
use rust_htslib::bam::{pileup::Pileup, HeaderView, self, Read};
use rust_lapper::{Interval, Lapper};
use mlua::Lua;
use perbase_lib::{
    par_granges::RegionProcessor,
    position::pileup_position::PileupPosition,
};
use crate::cached_faidx::CachedFaidx;
use crate::lua_filter::LuaReadFilter;

use std::path::PathBuf;

pub(crate) struct BasicProcessor {
    // An indexed bamfile to query for the region we were passed
    pub(crate) bamfile: PathBuf,
    pub(crate) expression: String,
    pub(crate) max_depth: u32,
    pub(crate) exclude_regions: Option<PathBuf>,
    pub(crate) mate_fix: bool,
    pub(crate) fasta_path: Option<PathBuf>,
    pub(crate) flanking: usize,
}

impl BasicProcessor {
    /// copied verbatim from perbase
    pub(crate) fn bed_to_intervals(
        header: &HeaderView,
        bed_file: &PathBuf,
        merge: bool,
    ) -> Result<Vec<Lapper<u32, ()>>> {
        let mut bed_reader = bed::Reader::from_file(bed_file)?;
        let mut intervals = vec![vec![]; header.target_count() as usize];
        for (i, record) in bed_reader.records().enumerate() {
            let record = record?;
            let tid = header
                .tid(record.chrom().as_bytes())
                .expect("Chromosome not found in BAM/CRAM header");
            let start = record
                .start()
                .try_into()
                .with_context(|| format!("BED record {} is invalid: unable to parse start", i))?;
            let stop = record
                .end()
                .try_into()
                .with_context(|| format!("BED record {} is invalid: unable to parse stop", i))?;
            if stop < start {
                return Err(anyhow!("BED record {} is invalid: stop < start", i));
            }
            intervals[tid as usize].push(Interval {
                start,
                stop,
                val: (),
            });
        }

        Ok(intervals
            .into_iter()
            .map(|ivs| {
                let mut lapper = Lapper::new(ivs);
                if merge {
                    lapper.merge_overlaps();
                }
                lapper
            })
            .collect())
    }
}

#[inline]
pub(crate) fn excluded(exclude_intervals: &Option<Vec<Lapper<u32, ()>>>, p: &Pileup) -> bool {
    match exclude_intervals {
        Some(ref intervals) => {
            let pos = p.pos();
            intervals[p.tid() as usize].find(pos, pos + 1).next().is_some()
        }
        None => false,
    }
}

impl RegionProcessor for BasicProcessor {
    type P = (PileupPosition, Option<String>);

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
        let result: Vec<(PileupPosition, Option<String>)> = p
            .flat_map(|p| {
                let pileup = p.expect("Extracted a pileup");
                // Verify that we are within the bounds of the chunk we are iterating on
                // Since pileup will pull reads that overhang edges.
                if pileup.pos() >= start
                    && pileup.pos() < stop
                    // and check if this position is excluded.
                    && !excluded(&exclude_intervals, &pileup)
                {
                    let pos = if self.mate_fix {
                        PileupPosition::from_pileup_mate_aware(
                            pileup, &header, &rf, None,
                        )
                    } else {
                        PileupPosition::from_pileup(pileup, &header, &rf, None)
                    };

                    let ref_seq = if let Some(fai) = &mut fai {
                        if self.flanking == 0 {
                            // When no flanking bases are requested, just get the single base
                            fai.fetch_seq_string(chrom, pos.pos as usize, pos.pos as usize + 1)
                                .unwrap_or_else(|_| ".".to_string())
                                .into()
                        } else {
                            let start = if pos.pos >= self.flanking as u32 {
                                pos.pos - self.flanking as u32
                            } else {
                                0
                            };
                            let end = pos.pos + self.flanking as u32 + 1;
                            
                            let left_padding = if pos.pos < self.flanking as u32 {
                                ".".repeat(self.flanking - pos.pos as usize)
                            } else {
                                String::new()
                            };
                            
                            let seq = fai
                                .fetch_seq_string(chrom, start as usize, end as usize)
                                .unwrap_or_else(|_| ".".repeat(2 * self.flanking + 1));
                            
                            let right_padding = if seq.len() < 2 * self.flanking + 1 {
                                ".".repeat(2 * self.flanking + 1 - seq.len())
                            } else {
                                String::new()
                            };
                            
                            Some(format!("{}{}{}", left_padding, seq, right_padding))
                        }
                    } else {
                        None
                    };

                    Some((pos, ref_seq))
                } else {
                    None
                }
            })
            .collect();
        result
    }
}
