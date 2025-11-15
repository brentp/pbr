use anyhow::{anyhow, Context, Result};
use bio::io::bed;
use rust_htslib::bam::{pileup::Pileup, HeaderView};
use rust_lapper::{Interval, Lapper};

use std::collections::HashSet;
use std::path::PathBuf;

pub(crate) struct BasicProcessor {
    // An indexed bamfile to query for the region we were passed
    pub(crate) bamfile: PathBuf,
    pub(crate) expression: String,
    pub(crate) max_depth: u32,
    pub(crate) exclude_regions: Option<PathBuf>,
    pub(crate) mate_fix: bool,
    pub(crate) fasta_path: Option<PathBuf>,
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
        let mut missing_chroms: HashSet<String> = HashSet::new();
        for (i, record) in bed_reader.records().enumerate() {
            let record = record?;
            let chrom = record.chrom();
            let Some(tid) = header.tid(chrom.as_bytes()) else {
                if missing_chroms.insert(chrom.to_string()) {
                    eprintln!(
                        "Chromosome {} not found in BAM/CRAM header; skipping intervals",
                        chrom
                    );
                }
                continue;
            };
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
            if p.tid() as usize >= intervals.len() {
                return false;
            }
            let ivs = &intervals[p.tid() as usize];
            let pos = p.pos();
            ivs.count(pos, pos + 1) > 0
        }
        None => false,
    }
}
