pub use rust_htslib::errors::{Error, Result};
use rust_htslib::faidx;
use std::path::Path;

/// CachedFaidx uses rust-htslib faidx reader
/// and caches the results to reduce disk access.
/// It does not do anything smart but should work well for
/// single consecutive bases as used in pbr.
pub struct CachedFaidx {
    faidx: faidx::Reader,
    cache: Vec<u8>,
    chrom: String,
    start: usize,
}

impl CachedFaidx {
    pub fn new<P: AsRef<Path>>(fasta_path: P) -> Result<Self> {
        let faidx = faidx::Reader::from_path(fasta_path)?;
        let cache = vec![0; 1000];
        Ok(CachedFaidx {
            faidx,
            cache,
            chrom: String::new(),
            start: 0,
        })
    }

    //pub fn n_seqs(&self) -> u64 {
    //    self.faidx.n_seqs()
    //}

    fn fetch_into_cache<N: AsRef<str>>(
        &mut self,
        chrom: N,
        start: usize,
        end: usize,
    ) -> Result<()> {
        let r = self.faidx.fetch_seq(chrom.as_ref(), start, end)?;
        self.chrom = String::from(chrom.as_ref());
        self.start = start;
        self.cache.clear();
        self.cache.extend_from_slice(r);
        Ok(())
    }

    #[allow(dead_code)]
    pub fn fetch_seq_string<N: AsRef<str> + std::cmp::PartialEq>(
        &mut self,
        chrom: N,
        start: usize,
        end: usize,
    ) -> Result<String> {
        let bytes = self.fetch_seq(chrom, start, end)?;
        Ok(std::str::from_utf8(bytes).unwrap().to_owned())
    }

    pub fn fetch_seq<N: AsRef<str> + std::cmp::PartialEq>(
        &mut self,
        chrom: N,
        start: usize,
        end: usize,
    ) -> Result<&[u8]> {
        if chrom.as_ref() == self.chrom
            && start >= self.start
            && end < self.start + self.cache.len()
        {
            let cstart = start - self.start;
            let cend = end - self.start;
            return Ok(&self.cache[cstart..cend + 1]);
        }
        self.fetch_into_cache(chrom, start, std::cmp::max(end, start + 1000))?;
        Ok(&self.cache[0..std::cmp::min(self.cache.len(), (end - start) + 1)])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn open_reader() -> CachedFaidx {
        CachedFaidx::new(format!("{}/test/test_cram.fa", env!("CARGO_MANIFEST_DIR")))
            .ok()
            .unwrap()
    }
    #[test]
    fn faidx_open() {
        open_reader();
    }

    #[test]
    fn faidx_read_chr_first_base() {
        let mut r = open_reader();

        let bseq = r.fetch_seq("chr1", 0, 0).unwrap();
        assert_eq!(bseq.len(), 1);
        assert_eq!(bseq, b"G");
    }

    #[test]
    fn faidx_read_chr_start() {
        let mut r = open_reader();

        let bseq = r.fetch_seq("chr1", 0, 9).unwrap();
        assert_eq!(bseq.len(), 10);
        assert_eq!(bseq, b"GGGCACAGCC");
    }

    #[test]
    fn faidx_read_chr_between() {
        let mut r = open_reader();

        let bseq = r.fetch_seq("chr1", 4, 14).unwrap();
        assert_eq!(bseq.len(), 11);
        assert_eq!(bseq, b"ACAGCCTCACC");

        let seq = r.fetch_seq_string("chr1", 4, 14).unwrap();
        assert_eq!(seq.len(), 11);
        assert_eq!(seq, "ACAGCCTCACC");
    }

    #[test]
    fn faidx_read_chr_end() {
        let mut r = open_reader();

        let bseq = r.fetch_seq("chr1", 110, 120).unwrap();
        assert_eq!(bseq.len(), 10);
        assert_eq!(bseq, b"CCCCTCCGTG");
    }

    #[test]
    fn faidx_read_twice_bytes() {
        let mut r = open_reader();
        let seq = r.fetch_seq("chr1", 110, 120).unwrap();
        assert_eq!(seq.len(), 10);
        assert_eq!(seq, b"CCCCTCCGTG");

        let seq = r.fetch_seq("chr1", 5, 9).unwrap();
        assert_eq!(seq.len(), 5);
        assert_eq!(seq, b"CAGCC");
    }

    #[test]
    fn faidx_position_too_large() {
        let mut r = open_reader();
        let position_too_large = i64::MAX as usize;
        let res = r.fetch_seq("chr1", position_too_large, position_too_large + 1);
        assert_eq!(res, Err(Error::FaidxPositionTooLarge));
    }

    #[test]
    fn open_many_readers() {
        for _ in 0..500_000 {
            let reader = open_reader();
            drop(reader);
        }
    }
}
