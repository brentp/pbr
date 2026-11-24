#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod cached_faidx;
mod processor;

use anyhow::Result;
use cached_faidx::CachedFaidx;
use clap::Parser;
use processor::{excluded, BasicProcessor};

use mlua::prelude::*;
use mlua::{AnyUserData, Compiler, Function, UserData, UserDataFields, UserDataMethods, Value};
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
use std::cell::Cell;
use std::path::PathBuf;
use std::sync::Arc;

struct ReadSlot {
    record: Cell<*const Record>,
    qpos: Cell<usize>,
}

// SAFETY: ReadSlot is strictly confined to the single thread owning its Lua VM.
// The raw pointer is only non-null during a synchronous filter_read call and is
// always cleared by ReadSlotGuard on return or unwind.
unsafe impl Send for ReadSlot {}
unsafe impl Sync for ReadSlot {}

struct ReadSlotGuard<'a>(&'a ReadSlot);

impl Drop for ReadSlotGuard<'_> {
    #[inline(always)]
    fn drop(&mut self) {
        self.0.record.set(std::ptr::null());
        self.0.qpos.set(usize::MAX);
    }
}

#[derive(Clone)]
struct LuaRead(Arc<ReadSlot>);

impl LuaRead {
    #[inline(always)]
    fn with_record<T>(&self, f: impl FnOnce(&Record) -> mlua::Result<T>) -> mlua::Result<T> {
        let ptr = self.0.record.get();
        if ptr.is_null() {
            return Err(LuaError::RuntimeError(
                "read reference is no longer valid".to_string(),
            ));
        }
        // SAFETY: ptr points to a valid &Record for the duration of filter_read.
        f(unsafe { &*ptr })
    }
}

impl UserData for LuaRead {
    fn add_fields<F: UserDataFields<Self>>(fields: &mut F) {
        fields.add_field_method_get("mapping_quality", |_, this| {
            this.with_record(|r| Ok(r.mapq()))
        });
        fields.add_field_method_get("flags", |_, this| this.with_record(|r| Ok(r.flags())));
        fields.add_field_method_get("tid", |_, this| this.with_record(|r| Ok(r.tid())));
        fields.add_field_method_get("start", |_, this| this.with_record(|r| Ok(r.pos())));
        fields.add_field_method_get("stop", |_, this| {
            this.with_record(|r| Ok(r.cigar().end_pos()))
        });
        fields.add_field_method_get("length", |_, this| this.with_record(|r| Ok(r.seq_len())));
        fields.add_field_method_get("insert_size", |_, this| {
            this.with_record(|r| Ok(r.insert_size()))
        });
        fields.add_field_method_get("qname", |_, this| {
            this.with_record(|r| {
                let q = r.qname();
                Ok(std::str::from_utf8(q).unwrap_or("").to_string())
            })
        });
        fields.add_field_method_get("sequence", |_, this| {
            this.with_record(|r| {
                let seq = r.seq();
                Ok(std::str::from_utf8(&seq.as_bytes())
                    .unwrap_or("")
                    .to_string())
            })
        });
        fields.add_field_method_get("bq", |_, this| {
            let qpos = this.0.qpos.get();
            if qpos == usize::MAX {
                return Ok(-1);
            }
            this.with_record(|r| Ok(r.qual()[qpos] as i32))
        });
        fields.add_field_method_get("distance_from_5prime", |_, this| {
            let qpos = this.0.qpos.get();
            if qpos == usize::MAX {
                return Ok(-1);
            }
            this.with_record(|r| {
                if r.is_reverse() {
                    Ok(r.seq_len() as i32 - qpos as i32)
                } else {
                    Ok(qpos as i32)
                }
            })
        });
        fields.add_field_method_get("distance_from_3prime", |_, this| {
            let qpos = this.0.qpos.get();
            if qpos == usize::MAX {
                return Ok(usize::MAX);
            }
            this.with_record(|r| {
                if r.is_reverse() {
                    Ok(qpos)
                } else {
                    Ok(r.seq_len() - qpos)
                }
            })
        });
        fields.add_field_method_get("indel_count", |_, this| {
            this.with_record(|r| {
                let cigar = r.cigar();
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
            })
        });
        fields.add_field_method_get("soft_clips_3_prime", |_, this| {
            this.with_record(|r| {
                let cigar = r.cigar();
                if r.is_reverse() {
                    Ok(cigar.leading_softclips())
                } else {
                    Ok(cigar.trailing_softclips())
                }
            })
        });
        fields.add_field_method_get("soft_clips_5_prime", |_, this| {
            this.with_record(|r| {
                let cigar = r.cigar();
                if r.is_reverse() {
                    Ok(cigar.trailing_softclips())
                } else {
                    Ok(cigar.leading_softclips())
                }
            })
        });
        fields.add_field_method_get("average_base_quality", |_, this| {
            this.with_record(|r| {
                let qual = r.qual();
                let sum = qual.iter().map(|q| *q as u64).sum::<u64>();
                let count = qual.len();
                Ok(sum as f64 / count as f64)
            })
        });
    }

    fn add_methods<M: UserDataMethods<Self>>(methods: &mut M) {
        methods.add_method("qpos", |_, this, ()| Ok(this.0.qpos.get()));

        methods.add_method("n_proportion_3_prime", |_, this, n_bases: usize| {
            this.with_record(|r| {
                let seq = r.seq();
                let mut count = 0;
                let reverse = r.is_reverse();
                for i in 0..n_bases {
                    let base =
                        seq[if reverse { i } else { seq.len() - 1 - i }].to_ascii_uppercase();
                    if base == b'N' {
                        count += 1;
                    }
                }
                Ok(count as f64 / n_bases as f64)
            })
        });

        methods.add_method("n_proportion_5_prime", |_, this, n_bases: usize| {
            this.with_record(|r| {
                let seq = r.seq();
                let mut count = 0;
                let reverse = r.is_reverse();
                for i in 0..n_bases {
                    let base =
                        seq[if reverse { seq.len() - 1 - i } else { i }].to_ascii_uppercase();
                    if base == b'N' {
                        count += 1;
                    }
                }
                Ok(count as f64 / n_bases as f64)
            })
        });

        methods.add_method("tag", |lua, this, tag: String| {
            this.with_record(|r| {
                let tag = tag.as_bytes();
                let aux = r.aux(tag).map_err(LuaError::external)?;
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
        });
    }
}

struct LuaReadFilter {
    filter_func: Function,
    read_slot: Arc<ReadSlot>,
    read_ud: AnyUserData,
}

impl LuaReadFilter {
    // Create a new LuaReadFilter instance with the given expression
    fn new(expression: &str, lua: &Lua) -> Result<Self> {
        lua.set_compiler(Compiler::new().set_optimization_level(2));
        let wrapped = format!("return function(read)\n{}\nend", expression);
        let filter_func: Function = lua.load(&wrapped).eval()?;
        let read_slot = Arc::new(ReadSlot {
            record: Cell::new(std::ptr::null()),
            qpos: Cell::new(usize::MAX),
        });
        let read_ud = lua.create_userdata(LuaRead(Arc::clone(&read_slot)))?;
        Ok(Self {
            filter_func,
            read_slot,
            read_ud,
        })
    }
}

struct PileSlot {
    pile: Cell<*const PileupPosition>,
}

// SAFETY: PileSlot is strictly confined to the single thread owning its Lua VM.
// The raw pointer is only non-null during a synchronous pile expression call and
// is always cleared by PileSlotGuard on return or unwind.
unsafe impl Send for PileSlot {}
unsafe impl Sync for PileSlot {}

struct PileSlotGuard<'a>(&'a PileSlot);

impl Drop for PileSlotGuard<'_> {
    #[inline(always)]
    fn drop(&mut self) {
        self.0.pile.set(std::ptr::null());
    }
}

#[derive(Clone)]
struct LuaPile(Arc<PileSlot>);

impl LuaPile {
    #[inline(always)]
    fn with_pile<T>(&self, f: impl FnOnce(&PileupPosition) -> mlua::Result<T>) -> mlua::Result<T> {
        let ptr = self.0.pile.get();
        if ptr.is_null() {
            return Err(LuaError::RuntimeError(
                "pile reference is no longer valid".to_string(),
            ));
        }
        // SAFETY: ptr points to a valid &PileupPosition for the duration of the call.
        f(unsafe { &*ptr })
    }
}

impl UserData for LuaPile {
    fn add_fields<F: UserDataFields<Self>>(fields: &mut F) {
        fields.add_field_method_get("depth", |_, this| this.with_pile(|p| Ok(p.depth)));
        fields.add_field_method_get("a", |_, this| this.with_pile(|p| Ok(p.a)));
        fields.add_field_method_get("c", |_, this| this.with_pile(|p| Ok(p.c)));
        fields.add_field_method_get("g", |_, this| this.with_pile(|p| Ok(p.g)));
        fields.add_field_method_get("t", |_, this| this.with_pile(|p| Ok(p.t)));
        fields.add_field_method_get("n", |_, this| this.with_pile(|p| Ok(p.n)));
        fields.add_field_method_get("fail", |_, this| this.with_pile(|p| Ok(p.fail)));
        fields.add_field_method_get("ins", |_, this| this.with_pile(|p| Ok(p.ins)));
        fields.add_field_method_get("del", |_, this| this.with_pile(|p| Ok(p.del)));
        fields.add_field_method_get("ref_skip", |_, this| this.with_pile(|p| Ok(p.ref_skip)));
        fields.add_field_method_get("pos", |_, this| this.with_pile(|p| Ok(p.pos)));
    }
}

struct LuaPileFilter {
    filter_func: Function,
    pile_slot: Arc<PileSlot>,
    pile_ud: AnyUserData,
}

impl LuaPileFilter {
    fn new(expression: &str, lua: &Lua) -> Result<Self> {
        lua.set_compiler(Compiler::new().set_optimization_level(2));
        let wrapped = format!("return function(pile)\n{}\nend", expression);
        let filter_func: Function = lua.load(&wrapped).eval()?;
        let pile_slot = Arc::new(PileSlot {
            pile: Cell::new(std::ptr::null()),
        });
        let pile_ud = lua.create_userdata(LuaPile(Arc::clone(&pile_slot)))?;
        Ok(Self {
            filter_func,
            pile_slot,
            pile_ud,
        })
    }

    #[inline]
    fn eval(&self, pile: &PileupPosition) -> mlua::Result<bool> {
        self.pile_slot.pile.set(pile as *const PileupPosition);
        let _guard = PileSlotGuard(&self.pile_slot);
        self.filter_func.call::<bool>(&self.pile_ud)
    }
}

impl ReadFilter for LuaReadFilter {
    /// Filter reads based user expression.
    #[inline]
    fn filter_read(&self, read: &Record, alignment: Option<&Alignment>) -> bool {
        self.read_slot.record.set(read as *const Record);
        self.read_slot
            .qpos
            .set(alignment.and_then(|a| a.qpos()).unwrap_or(usize::MAX));
        let _guard = ReadSlotGuard(&self.read_slot);

        match self.filter_func.call::<bool>(&self.read_ud) {
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

        let string_count = lua
            .create_function(|_, (haystack, needle): (String, String)| {
                assert!(needle.len() == 1);
                let needle = needle.chars().next().unwrap();
                Ok(haystack.chars().filter(|c| *c == needle).count())
            })
            .expect("eror creating function");
        lua.globals()
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
#[command(author, version, about, long_about = None, arg_required_else_help = true)]
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
    println!("# pbr version {}", env!("CARGO_PKG_VERSION"));

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
    let pile_filter = if let Some(expression) = opts.pile_expression {
        Some(LuaPileFilter::new(expression.as_str(), &pile_lua)?)
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
            if let Some(pile_filter) = &pile_filter {
                match pile_filter.eval(p) {
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
    fn test_pileup_position() -> Result<()> {
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
        for (expected, expression) in [
            (false, "pile.g > 3"),
            (true, "pile.a > 0"),
            (false, "pile.a > 10"),
            (false, "pile.ref_skip == 100"),
            (true, "pile.ref_skip == 9"),
        ] {
            let pf = LuaPileFilter::new(&(String::from("return ") + expression), &lua)?;
            let result = pf.eval(&pileup_position)?;
            assert_eq!(result, expected, "expression: {}", expression);
        }
        Ok(())
    }
}
