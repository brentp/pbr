#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

mod cached_faidx;
mod expr;
mod processor;

use anyhow::Result;
use cached_faidx::CachedFaidx;
use clap::Parser;
use expr::{
    call_fast_bool, compile_lua_filter, eval_fast_pile_field, eval_fast_read_field,
    is_read_only_expression, partition_read_only_clauses, try_compile_native_read_expr,
    FastPileField, FastReadField, NativeExpr, FAST_PILE_FIELDS, FAST_READ_FIELDS,
};
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
use std::cell::{Cell, RefCell};
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};
use std::path::PathBuf;
use std::sync::Arc;

/// Fast non-cryptographic multiplicative hasher specialized for aligned `usize` pointer
/// keys in `ReadCache`, avoiding `SipHash` overhead on every pileup alignment lookup.
#[derive(Default)]
struct PtrHasher(usize);

impl Hasher for PtrHasher {
    #[inline(always)]
    fn finish(&self) -> u64 {
        self.0 as u64
    }
    #[inline(always)]
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = self
                .0
                .wrapping_mul(0x517cc1b727220a95)
                .wrapping_add(b as usize);
        }
    }
    #[inline(always)]
    fn write_usize(&mut self, i: usize) {
        self.0 = (i >> 3).wrapping_mul(0x517cc1b727220a95);
    }
}

/// Per-region active-read cache mapping `bam1_t.data` pointer -> `(read_pos, passed)`.
/// Exists so read-level (position-independent) filters evaluate once per read (~150x fewer
/// calls) rather than at every covered base; storing `read_pos` guards against htslib
/// buffer reuse across excluded positions while `is_tail()` eviction keeps size O(depth).
type ReadCache = HashMap<usize, (i64, bool), BuildHasherDefault<PtrHasher>>;

/// Reusable thread-confined `Cell`-backed pointer slot that lets a single persistent Lua
/// userdata (`LuaRead` or `LuaPile`) be rebound to each `&Record`/`&PileupPosition` in O(1),
/// avoiding per-call `lua.scope` allocations, `named_user_value` table lookups, and GC churn.
struct PtrSlot<T> {
    ptr: Cell<*const T>,
    qpos: Cell<usize>,
}

// SAFETY: PtrSlot is strictly confined to the single thread owning its Lua VM;
// the raw pointer is only non-null during a synchronous filter call.
unsafe impl<T> Send for PtrSlot<T> {}
unsafe impl<T> Sync for PtrSlot<T> {}

impl<T> PtrSlot<T> {
    fn new() -> Arc<Self> {
        Arc::new(Self {
            ptr: Cell::new(std::ptr::null()),
            qpos: Cell::new(usize::MAX),
        })
    }

    #[inline(always)]
    fn bind(&self, val: &T, qpos: usize) -> PtrSlotGuard<'_, T> {
        self.ptr.set(val as *const T);
        self.qpos.set(qpos);
        PtrSlotGuard(self)
    }

    #[inline(always)]
    fn with_ref<R>(&self, f: impl FnOnce(&T) -> mlua::Result<R>) -> mlua::Result<R> {
        let ptr = self.ptr.get();
        if ptr.is_null() {
            return Err(LuaError::RuntimeError(
                "userdata reference is no longer valid".to_string(),
            ));
        }
        // SAFETY: ptr points to a live &T for the duration of the bound guard.
        f(unsafe { &*ptr })
    }
}

/// RAII guard returned by `PtrSlot::bind` that clears the raw pointer to null on drop
/// (including unwinds), preventing dangling-pointer access if a Lua script stashes `read`
/// or `pile` past the lifetime of a single filter call.
struct PtrSlotGuard<'a, T>(&'a PtrSlot<T>);

impl<T> Drop for PtrSlotGuard<'_, T> {
    #[inline(always)]
    fn drop(&mut self) {
        self.0.ptr.set(std::ptr::null());
        self.0.qpos.set(usize::MAX);
    }
}

/// Persistent `mlua::UserData` handle exposing BAM `Record` fields and methods (`read.*`,
/// `read:tag(...)`, etc.) to Lua through a `PtrSlot<Record>` whenever an expression accesses
/// attributes not handled by direct fast-field parameters or `NativeExpr`.
#[derive(Clone)]
struct LuaRead(Arc<PtrSlot<Record>>);

impl UserData for LuaRead {
    fn add_fields<F: UserDataFields<Self>>(fields: &mut F) {
        fields.add_field_method_get("mapping_quality", |_, this| {
            this.0.with_ref(|r| Ok(r.mapq()))
        });
        fields.add_field_method_get("flags", |_, this| this.0.with_ref(|r| Ok(r.flags())));
        fields.add_field_method_get("tid", |_, this| this.0.with_ref(|r| Ok(r.tid())));
        fields.add_field_method_get("start", |_, this| this.0.with_ref(|r| Ok(r.pos())));
        fields.add_field_method_get("stop", |_, this| {
            this.0.with_ref(|r| Ok(r.cigar().end_pos()))
        });
        fields.add_field_method_get("length", |_, this| this.0.with_ref(|r| Ok(r.seq_len())));
        fields.add_field_method_get("insert_size", |_, this| {
            this.0.with_ref(|r| Ok(r.insert_size()))
        });
        fields.add_field_method_get("qname", |_, this| {
            this.0.with_ref(|r| {
                let q = r.qname();
                Ok(std::str::from_utf8(q).unwrap_or("").to_string())
            })
        });
        fields.add_field_method_get("sequence", |_, this| {
            this.0.with_ref(|r| {
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
            this.0.with_ref(|r| Ok(r.qual()[qpos] as i32))
        });
        fields.add_field_method_get("distance_from_5prime", |_, this| {
            let qpos = this.0.qpos.get();
            if qpos == usize::MAX {
                return Ok(-1);
            }
            this.0.with_ref(|r| {
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
            this.0.with_ref(|r| {
                if r.is_reverse() {
                    Ok(qpos)
                } else {
                    Ok(r.seq_len() - qpos)
                }
            })
        });
        fields.add_field_method_get("indel_count", |_, this| {
            this.0.with_ref(|r| {
                let count = r
                    .cigar()
                    .iter()
                    .filter(|op| matches!(op, Cigar::Ins(_) | Cigar::Del(_)))
                    .count();
                Ok(count)
            })
        });
        fields.add_field_method_get("soft_clips_3_prime", |_, this| {
            this.0.with_ref(|r| {
                let cigar = r.cigar();
                if r.is_reverse() {
                    Ok(cigar.leading_softclips())
                } else {
                    Ok(cigar.trailing_softclips())
                }
            })
        });
        fields.add_field_method_get("soft_clips_5_prime", |_, this| {
            this.0.with_ref(|r| {
                let cigar = r.cigar();
                if r.is_reverse() {
                    Ok(cigar.trailing_softclips())
                } else {
                    Ok(cigar.leading_softclips())
                }
            })
        });
        fields.add_field_method_get("average_base_quality", |_, this| {
            this.0.with_ref(|r| {
                let qual = r.qual();
                let sum = qual.iter().map(|q| *q as u64).sum::<u64>();
                Ok(sum as f64 / qual.len() as f64)
            })
        });
    }

    fn add_methods<M: UserDataMethods<Self>>(methods: &mut M) {
        methods.add_method("qpos", |_, this, ()| Ok(this.0.qpos.get()));

        methods.add_method("n_proportion_3_prime", |_, this, n_bases: usize| {
            this.0.with_ref(|r| {
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
            this.0.with_ref(|r| {
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
            this.0.with_ref(|r| {
                let aux = r.aux(tag.as_bytes()).map_err(LuaError::external)?;
                let lua_val: Value =
                    match aux {
                        Aux::Char(v) => Value::String(lua.create_string(&[v])?),
                        Aux::I8(v) => Value::Number(v as f64),
                        Aux::U8(v) => Value::Number(v as f64),
                        Aux::I16(v) => Value::Number(v as f64),
                        Aux::U16(v) => Value::Number(v as f64),
                        Aux::I32(v) => Value::Number(v as f64),
                        Aux::U32(v) => Value::Number(v as f64),
                        Aux::Float(v) => Value::Number(v as f64),
                        Aux::Double(v) => Value::Number(v as f64),
                        Aux::String(v) => Value::String(lua.create_string(v)?),
                        Aux::ArrayFloat(v) => Value::Table(lua.create_sequence_from(
                            (0..v.len()).map(|i| v.get(i).unwrap_or(f32::NAN)),
                        )?),
                        Aux::ArrayI32(v) => Value::Table(lua.create_sequence_from(
                            (0..v.len()).map(|i| v.get(i).unwrap_or(i32::MIN)),
                        )?),
                        Aux::ArrayI8(v) => Value::Table(lua.create_sequence_from(
                            (0..v.len()).map(|i| v.get(i).unwrap_or(i8::MIN)),
                        )?),
                        Aux::ArrayU8(v) => Value::Table(lua.create_sequence_from(
                            (0..v.len()).map(|i| v.get(i).unwrap_or(u8::MIN)),
                        )?),
                        Aux::ArrayU16(v) => Value::Table(lua.create_sequence_from(
                            (0..v.len()).map(|i| v.get(i).unwrap_or(u16::MIN)),
                        )?),
                        Aux::ArrayU32(v) => Value::Table(lua.create_sequence_from(
                            (0..v.len()).map(|i| v.get(i).unwrap_or(u32::MIN)),
                        )?),
                        Aux::ArrayI16(v) => Value::Table(lua.create_sequence_from(
                            (0..v.len()).map(|i| v.get(i).unwrap_or(i16::MIN)),
                        )?),
                        Aux::HexByteArray(v) => {
                            let lstr = String::from_utf8_lossy(v.as_bytes());
                            Value::String(lua.create_string(lstr.as_ref())?)
                        }
                    };
                Ok(Some(lua_val))
            })
        });
    }
}

/// A compiled read-filter expression backed by either a pure-Rust `NativeExpr` or a
/// compiled Luau `Function` with extracted `FastReadField` parameters. Exists so callers
/// transparently get native Rust speed for common predicates and full Luau generality for
/// everything else.
enum CompiledReadExpr {
    Native(NativeExpr),
    Lua {
        func: Function,
        fast_fields: Vec<FastReadField>,
        needs_userdata: bool,
    },
}

impl CompiledReadExpr {
    fn compile(expression: &str, lua: &Lua) -> Result<Self> {
        if let Some(native) = try_compile_native_read_expr(expression) {
            return Ok(CompiledReadExpr::Native(native));
        }
        let (func, fast_fields, needs_userdata) =
            compile_lua_filter(expression, "read", FAST_READ_FIELDS, lua)?;
        Ok(CompiledReadExpr::Lua {
            func,
            fast_fields,
            needs_userdata,
        })
    }

    #[inline]
    fn eval(
        &self,
        read: &Record,
        qpos: usize,
        read_slot: &PtrSlot<Record>,
        read_ud: &AnyUserData,
    ) -> bool {
        match self {
            CompiledReadExpr::Native(expr) => expr.eval(read, qpos),
            CompiledReadExpr::Lua {
                func,
                fast_fields,
                needs_userdata,
            } => {
                let mut vals = [0.0f64; 9];
                for (idx, &field) in fast_fields.iter().enumerate() {
                    vals[idx] = eval_fast_read_field(field, read, qpos);
                }
                let res = if *needs_userdata {
                    let _guard = read_slot.bind(read, qpos);
                    call_fast_bool(func, Some(read_ud), &vals, fast_fields.len())
                } else {
                    call_fast_bool(func, None, &vals, fast_fields.len())
                };
                match res {
                    Ok(r) => r,
                    Err(e) => {
                        eprintln!("Error evaluating expression: {}", e);
                        false
                    }
                }
            }
        }
    }
}

/// `perbase::ReadFilter` implementation that evaluates the user's Lua read expression
/// using three fast paths: whole-expression native Rust evaluation, per-read caching of
/// read-only expressions (or hoisted top-level `and`-clauses), and per-base evaluation
/// with direct stack-passed arguments and a persistent `LuaRead` userdata.
struct LuaReadFilter {
    read_only_filter: Option<CompiledReadExpr>,
    base_filter: Option<CompiledReadExpr>,
    read_slot: Arc<PtrSlot<Record>>,
    read_ud: AnyUserData,
    read_cache: RefCell<ReadCache>,
}

impl LuaReadFilter {
    fn new(expression: &str, lua: &Lua) -> Result<Self> {
        lua.set_compiler(Compiler::new().set_optimization_level(2));

        // If the entire expression compiles to native Rust, evaluate directly with zero
        // Lua calls and zero HashMap cache overhead.
        let (read_only_filter, base_filter) =
            if let Some(native) = try_compile_native_read_expr(expression) {
                (None, Some(CompiledReadExpr::Native(native)))
            } else if is_read_only_expression(expression) {
                (Some(CompiledReadExpr::compile(expression, lua)?), None)
            } else if let Some((read_only, per_base)) = partition_read_only_clauses(expression) {
                let compile_and = |clauses: Vec<&str>| {
                    let joined = clauses
                        .iter()
                        .map(|c| format!("({})", c))
                        .collect::<Vec<_>>()
                        .join(" and ");
                    CompiledReadExpr::compile(&format!("return {}", joined), lua)
                };
                (Some(compile_and(read_only)?), Some(compile_and(per_base)?))
            } else {
                (None, Some(CompiledReadExpr::compile(expression, lua)?))
            };

        let read_slot = PtrSlot::new();
        let read_ud = lua.create_userdata(LuaRead(Arc::clone(&read_slot)))?;
        Ok(Self {
            read_only_filter,
            base_filter,
            read_slot,
            read_ud,
            read_cache: RefCell::new(ReadCache::default()),
        })
    }

    #[inline]
    fn eval_read_only(
        &self,
        expr: &CompiledReadExpr,
        read: &Record,
        alignment: Option<&Alignment>,
    ) -> bool {
        let Some(aln) = alignment else {
            return expr.eval(read, usize::MAX, &self.read_slot, &self.read_ud);
        };
        let key = read.inner().data as usize;
        let pos = read.pos();
        let is_tail = aln.is_tail();
        if !aln.is_head() {
            let mut cache = self.read_cache.borrow_mut();
            if let Some(&(cached_pos, cached_res)) = cache.get(&key) {
                if cached_pos == pos {
                    if is_tail {
                        cache.remove(&key);
                    }
                    return cached_res;
                }
            }
        }
        let res = expr.eval(read, usize::MAX, &self.read_slot, &self.read_ud);
        let mut cache = self.read_cache.borrow_mut();
        if is_tail {
            cache.remove(&key);
        } else {
            cache.insert(key, (pos, res));
        }
        res
    }
}

impl ReadFilter for LuaReadFilter {
    #[inline]
    fn filter_read(&self, read: &Record, alignment: Option<&Alignment>) -> bool {
        if let Some(ro) = &self.read_only_filter {
            if !self.eval_read_only(ro, read, alignment) {
                return false;
            }
        }
        if let Some(base) = &self.base_filter {
            let qpos = alignment.and_then(|a| a.qpos()).unwrap_or(usize::MAX);
            base.eval(read, qpos, &self.read_slot, &self.read_ud)
        } else {
            true
        }
    }
}

/// Persistent `mlua::UserData` handle exposing `PileupPosition` fields (`pile.*`) to Lua
/// through a `PtrSlot<PileupPosition>`, used as a fallback if `--pile-expression` references
/// `pile` in a way that cannot be rewritten to direct primitive parameters.
#[derive(Clone)]
struct LuaPile(Arc<PtrSlot<PileupPosition>>);

impl UserData for LuaPile {
    fn add_fields<F: UserDataFields<Self>>(fields: &mut F) {
        fields.add_field_method_get("depth", |_, this| this.0.with_ref(|p| Ok(p.depth)));
        fields.add_field_method_get("a", |_, this| this.0.with_ref(|p| Ok(p.a)));
        fields.add_field_method_get("c", |_, this| this.0.with_ref(|p| Ok(p.c)));
        fields.add_field_method_get("g", |_, this| this.0.with_ref(|p| Ok(p.g)));
        fields.add_field_method_get("t", |_, this| this.0.with_ref(|p| Ok(p.t)));
        fields.add_field_method_get("n", |_, this| this.0.with_ref(|p| Ok(p.n)));
        fields.add_field_method_get("fail", |_, this| this.0.with_ref(|p| Ok(p.fail)));
        fields.add_field_method_get("ins", |_, this| this.0.with_ref(|p| Ok(p.ins)));
        fields.add_field_method_get("del", |_, this| this.0.with_ref(|p| Ok(p.del)));
        fields.add_field_method_get("ref_skip", |_, this| this.0.with_ref(|p| Ok(p.ref_skip)));
        fields.add_field_method_get("pos", |_, this| this.0.with_ref(|p| Ok(p.pos)));
    }
}

/// Compiled `--pile-expression` filter that rewrites `pile.<field>` accesses to stack-passed
/// `f64` parameters and reuses a single persistent `LuaPile` userdata across pileup columns.
struct LuaPileFilter {
    filter_func: Function,
    pile_slot: Arc<PtrSlot<PileupPosition>>,
    pile_ud: AnyUserData,
    fast_fields: Vec<FastPileField>,
    needs_userdata: bool,
}

impl LuaPileFilter {
    fn new(expression: &str, lua: &Lua) -> Result<Self> {
        lua.set_compiler(Compiler::new().set_optimization_level(2));
        let (filter_func, fast_fields, needs_userdata) =
            compile_lua_filter(expression, "pile", FAST_PILE_FIELDS, lua)?;
        let pile_slot = PtrSlot::new();
        let pile_ud = lua.create_userdata(LuaPile(Arc::clone(&pile_slot)))?;
        Ok(Self {
            filter_func,
            pile_slot,
            pile_ud,
            fast_fields,
            needs_userdata,
        })
    }

    #[inline]
    fn eval(&self, pile: &PileupPosition) -> mlua::Result<bool> {
        let mut vals = [0.0f64; 11];
        for (idx, &field) in self.fast_fields.iter().enumerate() {
            vals[idx] = eval_fast_pile_field(field, pile);
        }
        if self.needs_userdata {
            let _guard = self.pile_slot.bind(pile, usize::MAX);
            call_fast_bool(
                &self.filter_func,
                Some(&self.pile_ud),
                &vals,
                self.fast_fields.len(),
            )
        } else {
            call_fast_bool(&self.filter_func, None, &vals, self.fast_fields.len())
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

        let rf = LuaReadFilter::new(&self.expression, &lua).unwrap_or_else(|_| {
            panic!(
                "error creating lua read filter with expression {}",
                &self.expression
            )
        });

        let exclude_intervals = self.exclude_regions.as_ref().map(|regions_bed| {
            Self::bed_to_intervals(&header, regions_bed, true).expect("BED file")
        });

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

/// Command-line arguments parsed by `clap` specifying the input BAM/CRAM, Lua read and
/// optional pileup filter expressions, thread count, max depth, reference FASTA, BED
/// include/exclude regions, and mate-overlap deduplication.
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
        let mut header = Header::new();
        let mut sq = HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", "chr1");
        sq.push_tag(b"LN", &1000000u32);
        header.push_record(&sq);
        let header_view = HeaderView::from_header(&header);

        let record = Record::from_sam(
            &header_view,
            b"test_read\t0\tchr1\t100\t30\t4M\t*\t0\t0\tACGT\t&&&&\tRG:Z:test",
        )
        .expect("Failed to create record from SAM");

        let tmp = NamedTempFile::new()?;
        let path = tmp.path();
        {
            let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam)?;
            writer.write(&record)?;
        }
        bam::index::build(path, None, bam::index::Type::Bai, 1)?;

        let mut reader = IndexedReader::from_path(path)?;
        reader.fetch(("chr1", 100, 101))?;
        let mut pileups = reader.pileup();
        let pileup: Pileup = pileups.next().unwrap().expect("Failed to get pileup");
        let alignment = pileup
            .alignments()
            .next()
            .expect("No alignment found in pileup");

        let lua = Lua::new();
        let rf = LuaReadFilter::new(
            "return read.bq > 0 and read.distance_from_5prime == 0 and read.distance_from_3prime > 0",
            &lua,
        )?;

        let result = rf.filter_read(&alignment.record(), Some(&alignment));
        assert!(result);

        // Split into a cached read-level clause and a per-base clause; the trailing
        // comments must not break compilation of the re-joined clauses.
        let rf = LuaReadFilter::new(
            "return read:tag('RG') == 'test' -- rg\n and read.bq > 0 -- bq\n",
            &lua,
        )?;
        assert!(rf.read_only_filter.is_some() && rf.base_filter.is_some());
        assert!(rf.filter_read(&alignment.record(), Some(&alignment)));

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
