//! Expression analysis, fast-field argument rewriting, read-level clause splitting,
//! and native Rust fast-path compilation for `pbr` Lua filter expressions.
//!
//! Expressions are parsed with full_moon's Luau parser; anything the analyses don't
//! understand is handed to Luau unchanged.

use anyhow::Result;
use full_moon::ast::{
    Assignment, Ast, BinOp, Call, CompoundAssignment, Expression, FunctionArgs, FunctionCall,
    Index, LastStmt, Prefix, Suffix, UnOp, Var, VarExpression,
};
use full_moon::node::Node;
use full_moon::tokenizer::{Symbol, TokenReference, TokenType};
use full_moon::visitors::{Visit, Visitor};
use mlua::{AnyUserData, Function, Lua};
use perbase_lib::position::pileup_position::PileupPosition;
use rust_htslib::bam::record::Record;

/// O(1) numeric fields on `read` that can be extracted in Rust and passed directly as
/// stack arguments to Lua (or evaluated by `NativeExpr`), avoiding per-field Lua
/// `__index` userdata dispatch overhead on every read-base.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) enum FastReadField {
    MappingQuality,
    Flags,
    Tid,
    Start,
    Length,
    InsertSize,
    Bq,
    DistanceFrom5Prime,
    DistanceFrom3Prime,
}

pub(crate) const FAST_READ_FIELDS: &[(&str, FastReadField)] = &[
    ("mapping_quality", FastReadField::MappingQuality),
    ("flags", FastReadField::Flags),
    ("tid", FastReadField::Tid),
    ("start", FastReadField::Start),
    ("length", FastReadField::Length),
    ("insert_size", FastReadField::InsertSize),
    ("bq", FastReadField::Bq),
    ("distance_from_5prime", FastReadField::DistanceFrom5Prime),
    ("distance_from_3prime", FastReadField::DistanceFrom3Prime),
];

/// O(1) numeric fields on `pile` that can be extracted in Rust and passed directly as
/// stack arguments to `--pile-expression`, avoiding Lua `__index` userdata lookup overhead.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub(crate) enum FastPileField {
    Depth,
    A,
    C,
    G,
    T,
    N,
    Fail,
    Ins,
    Del,
    RefSkip,
    Pos,
}

pub(crate) const FAST_PILE_FIELDS: &[(&str, FastPileField)] = &[
    ("depth", FastPileField::Depth),
    ("a", FastPileField::A),
    ("c", FastPileField::C),
    ("g", FastPileField::G),
    ("t", FastPileField::T),
    ("n", FastPileField::N),
    ("fail", FastPileField::Fail),
    ("ins", FastPileField::Ins),
    ("del", FastPileField::Del),
    ("ref_skip", FastPileField::RefSkip),
    ("pos", FastPileField::Pos),
];

const READ_ONLY_FIELDS: &[&str] = &[
    "mapping_quality",
    "flags",
    "tid",
    "start",
    "stop",
    "length",
    "insert_size",
    "qname",
    "sequence",
    "indel_count",
    "soft_clips_3_prime",
    "soft_clips_5_prime",
    "average_base_quality",
];

const READ_ONLY_METHODS: &[&str] = &["n_proportion_3_prime", "n_proportion_5_prime", "tag"];

#[inline(always)]
pub(crate) fn eval_fast_read_field(field: FastReadField, read: &Record, qpos: usize) -> f64 {
    match field {
        FastReadField::MappingQuality => read.mapq() as f64,
        FastReadField::Flags => read.flags() as f64,
        FastReadField::Tid => read.tid() as f64,
        FastReadField::Start => read.pos() as f64,
        FastReadField::Length => read.seq_len() as f64,
        FastReadField::InsertSize => read.insert_size() as f64,
        FastReadField::Bq => {
            if qpos == usize::MAX {
                -1.0
            } else {
                read.qual()[qpos] as f64
            }
        }
        FastReadField::DistanceFrom5Prime => {
            if qpos == usize::MAX {
                -1.0
            } else if read.is_reverse() {
                (read.seq_len() as i32 - qpos as i32) as f64
            } else {
                qpos as f64
            }
        }
        FastReadField::DistanceFrom3Prime => {
            if qpos == usize::MAX {
                usize::MAX as f64
            } else if read.is_reverse() {
                qpos as f64
            } else {
                (read.seq_len() - qpos) as f64
            }
        }
    }
}

#[inline(always)]
pub(crate) fn eval_fast_pile_field(field: FastPileField, pile: &PileupPosition) -> f64 {
    match field {
        FastPileField::Depth => pile.depth as f64,
        FastPileField::A => pile.a as f64,
        FastPileField::C => pile.c as f64,
        FastPileField::G => pile.g as f64,
        FastPileField::T => pile.t as f64,
        FastPileField::N => pile.n as f64,
        FastPileField::Fail => pile.fail as f64,
        FastPileField::Ins => pile.ins as f64,
        FastPileField::Del => pile.del as f64,
        FastPileField::RefSkip => pile.ref_skip as f64,
        FastPileField::Pos => pile.pos as f64,
    }
}

// ---------------------------------------------------------------------------
// Expression analysis on the full_moon (Luau) AST
// ---------------------------------------------------------------------------

/// Parse `expr` as a Luau chunk. Syntax errors return `None`, so callers fall back to
/// handing the original text to Luau, which reports the error.
fn parse(expr: &str) -> Option<Ast> {
    full_moon::parse(expr).ok()
}

/// The expression of a chunk that is exactly `return <expr>`.
fn single_return(ast: &Ast) -> Option<&Expression> {
    let block = ast.nodes();
    if block.stmts().next().is_some() {
        return None;
    }
    let LastStmt::Return(ret) = block.last_stmt()? else {
        return None;
    };
    let mut returns = ret.returns().iter();
    let expr = returns.next()?;
    returns.next().is_none().then_some(expr)
}

fn ident(t: &TokenReference) -> Option<&str> {
    match t.token_type() {
        TokenType::Identifier { identifier } => Some(identifier.as_str()),
        _ => None,
    }
}

/// `[start, end)` byte span of a node in the source, excluding surrounding comments and
/// whitespace.
fn byte_span(node: &impl Node) -> Option<(usize, usize)> {
    let (start, end) = node.range()?;
    Some((start.bytes(), end.bytes()))
}

/// How an expression uses one object name (`read` or `pile`).
struct ObjUses<'a> {
    obj: &'a str,
    /// `obj.<field>` accesses: field name and the byte span of `obj.<field>`.
    fields: Vec<(String, usize, usize)>,
    /// `obj:<method>(...)` calls.
    methods: Vec<String>,
    /// `obj` used any other way: bare, `obj[...]`, `obj(...)`, `(obj).x`.
    other: bool,
    /// `obj` or one of its fields is assigned to.
    assigned: bool,
    /// `print` or `io` is referenced, so evaluation has side effects.
    side_effects: bool,
}

impl<'a> ObjUses<'a> {
    fn of(node: &impl Visit, obj: &'a str) -> Self {
        let mut uses = ObjUses {
            obj,
            fields: Vec::new(),
            methods: Vec::new(),
            other: false,
            assigned: false,
            side_effects: false,
        };
        node.visit(&mut uses);
        uses
    }

    fn name(&mut self, t: &TokenReference) -> bool {
        match ident(t) {
            Some("print" | "io") => {
                self.side_effects = true;
                false
            }
            Some(n) => n == self.obj,
            None => false,
        }
    }

    fn prefixed<'s>(&mut self, prefix: &Prefix, mut suffixes: impl Iterator<Item = &'s Suffix>) {
        let Prefix::Name(name) = prefix else {
            return;
        };
        if !self.name(name) {
            return;
        }
        match suffixes.next() {
            Some(Suffix::Index(Index::Dot { name: field, .. })) => {
                match (ident(field), byte_span(name), byte_span(field)) {
                    (Some(f), Some((start, _)), Some((_, end))) => {
                        self.fields.push((f.to_string(), start, end))
                    }
                    _ => self.other = true,
                }
            }
            Some(Suffix::Call(Call::MethodCall(m))) => match ident(m.name()) {
                Some(m) => self.methods.push(m.to_string()),
                None => self.other = true,
            },
            _ => self.other = true,
        }
    }

    fn assign_target(&mut self, var: &Var) {
        let target = match var {
            Var::Name(t) => Some(t),
            Var::Expression(v) => match v.prefix() {
                Prefix::Name(t) => Some(t),
                _ => None,
            },
            _ => None,
        };
        if target.and_then(ident) == Some(self.obj) {
            self.assigned = true;
        }
    }

    /// True iff every use of `read` is a known position-independent field or method and
    /// there are no side effects.
    fn read_only(&self) -> bool {
        !self.other
            && !self.assigned
            && !self.side_effects
            && self
                .fields
                .iter()
                .all(|(f, ..)| READ_ONLY_FIELDS.contains(&f.as_str()))
            && self
                .methods
                .iter()
                .all(|m| READ_ONLY_METHODS.contains(&m.as_str()))
    }
}

impl Visitor for ObjUses<'_> {
    fn visit_var(&mut self, var: &Var) {
        if let Var::Name(t) = var {
            if self.name(t) {
                self.other = true;
            }
        }
    }

    fn visit_var_expression(&mut self, v: &VarExpression) {
        self.prefixed(v.prefix(), v.suffixes());
    }

    fn visit_function_call(&mut self, c: &FunctionCall) {
        self.prefixed(c.prefix(), c.suffixes());
    }

    fn visit_assignment(&mut self, a: &Assignment) {
        for var in a.variables() {
            self.assign_target(var);
        }
    }

    fn visit_compound_assignment(&mut self, a: &CompoundAssignment) {
        self.assign_target(a.lhs());
    }
}

// ---------------------------------------------------------------------------
// Fast-field parameter rewriting + zero-alloc tuple invocation
// ---------------------------------------------------------------------------

/// Rewrite direct `obj.<fast_field>` accesses to local parameters `_pbr_<field>` so
/// Luau operates on stack-passed `f64`s instead of `__index` userdata callbacks.
/// Returns the rewritten source, the fields used (in parameter order), and whether `obj`
/// itself is still needed.
pub(crate) fn rewrite_fast_fields<F: Copy + PartialEq>(
    expr: &str,
    obj_name: &str,
    table: &[(&'static str, F)],
) -> (String, Vec<(&'static str, F)>, bool) {
    let Some(ast) = parse(expr) else {
        return (expr.to_string(), Vec::new(), true);
    };
    let uses = ObjUses::of(ast.nodes(), obj_name);
    if uses.assigned {
        return (expr.to_string(), Vec::new(), true);
    }

    let mut fields = uses.fields;
    fields.sort_by_key(|&(_, start, _)| start);
    let mut out = String::with_capacity(expr.len() + 32);
    let mut used: Vec<(&'static str, F)> = Vec::new();
    let mut needs_userdata = uses.other || !uses.methods.is_empty();
    let mut cursor = 0;
    for (field, start, end) in fields {
        let Some(&(fname, fenum)) = table.iter().find(|(n, _)| *n == field) else {
            needs_userdata = true;
            continue;
        };
        if !used.iter().any(|(_, existing)| *existing == fenum) {
            used.push((fname, fenum));
        }
        out.push_str(&expr[cursor..start]);
        out.push_str("_pbr_");
        out.push_str(fname);
        cursor = end;
    }
    out.push_str(&expr[cursor..]);
    (out, used, needs_userdata)
}

/// Compile a Lua filter expression wrapped as `return function([obj,] _pbr_f1, ...) ... end`.
pub(crate) fn compile_lua_filter<F: Copy + PartialEq>(
    expression: &str,
    obj_name: &str,
    table: &[(&'static str, F)],
    lua: &Lua,
) -> Result<(Function, Vec<F>, bool)> {
    let (rewritten, used_fields, needs_userdata) = rewrite_fast_fields(expression, obj_name, table);
    let mut params = Vec::with_capacity(used_fields.len() + 1);
    if needs_userdata {
        params.push(obj_name.to_string());
    }
    let mut fast_fields = Vec::with_capacity(used_fields.len());
    for (name, field) in used_fields {
        params.push(format!("_pbr_{}", name));
        fast_fields.push(field);
    }
    let wrapped = format!("return function({})\n{}\nend", params.join(", "), rewritten);
    let func: Function = lua.load(&wrapped).eval()?;
    Ok((func, fast_fields, needs_userdata))
}

#[inline(always)]
pub(crate) fn call_fast_bool(
    func: &Function,
    ud: Option<&AnyUserData>,
    v: &[f64],
    n: usize,
) -> mlua::Result<bool> {
    macro_rules! dispatch {
        ($($len:literal => ($($idx:tt),*)),* $(,)?) => {
            match (ud, n) {
                (None, 0) => func.call(()),
                (Some(u), 0) => func.call(u),
                $(
                    (None, $len) => func.call(($(v[$idx],)*)),
                    (Some(u), $len) => func.call((u, $(v[$idx],)*)),
                )*
                _ => unreachable!("more fast fields than supported"),
            }
        };
    }
    dispatch!(
        1 => (0),
        2 => (0, 1),
        3 => (0, 1, 2),
        4 => (0, 1, 2, 3),
        5 => (0, 1, 2, 3, 4),
        6 => (0, 1, 2, 3, 4, 5),
        7 => (0, 1, 2, 3, 4, 5, 6),
        8 => (0, 1, 2, 3, 4, 5, 6, 7),
        9 => (0, 1, 2, 3, 4, 5, 6, 7, 8),
        10 => (0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
        11 => (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    )
}

// ---------------------------------------------------------------------------
// Read-only detection & top-level `and`-clause partitioning
// ---------------------------------------------------------------------------

/// True iff every `read` access in `expr` is a known position-independent field/method
/// and the expression does not reference `print` or `io`.
pub(crate) fn is_read_only_expression(expr: &str) -> bool {
    parse(expr).is_some_and(|ast| ObjUses::of(ast.nodes(), "read").read_only())
}

fn flatten_and<'e>(expr: &'e Expression, out: &mut Vec<&'e Expression>) {
    match expr {
        Expression::BinaryOperator {
            lhs,
            binop: BinOp::And(_),
            rhs,
        } => {
            flatten_and(lhs, out);
            flatten_and(rhs, out);
        }
        _ => out.push(expr),
    }
}

/// If `expr` is `return <c1> and <c2> ...` with no side effects, partition its top-level
/// `and` clauses into (read-only, position-dependent) source slices. Returns `None` unless
/// both groups are non-empty.
pub(crate) fn partition_read_only_clauses(expr: &str) -> Option<(Vec<&str>, Vec<&str>)> {
    let ast = parse(expr)?;
    let ret = single_return(&ast)?;
    if ObjUses::of(ret, "read").side_effects {
        return None;
    }
    let mut clauses = Vec::new();
    flatten_and(ret, &mut clauses);
    let (mut read_only, mut per_base) = (Vec::new(), Vec::new());
    for clause in clauses {
        let (start, end) = byte_span(clause)?;
        if ObjUses::of(clause, "read").read_only() {
            read_only.push(&expr[start..end]);
        } else {
            per_base.push(&expr[start..end]);
        }
    }
    (!read_only.is_empty() && !per_base.is_empty()).then_some((read_only, per_base))
}

// ---------------------------------------------------------------------------
// Rust-native fast path for simple boolean / numeric / bit32 read expressions
// ---------------------------------------------------------------------------

/// Numeric expression node in the Rust-native fast-path AST (constants, `FastReadField`
/// loads, arithmetic, and `bit32` ops). Exists to evaluate common numeric read predicates
/// directly in Rust without crossing the Rust→Lua FFI boundary on every base.
#[derive(Debug, Clone, PartialEq)]
enum NativeNum {
    Const(f64),
    Field(FastReadField),
    Neg(Box<NativeNum>),
    Add(Box<NativeNum>, Box<NativeNum>),
    Sub(Box<NativeNum>, Box<NativeNum>),
    Mul(Box<NativeNum>, Box<NativeNum>),
    Div(Box<NativeNum>, Box<NativeNum>),
    BitAnd(Vec<NativeNum>),
    BitOr(Vec<NativeNum>),
    BitXor(Vec<NativeNum>),
    BitNot(Box<NativeNum>),
}

impl NativeNum {
    #[inline]
    fn eval(&self, read: &Record, qpos: usize) -> f64 {
        match self {
            NativeNum::Const(v) => *v,
            NativeNum::Field(f) => eval_fast_read_field(*f, read, qpos),
            NativeNum::Neg(a) => -a.eval(read, qpos),
            NativeNum::Add(a, b) => a.eval(read, qpos) + b.eval(read, qpos),
            NativeNum::Sub(a, b) => a.eval(read, qpos) - b.eval(read, qpos),
            NativeNum::Mul(a, b) => a.eval(read, qpos) * b.eval(read, qpos),
            NativeNum::Div(a, b) => a.eval(read, qpos) / b.eval(read, qpos),
            NativeNum::BitAnd(args) => args
                .iter()
                .fold(u32::MAX, |acc, a| acc & (a.eval(read, qpos) as i64) as u32)
                as f64,
            NativeNum::BitOr(args) => args
                .iter()
                .fold(0u32, |acc, a| acc | (a.eval(read, qpos) as i64) as u32)
                as f64,
            NativeNum::BitXor(args) => args
                .iter()
                .fold(0u32, |acc, a| acc ^ (a.eval(read, qpos) as i64) as u32)
                as f64,
            NativeNum::BitNot(a) => (!((a.eval(read, qpos) as i64) as u32)) as f64,
        }
    }
}

/// Numeric comparison operator used by `NativeBool::Cmp`, evaluated with IEEE-754 `f64`
/// semantics to match Lua comparisons.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CmpOp {
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
}

/// Boolean expression node in the Rust-native fast-path AST supporting short-circuiting
/// `and`/`or`, `not`, and numeric comparisons.
#[derive(Debug, Clone, PartialEq)]
enum NativeBool {
    Const(bool),
    Not(Box<NativeBool>),
    And(Box<NativeBool>, Box<NativeBool>),
    Or(Box<NativeBool>, Box<NativeBool>),
    Cmp(CmpOp, NativeNum, NativeNum),
}

/// Opaque compiled native Rust representation of a simple boolean/numeric/`bit32` read
/// expression. Wraps `NativeBool` so AST internals stay private to `expr` while letting
/// callers run hot per-base filters with zero Lua calls.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct NativeExpr(NativeBool);

impl NativeExpr {
    #[inline]
    pub(crate) fn eval(&self, read: &Record, qpos: usize) -> bool {
        self.0.eval(read, qpos)
    }
}

impl NativeBool {
    #[inline]
    fn eval(&self, read: &Record, qpos: usize) -> bool {
        match self {
            NativeBool::Const(b) => *b,
            NativeBool::Not(a) => !a.eval(read, qpos),
            NativeBool::And(a, b) => a.eval(read, qpos) && b.eval(read, qpos),
            NativeBool::Or(a, b) => a.eval(read, qpos) || b.eval(read, qpos),
            NativeBool::Cmp(op, a, b) => {
                let (av, bv) = (a.eval(read, qpos), b.eval(read, qpos));
                match op {
                    CmpOp::Eq => av == bv,
                    CmpOp::Ne => av != bv,
                    CmpOp::Lt => av < bv,
                    CmpOp::Le => av <= bv,
                    CmpOp::Gt => av > bv,
                    CmpOp::Ge => av >= bv,
                }
            }
        }
    }
}

/// Recursive-descent parser over `&[Token]` that compiles the supported boolean/numeric/
/// `bit32` subset into a `NativeBool` (constant-folding pure-constant `bit32` calls) and
/// Lua number literal → `f64` (decimal, `0x` hex, Luau `0b` binary and `_` separators).
fn parse_number(t: &TokenReference) -> Option<f64> {
    let TokenType::Number { text } = t.token_type() else {
        return None;
    };
    let s = text.as_str().replace('_', "").to_ascii_lowercase();
    if let Some(hex) = s.strip_prefix("0x") {
        u64::from_str_radix(hex, 16).ok().map(|v| v as f64)
    } else if let Some(bin) = s.strip_prefix("0b") {
        u64::from_str_radix(bin, 2).ok().map(|v| v as f64)
    } else {
        s.parse().ok()
    }
}

#[inline]
fn to_u32(v: f64) -> u32 {
    // Matches Luau's bit32 argument conversion: (unsigned)(long long)n.
    (v as i64) as u32
}

/// Compile a boolean-valued Luau expression to `NativeBool`, or `None` if it uses anything
/// outside the supported subset (execution then falls back to Luau). Precedence comes from
/// full_moon, so e.g. `not a == b` is `(not a) == b` and is rejected here.
fn native_bool(e: &Expression) -> Option<NativeBool> {
    Some(match e {
        Expression::Parentheses { expression, .. } => native_bool(expression)?,
        Expression::UnaryOperator {
            unop: UnOp::Not(_),
            expression,
        } => NativeBool::Not(Box::new(native_bool(expression)?)),
        Expression::Symbol(t) => match t.token_type() {
            TokenType::Symbol {
                symbol: Symbol::True,
            } => NativeBool::Const(true),
            TokenType::Symbol {
                symbol: Symbol::False,
            } => NativeBool::Const(false),
            _ => return None,
        },
        Expression::BinaryOperator { lhs, binop, rhs } => {
            let op = match binop {
                BinOp::And(_) => {
                    return Some(NativeBool::And(
                        Box::new(native_bool(lhs)?),
                        Box::new(native_bool(rhs)?),
                    ))
                }
                BinOp::Or(_) => {
                    return Some(NativeBool::Or(
                        Box::new(native_bool(lhs)?),
                        Box::new(native_bool(rhs)?),
                    ))
                }
                BinOp::TwoEqual(_) => CmpOp::Eq,
                BinOp::TildeEqual(_) => CmpOp::Ne,
                BinOp::LessThan(_) => CmpOp::Lt,
                BinOp::LessThanEqual(_) => CmpOp::Le,
                BinOp::GreaterThan(_) => CmpOp::Gt,
                BinOp::GreaterThanEqual(_) => CmpOp::Ge,
                _ => return None,
            };
            NativeBool::Cmp(op, native_num(lhs)?, native_num(rhs)?)
        }
        _ => return None,
    })
}

fn native_num(e: &Expression) -> Option<NativeNum> {
    Some(match e {
        Expression::Number(t) => NativeNum::Const(parse_number(t)?),
        Expression::Parentheses { expression, .. } => native_num(expression)?,
        Expression::UnaryOperator {
            unop: UnOp::Minus(_),
            expression,
        } => match native_num(expression)? {
            NativeNum::Const(v) => NativeNum::Const(-v),
            other => NativeNum::Neg(Box::new(other)),
        },
        Expression::BinaryOperator { lhs, binop, rhs } => {
            let (a, b) = (Box::new(native_num(lhs)?), Box::new(native_num(rhs)?));
            match binop {
                BinOp::Plus(_) => NativeNum::Add(a, b),
                BinOp::Minus(_) => NativeNum::Sub(a, b),
                BinOp::Star(_) => NativeNum::Mul(a, b),
                BinOp::Slash(_) => NativeNum::Div(a, b),
                _ => return None,
            }
        }
        Expression::Var(Var::Expression(v)) => native_read_field(v)?,
        Expression::FunctionCall(c) => native_bit32_call(c)?,
        _ => return None,
    })
}

/// `read.<fast field>` with nothing after it.
fn native_read_field(v: &VarExpression) -> Option<NativeNum> {
    let Prefix::Name(obj) = v.prefix() else {
        return None;
    };
    if ident(obj) != Some("read") {
        return None;
    }
    let mut suffixes = v.suffixes();
    let (Some(Suffix::Index(Index::Dot { name, .. })), None) = (suffixes.next(), suffixes.next())
    else {
        return None;
    };
    let field = ident(name)?;
    FAST_READ_FIELDS
        .iter()
        .find(|(n, _)| *n == field)
        .map(|&(_, f)| NativeNum::Field(f))
}

/// `bit32.band/bor/bxor/bnot(...)`, constant-folded when every argument is a constant.
fn native_bit32_call(c: &FunctionCall) -> Option<NativeNum> {
    let Prefix::Name(module) = c.prefix() else {
        return None;
    };
    if ident(module) != Some("bit32") {
        return None;
    }
    let mut suffixes = c.suffixes();
    let (
        Some(Suffix::Index(Index::Dot { name, .. })),
        Some(Suffix::Call(Call::AnonymousCall(call_args))),
        None,
    ) = (suffixes.next(), suffixes.next(), suffixes.next())
    else {
        return None;
    };
    let FunctionArgs::Parentheses { arguments, .. } = call_args.as_ref() else {
        return None;
    };
    let mut args = arguments
        .iter()
        .map(native_num)
        .collect::<Option<Vec<_>>>()?;
    let consts: Option<Vec<u32>> = args
        .iter()
        .map(|a| match a {
            NativeNum::Const(v) => Some(to_u32(*v)),
            _ => None,
        })
        .collect();

    match (ident(name)?, consts) {
        ("band", Some(vs)) => Some(NativeNum::Const(
            vs.into_iter().fold(u32::MAX, |a, b| a & b) as f64,
        )),
        ("band", None) => Some(NativeNum::BitAnd(args)),
        ("bor", Some(vs)) => Some(NativeNum::Const(
            vs.into_iter().fold(0u32, |a, b| a | b) as f64
        )),
        ("bor", None) => Some(NativeNum::BitOr(args)),
        ("bxor", Some(vs)) => Some(NativeNum::Const(
            vs.into_iter().fold(0u32, |a, b| a ^ b) as f64
        )),
        ("bxor", None) => Some(NativeNum::BitXor(args)),
        ("bnot", Some(vs)) if vs.len() == 1 => Some(NativeNum::Const((!vs[0]) as f64)),
        ("bnot", None) if args.len() == 1 => Some(NativeNum::BitNot(Box::new(args.pop()?))),
        _ => None,
    }
}

pub(crate) fn try_compile_native_read_expr(expr: &str) -> Option<NativeExpr> {
    let ast = parse(expr)?;
    native_bool(single_return(&ast)?).map(NativeExpr)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rewrite_fast_fields() {
        let (rw, fields, needs_ud) = rewrite_fast_fields(
            r#"-- read.tid comment
return read.mapping_quality > 10 and read.bq > 20 and read.bq < 60 and "read.flags" ~= "" and foo.read.start ~= 0"#,
            "read",
            FAST_READ_FIELDS,
        );
        assert!(!needs_ud);
        assert_eq!(
            fields,
            vec![
                ("mapping_quality", FastReadField::MappingQuality),
                ("bq", FastReadField::Bq),
            ]
        );
        assert!(rw.contains("-- read.tid comment"));
        assert!(rw.contains("_pbr_mapping_quality > 10"));
        assert!(rw.contains("_pbr_bq > 20 and _pbr_bq < 60"));
        assert!(rw.contains("\"read.flags\""));
        assert!(rw.contains("foo.read.start"));

        let (rw2, fields2, needs_ud2) = rewrite_fast_fields(
            "return read.mapping_quality > 10 and read:tag('NM') < 3",
            "read",
            FAST_READ_FIELDS,
        );
        assert!(needs_ud2);
        assert_eq!(
            fields2,
            vec![("mapping_quality", FastReadField::MappingQuality)]
        );
        assert!(rw2.contains("_pbr_mapping_quality > 10 and read:tag('NM') < 3"));

        // Non-fast fields keep the userdata; assignments disable rewriting entirely.
        let (_, _, needs_ud3) =
            rewrite_fast_fields("return read.qname ~= ''", "read", FAST_READ_FIELDS);
        assert!(needs_ud3);
        let (rw4, fields4, needs_ud4) =
            rewrite_fast_fields("read.bq = 5; return read.bq > 1", "read", FAST_READ_FIELDS);
        assert!(needs_ud4 && fields4.is_empty());
        assert_eq!(rw4, "read.bq = 5; return read.bq > 1");
    }

    #[test]
    fn test_is_read_only_expression() {
        assert!(is_read_only_expression(
            "return read.mapping_quality > 10 and read.indel_count == 0 and read:tag('NM') < 3"
        ));
        assert!(!is_read_only_expression(
            "return read.mapping_quality > 10 and read.bq > 20"
        ));
        assert!(!is_read_only_expression(
            "return read.distance_from_5prime > 10"
        ));
        assert!(!is_read_only_expression(
            "print(read.qname); return read.mapping_quality > 10"
        ));
        assert!(!is_read_only_expression("return read['bq'] > 20"));
        assert!(!is_read_only_expression("return f(read)"));
        assert!(!is_read_only_expression("return read:qpos() > 3"));
        assert!(!is_read_only_expression("return read.flags >"));
    }

    #[test]
    fn test_partition_read_only_clauses() {
        let (ro, base) = partition_read_only_clauses(
            "return read.mapping_quality > 10 and read.bq > 20 and (read.flags == 0 or read.flags == 99)",
        )
        .unwrap();
        assert_eq!(
            ro,
            vec![
                "read.mapping_quality > 10",
                "(read.flags == 0 or read.flags == 99)"
            ]
        );
        assert_eq!(base, vec!["read.bq > 20"]);

        // Comments around clauses must not be captured: clauses are re-joined on one line.
        let (ro, base) = partition_read_only_clauses(
            "return read:tag('NM') < 3 -- nm\n  and read.bq > 20 -- bq\n",
        )
        .unwrap();
        assert_eq!(
            (ro, base),
            (vec!["read:tag('NM') < 3"], vec!["read.bq > 20"])
        );

        assert!(
            partition_read_only_clauses("return read.bq > 20 or read.mapping_quality > 10")
                .is_none()
        );
        assert!(partition_read_only_clauses("local x = 1; return read.bq > x").is_none());
        assert!(partition_read_only_clauses("return print(read.qname) and read.bq > 20").is_none());
        assert!(partition_read_only_clauses("return read.bq > 20 and").is_none());
    }

    #[test]
    fn test_native_compiler() {
        assert!(try_compile_native_read_expr(
            "return read.mapping_quality > 10 and read.bq > 20 and read.distance_from_5prime > 10 and read.distance_from_3prime > 10 and bit32.band(read.flags, bit32.bor(4, 256, 512, 1024)) == 0"
        )
        .is_some());
        assert!(try_compile_native_read_expr(
            "return (read.bq + 5) >= 25 and not (read.insert_size < -100 or read.length == 0)"
        )
        .is_some());
        assert!(try_compile_native_read_expr(
            "return read.mapping_quality > 10 and read:tag('NM') < 3"
        )
        .is_none());
        assert!(try_compile_native_read_expr("return read.bq").is_none());

        // `not` binds tighter than comparisons in Lua; leave these to Luau.
        assert!(try_compile_native_read_expr("return not read.flags == 0").is_none());
        assert!(try_compile_native_read_expr("return not (read.bq + 5) > 3").is_none());
        assert!(try_compile_native_read_expr("return not true == false").is_none());
        assert!(try_compile_native_read_expr("return not not (read.flags == 0)").is_some());
        // Unsupported operators fall back to Luau.
        assert!(try_compile_native_read_expr("return read.flags % 2 == 0").is_none());

        assert_eq!(
            try_compile_native_read_expr("return read.flags == 0x1_0 and 0b11 == 3"),
            try_compile_native_read_expr("return read.flags == 16 and 3 == 3"),
        );
        assert_eq!(
            try_compile_native_read_expr("return bit32.bor(4, 256) == 260"),
            try_compile_native_read_expr("return 260 == 260"),
        );
    }
}
