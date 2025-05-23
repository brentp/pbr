# pbr

[![Rust](https://github.com/brentp/pbr/actions/workflows/rust.yml/badge.svg)](https://github.com/brentp/pbr/actions/workflows/rust.yml)
drunk on [perbase](https://github.com/sstadick/perbase) pileups and [lua](https://github.com/khvzak/mlua/) expressions.

An example use is in calculating mutation rates. Variant callers have filters such as these
that affect the number of mutations called; we also want to scale the denominator in similar
way so we can compare rates across samples.

## Expression

The following attributes are available on the `read` object in lua expressions:

```
mapping_quality
flags # integer. must use bit32 lua module (builtin here) to do operations
tid
start
stop
# where in the current read is the pileup given by qpos with convenience of distance_from_[left/right]_end
qpos
distance_from_5prime
distance_from_3prime
insert_size
qname
bq # base_quality at current site
length # length of the read sequence
sequence
n_proportion_5_prime(bases:number)
n_proportion_3_prime(bases:number)
indel_count
soft_clips_3_prime
soft_clips_5_prime
tag(name: string)
```

An example expression could be:

```lua
-- high mapping-q    and base-quality for this column
read.mapping_quality > 10 and read.bq > 20 \
--  and  not within 10 bases of left end      or     right end
    and read.distance_from_5prime > 10 and read.distance_from_3prime > 10 \
--  and  exclude read if unmapped, not primary, qc_fail, or duplicate.
    and bit32.band(read.flags, bit32.bor(4, 256, 512, 1024)) == 0 \
--  and exclude read if it has more than 5% N's in the sequence
    and string_count(read.sequence, 'N') < 0.05 * read.length
```

this runs as:

```
pbr $bam "return $expression" > out.pileup
```

where the $expression argument is the lua expression.

- Note that we can use, e.g. `print(read.qname, read.flags); return $expression)` to help with debugging.
- Note that the expression _must_ contain **'return'**

# Usage

```
pileups filtered with lua expressions

Usage: pbr [OPTIONS] <BAM_PATH> <EXPRESSION>

Arguments:
  <BAM_PATH>    Path to the bamfile
  <EXPRESSION>  Lua expression to evaluate

Options:
  -t, --threads <THREADS>                  Number of threads to use [default: 2]
  -m, --max-depth <MAX_DEPTH>              maximum depth in the pileup [default: 100000]
  -b, --bedfile <BEDFILE>                  optional path to the BED of include regions
  -f, --fasta <FASTA>                      optional path to the reference fasta file
  -e, --exclude <EXCLUDE>                  optional path to BED of exclude regions
      --mate-fix                           adjust depth to not double count overlapping mates
  -p, --pile-expression <PILE_EXPRESSION>  optional expression required for the pileup
  -h, --help                               Print help
  -V, --version                            Print version
```

## PileExpression

Note that the pile-expression is also a lua expression; it is applied to the Pileup (column) rather than to the reads.
The available attributes on the `pile` object are:

```
depth,a,c,g,t,n,fail,ins,del,ref_skip
```

An example --pile-expression would look like:

```
return pile.n / pile.depth < 0.05
```

To require that fewer than 5% of the reads in the pile are 'N'. Positions that do not pass this expression will **not** be printed.
