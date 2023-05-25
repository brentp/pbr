# pbr
drunk on [perbase](https://github.com/sstadick/perbase) pileups and lua expressions.

## Expression

The following attributes are available on the `read` object in lua expressions:

```        
mapping_quality
flags # integer. must use bit32 lua module (builtin here) to do operations
tid
pos
start
stop
# where in the current read is the pileup given by qpos with convenience of distance_from_[left/right]_end
qpos
distance_from_left_end 
distance_from_right_end
insert_size
qname
base_qualities # list of base_qualities
bq # base_quality at current site
sequence # read sequence
cigar_string
length # length of the read sequence
```

An example expression could be:

```lua
-- high mapping-q    and base-quality for this column 
read.mapping_quality > 10 and read.bq > 20 \
--  and  not within 10 bases of left end      or     right end 
    and read.distance_from_left_end > 10 and read.distance_from_right_end > 10 \
--  and  exclude read if unmapped, not primary, qc_fail, or duplicate. 
    and bit32.band(read.flags, bit32.bor(4, 256, 512, 1024)) == 0 \
--  and exclude read if it has more than 5% N's in the sequence
    and string_count(read.sequence, 'N') < 0.05 * read.length
```


For now, this runs as:

```
 cargo run $bam "return $expression"
```
where the $expression argument is the lua expression.

+ Note that we can use, e.g. `print(read.qname, read.flags); return $expression)` to help with debugging.
+ Note that the expression *must* contain 'return'
