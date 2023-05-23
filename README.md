# pbr
drunk on [perbase](https://github.com/sstadick/perbase) pileups and lua expressions.

This should be able to:

1. Exclude called sites with ≥5% N frequency (no call rate)
2. Trim terminal ends of reads (15 bases 5' and 3')
3. Exclude coverage at sites overlapping homopolymers (≥4bp in length) and fraguracy
4. Exclude reads with ≥5% N’s on them

We can do these as follows:

1. Remove columns where %N > 0. (post-processing)
2. TBD
3. Use exclude regions
4. use lua expression: `string_count(read.sequence, 'N') < 0.05 * #read.sequence`


For now, this runs as:

```
 cargo run $bam "string_count(read.sequence, 'N') < 0.05 * #read.sequence and read.mapping_quality > 30"
```
where the string argument is the lua expression (and yes, #read.sequence gives length of the read sequence) for filtering reads.
