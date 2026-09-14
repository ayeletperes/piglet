**allele_diff_indices_parallel** - *Deprecated: use `[allele_diff_paired](allele_diff_paired.md)`*

Description
--------------------

Renamed, and **the counts change**. The old implementation ignored gaps and
ambiguous bases only on the germline side, so an aligned position where the germline
carried a base and the input carried a gap was counted as a mismatch. It also looped
to the length of the germline while indexing the input, reading past the end of the
input whenever the input was shorter -- undefined behaviour, and on real IGHV data it
returned counts derived partly from adjacent memory.


Usage
--------------------
```
allele_diff_indices_parallel(
germs,
inputs,
X = 0L,
parallel = FALSE,
return_count = FALSE
)
```

Arguments
-------------------

germs
:   A vector of strings representing germline sequences.

inputs
:   A vector of strings representing input sequences, paired with
`germs` element by element.

X
:   The position from which mismatches are counted, zero-based (default 0).

parallel
:   Ignored. Accepted so existing calls keep working.

return_count
:   Return the number of mismatches per pair rather than their
positions (default `FALSE`).




Value
-------------------

As `[allele_diff_paired](allele_diff_paired.md)`.


Details
-------------------

Calls now route to `[allele_diff_paired](allele_diff_paired.md)`, which ignores such positions on
either side and pads the shorter sequence. Against the old behaviour, 194 of 400 real
germline pairs differ, by a mean of 8 mismatches. Anything that depended on the old
numbers needs re-checking rather than re-running.




See also
-------------------

`[allele_diff_paired](allele_diff_paired.md)`






