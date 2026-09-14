**allele_diff_paired** - *Count or locate SNPs between paired germline and input sequences*

Description
--------------------

Compares each germline sequence in `germs` against the input sequence at the
same position in `inputs`, and returns either the positions of the mismatches
or how many there are. The two vectors are paired element by element, which is what
distinguishes this from `[allele_diff_indices](allele_diff_indices.md)`, where every sequence is
compared against the first.


Usage
--------------------
```
allele_diff_paired(
germs,
inputs,
X = 0L,
return_count = FALSE,
non_mismatch_chars_nullable = NULL
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

return_count
:   Return the number of mismatches per pair rather than their
positions (default `FALSE`).

non_mismatch_chars_nullable
:   Characters ignored on either side when comparing
(default `N`, `.`, `-`).




Value
-------------------

A list of integer vectors of one-based positions when
`return_count = FALSE`, or an integer vector of counts when `TRUE`.


Details
-------------------

A position is ignored when *either* sequence carries a non-mismatch character,
by default a gap or an ambiguous base (`N`, `.`, `-`). Sequences of
unequal length are padded to the longer with `N`, so the padded positions never
count as mismatches.



Examples
-------------------

```R
germs  <- c("ACGTACGT", "ACGTACGT")
inputs <- c("ACGTTCGT", "ACGTTCGA")

allele_diff_paired(germs, inputs, return_count = TRUE)   # 1, 2

```


```
[1] 1 2

```


```R
allele_diff_paired(germs, inputs)                        # list(5), c(5, 8)

```


```
[[1]]
[1] 5

[[2]]
[1] 5 8


```








