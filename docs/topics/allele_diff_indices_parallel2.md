**allele_diff_indices_parallel2** - *Deprecated: use `[allele_diff_paired](allele_diff_paired.md)`*

Description
--------------------

Renamed. The old name said `parallel`, which was an argument rather than a
behaviour and never had any effect (the package is not built with OpenMP), and
`indices`, which is wrong whenever `return_count = TRUE`.


Usage
--------------------
```
allele_diff_indices_parallel2(
germs,
inputs,
X = 0L,
parallel = FALSE,
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

parallel
:   Ignored. Accepted so existing calls keep working.

return_count
:   Return the number of mismatches per pair rather than their
positions (default `FALSE`).

non_mismatch_chars_nullable
:   Characters ignored on either side when comparing
(default `N`, `.`, `-`).




Value
-------------------

As `[allele_diff_paired](allele_diff_paired.md)`.


Details
-------------------

Results are unchanged: this is the same implementation under a clearer name.




See also
-------------------

`[allele_diff_paired](allele_diff_paired.md)`






