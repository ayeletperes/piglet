**allele_diff_indices_parallel** - *Calculate SNPs or their count for each germline-input sequence pair with optional parallel execution.*

Description
--------------------

Calculate SNPs or their count for each germline-input sequence pair with optional parallel execution.


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
:   A vector of strings representing input sequences.

X
:   The threshold index from which to return SNP indices or counts (default: 0).

parallel
:   A boolean flag to enable parallel processing (default: FALSE).

return_count
:   A boolean flag to return the count of mutations instead of their indices (default: FALSE).




Value
-------------------

A list of integer vectors (if return_count = FALSE) or a vector of integers (if return_count = TRUE).









