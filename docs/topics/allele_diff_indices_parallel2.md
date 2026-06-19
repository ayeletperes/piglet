**allele_diff_indices_parallel2** - *Calculate SNPs or their count for each germline-input sequence pair with optional parallel execution.*

Description
--------------------

This function compares germline sequences (`germs`) and input sequences (`inputs`)
and identifies single nucleotide polymorphisms (SNPs) or their counts, with optional parallel execution.
The comparison ignores specified non-mismatch characters (e.g., gaps or ambiguous bases).


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
:   A vector of strings representing input sequences.

X
:   The threshold index from which to return SNP indices or counts (default: 0).

parallel
:   A boolean flag to enable parallel processing (default: FALSE).

return_count
:   A boolean flag to return the count of mutations instead of their indices (default: FALSE).

non_mismatch_chars_nullable
:   A set of characters that are ignored when comparing sequences (default: 'N', '.', '-').




Value
-------------------

A list of integer vectors (if `return_count = FALSE`) or a vector of integers (if `return_count = TRUE`).



Examples
-------------------

```R
# Example usage
germs <- c("ATCG", "ATCC")
inputs <- c("ATTG", "ATTA")
X <- 0

# Return indices of SNPs
result_indices <- allele_diff_indices_parallel2(germs, inputs, X, 
parallel = TRUE, return_count = FALSE)
print(result_indices)  # list(c(4), c(3, 4))

```


```
[[1]]
[1] 3

[[2]]
[1] 3 4


```


```R

# Return counts of SNPs
result_counts <- allele_diff_indices_parallel2(germs, inputs, X, 
parallel = FALSE, return_count = TRUE)
print(result_counts)  # c(1, 2)

```


```
[1] 1 2

```








