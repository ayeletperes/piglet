**allele_diff_indices** - *Calculate differences between characters in columns of germs and return their indices as an int vector.*

Description
--------------------

Calculate differences between characters in columns of germs and return their indices as an int vector.


Usage
--------------------
```
allele_diff_indices(germs, X = 0L, non_mismatch_chars_nullable = NULL)
```

Arguments
-------------------

germs
:   A vector of strings representing germ sequences.

X
:   The threshold index from which to return differences as indices.

non_mismatch_chars_nullable
:   A set of characters that are ignored when comparing sequences (default: 'N', '.', '-').




Value
-------------------

A vector of integers containing indices of differing columns.



Examples
-------------------

```R
germs = c("ATCG", "ATCC") 
X = 3 
result = allele_diff_indices(germs, X)
# 1, 2, 3

```








