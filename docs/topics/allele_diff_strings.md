**allele_diff_strings** - *Calculate differences between characters in columns of germs and return them as a string vector.*

Description
--------------------

Calculate differences between characters in columns of germs and return them as a string vector.


Usage
--------------------
```
allele_diff_strings(germs, X = 0L, non_mismatch_chars_nullable = NULL)
```

Arguments
-------------------

germs
:   A vector of strings representing germ sequences.

X
:   The threshold index from which to return differences as strings.

non_mismatch_chars_nullable
:   A set of characters that are ignored when comparing sequences (default: 'N', '.', '-').




Value
-------------------

A vector of strings containing differences between characters in columns.



Examples
-------------------

```R
germs = c("ATCG", "ATCC") 
X = 3 
result = allele_diff_strings(germs, X) 
# "A2T", "T3C", "C2G"

```








