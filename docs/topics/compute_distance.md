**compute_distance** - *Compute distance matrix*

Description
--------------------

Compute a pairwise distance matrix between sequences using stringdist.


Usage
--------------------
```
compute_distance(
sequences,
method = c("hamming", "lv"),
trim_3prime = NULL,
quiet = TRUE,
return_type = c("dist", "matrix")
)
```

Arguments
-------------------

sequences
:   A named character vector of sequences

method
:   Distance method: "hamming" or "lv" (Levenshtein). Default is "hamming".

trim_3prime
:   Optional position to trim sequences from 3' end

quiet
:   Logical. Suppress messages. Default is TRUE.

return_type
:   One of "dist" (default) or "matrix"




Value
-------------------

A dist object or matrix of pairwise distances




See also
-------------------

`[igDistance](igDistance.md)` for more distance options






