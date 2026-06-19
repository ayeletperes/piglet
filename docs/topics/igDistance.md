**igDistance** - *Germline set alleles distance*

Description
--------------------

Calculates the distance between pairs of alleles based on their aligned germline sequences.
Supports multiple distance methods for different segment types.


Usage
--------------------
```
igDistance(
germline_set,
AA = FALSE,
method = c("decipher", "hamming", "lv"),
trim_3prime = NULL,
return_type = c("matrix", "dist"),
quiet = TRUE
)
```

Arguments
-------------------

germline_set
:   A character vector of aligned allele sequences. See details for curating options.

AA
:   Logical (FALSE by default). If TRUE, calculate the distance based on amino acid sequences.

method
:   Distance calculation method. One of:

+  `"decipher"`: Uses DECIPHER::DistanceMatrix (requires aligned sequences, best for V segments)
+  `"hamming"`: Hamming distance (requires equal length, sequences padded if needed)
+  `"lv"`: Levenshtein distance (handles variable length, good for D/J segments)


trim_3prime
:   Optional position to trim sequences from 3' end before distance calculation

return_type
:   One of "matrix" (default) or "dist" to return a dist object

quiet
:   Logical (TRUE by default). Suppress informative messages




Value
-------------------

A `matrix` or `dist` object of the computed distances between allele pairs.


Details
-------------------

The aligned IMGT IGHV allele germline set can be downloaded from the IMGT site
[https://www.imgt.org/](https://www.imgt.org/) under the section genedb.

For V segments, the "decipher" method is recommended as it handles alignment gaps properly.
For D and J segments which may have variable lengths, the "lv" (Levenshtein) method is appropriate.



Examples
-------------------

```R
data(HVGERM)
# Using DECIPHER method (default, for V segments)
d1 <- igDistance(HVGERM[1:10], method = "decipher")

# Using Hamming distance
d2 <- igDistance(HVGERM[1:10], method = "hamming")

# Using Levenshtein distance (good for D/J segments)
d3 <- igDistance(HVGERM[1:10], method = "lv")

```



See also
-------------------

`[ighvDistance](ighvDistance.md)` for backward compatibility wrapper






