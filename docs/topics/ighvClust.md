**ighvClust** - *Allele similarity clustering (deprecated)*

Description
--------------------

This function is deprecated. Use `[igClust](igClust.md)` instead.


Usage
--------------------
```
ighvClust(
germline_distance,
family_threshold = 75,
allele_cluster_threshold = 95,
cluster_method = "complete"
)
```

Arguments
-------------------

germline_distance
:   A germline set distance matrix created by `[igDistance](igDistance.md)`.

family_threshold
:   The threshold for family-level grouping.
For `distance_method = "decipher"`: similarity percentage (default 75).
For `"hamming"` / `"lv"`: maximum number of mismatches (must be >= allele_cluster_threshold).

allele_cluster_threshold
:   The threshold for allele-cluster-level grouping (hierarchical only).
For `distance_method = "decipher"`: similarity percentage (default 95).
For `"hamming"` / `"lv"`: maximum number of mismatches (must be <= family_threshold).

cluster_method
:   The linkage method for hierarchical clustering (used for family assignment in both methods). Default is "complete".




Value
-------------------

A named list with clustering results.




See also
-------------------

`[igClust](igClust.md)` for the current implementation






