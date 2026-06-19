**igClust** - *Allele similarity clustering*

Description
--------------------

Cluster the distance matrix to create allele clusters. Supports both
hierarchical clustering (default) and Leiden community detection.


Usage
--------------------
```
igClust(
germline_distance,
method = c("hierarchical", "leiden"),
distance_method = "decipher",
family_threshold = 75,
allele_cluster_threshold = 95,
cluster_method = "complete",
resolution = NULL,
target_clusters = NULL,
optimize_silhouette = TRUE,
ncores = 1,
quiet = FALSE
)
```

Arguments
-------------------

germline_distance
:   A germline set distance matrix created by `[igDistance](igDistance.md)`.

method
:   Clustering method. One of "hierarchical" (default) or "leiden".

distance_method
:   The distance method used to compute `germline_distance`.
One of `"decipher"` (default), `"hamming"`, or `"lv"`.
For `"decipher"`, thresholds are similarity percentages in `[0, 100]`;
for `"hamming"` and `"lv"`, thresholds are raw integer mismatch counts.

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

resolution
:   Resolution parameter for Leiden clustering. If NULL, will be optimized.

target_clusters
:   Target number of clusters for Leiden optimization. Default is NULL.

optimize_silhouette
:   Logical. Optimize resolution using silhouette score (Leiden only). Default is TRUE.

ncores
:   Number of cores for parallel processing (Leiden only). Default is 1.

quiet
:   Logical. Suppress messages. Default is FALSE.




Value
-------------------

A named list that includes:

+  `alleleClusterTable`: data.frame of allele clusters
+  `threshold`: list of threshold parameters
+  `hclustAlleleCluster`: hierarchical clustering object (both methods)
+  `communityObject`: community detection result (Leiden method)
+  `graphObject`: igraph object (Leiden method)
+  `silhouetteScore`: silhouette score (Leiden method)
+  `resolutionParameter`: resolution used (Leiden method)





See also
-------------------

`[igDistance](igDistance.md)`, `[inferAlleleClusters](inferAlleleClusters.md)`






