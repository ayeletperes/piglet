**optimize_resolution** - *Optimize resolution parameter using silhouette score*

Description
--------------------

Performs a grid search over resolution parameters and selects the one
that maximizes the silhouette score.


Usage
--------------------
```
optimize_resolution(
g,
distance_matrix,
target_clusters = 80,
resolution_range_low = 0.1,
resolution_range_high = 0.5,
max_steps = 20,
ncores = 1
)
```

Arguments
-------------------

g
:   An igraph graph object with weighted edges

distance_matrix
:   The distance matrix (as dist object) used for silhouette calculation

target_clusters
:   Target number of clusters for initial tuning. Default is 80.

resolution_range_low
:   Fractional range below tuned resolution. Default is 0.1.

resolution_range_high
:   Fractional range above tuned resolution. Default is 0.5.

max_steps
:   Maximum steps for initial tuning. Default is 20.

ncores
:   Number of cores for parallel processing. Default is 1.




Value
-------------------

A list containing:

+  `results`: data.frame with Resolution, ClusterCount, Silhouette
+  `partitions`: list of membership vectors for each resolution
+  `best_resolution`: optimal resolution parameter
+  `best_partition`: membership vector at optimal resolution
+  `best_clusters`: number of clusters at optimal resolution





See also
-------------------

`[detect_communities_leiden](detect_communities_leiden.md)`, `[igClust](igClust.md)`






