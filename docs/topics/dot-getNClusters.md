**.getNClusters** - *Find resolution for target cluster count*

Description
--------------------

Uses binary search to find a resolution parameter that produces approximately
the target number of clusters.


Usage
--------------------
```
.getNClusters(
g,
n_cluster,
range_min = 0,
range_max = 6,
max_steps = 20,
method = "leiden"
)
```

Arguments
-------------------

g
:   An igraph graph object with weighted edges

n_cluster
:   Target number of clusters

range_min
:   Minimum resolution to search. Default is 0.

range_max
:   Maximum resolution to search. Default is 6.

max_steps
:   Maximum number of search iterations. Default is 20.

method
:   Community detection method: "leiden" or "louvain". Default is "leiden".




Value
-------------------

A list containing:

+  `partition`: The community detection result
+  `clusters`: Number of clusters found
+  `best_resolution`: The resolution parameter used










