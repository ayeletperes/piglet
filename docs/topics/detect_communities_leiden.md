**detect_communities_leiden** - *Leiden community detection*

Description
--------------------

Performs community detection on a weighted graph using the Leiden algorithm
with CPM (Constant Potts Model) objective function.


Usage
--------------------
```
detect_communities_leiden(g, resolution = 1)
```

Arguments
-------------------

g
:   An igraph graph object with weighted edges

resolution
:   Resolution parameter for Leiden algorithm. Higher values
produce more communities. Default is 1.0.




Value
-------------------

An igraph communities object


Details
-------------------

The Leiden algorithm is a community detection method that optimizes a quality
function (here CPM). It guarantees connected communities and is generally
faster than Louvain while producing better quality partitions.



Examples
-------------------

```R
data(HVGERM)
d <- igDistance(HVGERM[1:10], method = "hamming")
g <- distance_to_graph(d)
comm <- detect_communities_leiden(g, resolution = 0.5)

```



See also
-------------------

`[distance_to_graph](distance_to_graph.md)`, `[optimize_resolution](optimize_resolution.md)`






