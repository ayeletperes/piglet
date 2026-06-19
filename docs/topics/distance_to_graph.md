**distance_to_graph** - *Convert distance matrix to weighted graph*

Description
--------------------

Converts a distance matrix to a weighted igraph object using a log transform
that spreads small distances and produces weights in [0,1].


Usage
--------------------
```
distance_to_graph(distance_matrix)
```

Arguments
-------------------

distance_matrix
:   A distance matrix or dist object




Value
-------------------

An igraph object with weighted edges


Details
-------------------

The transformation uses a log-based similarity measure:

1.  Normalize distances by max distance
1.  Apply -log transform to convert to similarity
1.  Normalize similarities to [0,1] range
1.  Create weighted undirected graph




Examples
-------------------

```R
data(HVGERM)
d <- igDistance(HVGERM[1:10], method = "hamming")
g <- distance_to_graph(d)

```



See also
-------------------

`[detect_communities_leiden](detect_communities_leiden.md)`, `[igClust](igClust.md)`






