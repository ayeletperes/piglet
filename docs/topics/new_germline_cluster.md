**new_germline_cluster** - *Create a GermlineCluster object*

Description
--------------------

`GermlineCluster` is an S3 class that stores the output of
`[inferAlleleClusters](inferAlleleClusters.md)`. It contains the allele cluster table,
clustering objects, and threshold parameters used for inference.


Usage
--------------------
```
new_germline_cluster(
germlineSet,
alleleClusterSet,
alleleClusterTable,
threshold,
hclustAlleleCluster = NULL,
clusteringMethod = "hierarchical",
communityObject = NULL,
graphObject = NULL,
distanceMatrix = NULL,
silhouetteScore = NA_real_,
resolutionParameter = NA_real_,
locus = "IGHV",
.familiesCut = NULL
)
```

Arguments
-------------------

germlineSet
:   The original germline set provided.

alleleClusterSet
:   The renamed germline set with allele clusters.

alleleClusterTable
:   The allele cluster table.

threshold
:   The threshold used for family and allele clusters.

hclustAlleleCluster
:   A hierarchical clustering object for the germline set,
or `NULL`.

clusteringMethod
:   The clustering method used, either `"hierarchical"`
or `"leiden"`.

communityObject
:   A community detection object for Leiden clustering, or `NULL`.

graphObject
:   An igraph graph object for Leiden clustering, or `NULL`.

distanceMatrix
:   The distance matrix used for clustering, or `NULL`.

silhouetteScore
:   The silhouette score for community detection.

resolutionParameter
:   The resolution parameter used for Leiden clustering.

locus
:   The locus identifier, for example `"IGHV"`, `"IGHD"`, `"IGHJ"`.

.familiesCut
:   The family-level cut data stored with the clustering result,
or `NULL`.




Value
-------------------

An object of class `"GermlineCluster"`.




See also
-------------------

`[inferAlleleClusters](inferAlleleClusters.md)`

`[GermlineCluster](GermlineCluster.md)`






