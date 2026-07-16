**alleleClusterNames** - *Allele similarity cluster naming scheme*

Description
--------------------

For a given cluster the function collapse similar sequences and renames the sequences based on the ASC name scheme


Usage
--------------------
```
alleleClusterNames(
cluster,
allele.cluster.table,
germ.dist,
chain,
segment,
family_prefix = TRUE,
retain_subgroup = FALSE
)
```

Arguments
-------------------

cluster
:   A vector with the cluster identifier - the family and allele cluster number.

allele.cluster.table
:   A data.frame with the list of all germline sequences and their clusters.

germ.dist
:   A matrix with the germline distance between the germline set sequences.

chain
:   A character with the chain identifier: IGH/IGL/IGK/TRB/TRA... (Currently only IGH is supported)

segment
:   A character with the segment identifier: IGHV/IGHD/IGHJ.... (Currently only IGHV is supported)

family_prefix
:   Logical. If TRUE (default), prepend "F" to the family number in ASC names (e.g. IGHVF1-G1*01). If FALSE, omit the "F" (e.g. IGHV1-G1*01).

retain_subgroup
:   Logical. If TRUE, retain the original IMGT subgroup in the ASC name instead of renumbering by the family clustering (e.g. IGKV1-12*01 stays in subgroup IGKV1-G..). When TRUE, `family_prefix` is ignored. Default FALSE.




Value
-------------------

A data.frame with the clusters renamed alleles based on the ASC scheme.









