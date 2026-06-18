# Build allele similarity clusters

This guide walks through producing an Allele Similarity Cluster (ASC) set from an
IGHV germline reference. For the concepts behind it, read the
[ASC lesson](../lessons/allele-similarity-clusters.md).

## 1. Get a gapped germline reference

ASC clustering works best with IUIS/IMGT alleles gapped by the IMGT scheme.
PIgLET ships the human IGHV reference as `HVGERM`; you can also pass a character
vector of sequences or a path to a FASTA file.

```r
library(piglet)
germline_set <- HVGERM
```

## 2. Run the clustering

```r
clusters <- inferAlleleClusters(
  germline_set             = germline_set,
  trim_3prime_side         = 318,   # trim position; NULL keeps full length
  mask_5prime_side         = 0,     # mask 5' nucleotides to mimic short libraries
  family_threshold         = 75,    # family-level similarity %
  allele_cluster_threshold = 95     # allele-cluster-level similarity %
)
```

Under the hood this runs `ighvDistance()` → `ighvClust()` →
`generateReferenceSet()`. You can call those individually for finer control.

## 3. Inspect and plot

```r
clusters          # print summary
summary(clusters)
plot(clusters)    # allele cluster dendrogram
```

The return value is a [`GermlineCluster`](../concepts/germline-cluster.md)
object.

## 4. Use the clustered names downstream

Rename V allele calls in an AIRR dataset to the cluster scheme with
`assignAlleleClusters()`, or convert an existing germline set to an ASC germline
set with `germlineASC()`.

## Full walkthrough

The complete, runnable narrative (with figures) lives in the
**[PIgLET vignette](../vignettes/PIgLET-vignette.md)**.
