# Quick start

This page runs both PIgLET tools end-to-end. It assumes you have already
[installed](installation.md) the package.

```r
library(piglet)
```

## 1. Build an allele similarity cluster (ASC) set

`inferAlleleClusters()` takes a gapped IGHV germline reference and groups
similar/duplicated alleles into clusters, producing an unambiguous renaming
scheme. The package ships the human IGHV reference as `HVGERM`.

```r
clusters <- inferAlleleClusters(
  germline_set            = HVGERM,
  trim_3prime_side        = 318,
  family_threshold        = 75,
  allele_cluster_threshold = 95
)

clusters          # prints a summary of the GermlineCluster object
plot(clusters)    # dendrogram of the clustering
```

The result is a [`GermlineCluster`](../concepts/germline-cluster.md) object that
holds the renamed germline sequences, the cluster assignments, and the threshold
values used.

## 2. Infer an allele-based genotype

`inferGenotypeAllele()` decides which alleles are genuinely present in a subject
by comparing each allele's frequency against a population-derived threshold
(supplied through an allele cluster table). It expects AIRR-format data with a
`v_call` column.

```r
genotype <- inferGenotypeAllele(
  data               = airr_data,             # AIRR-format rearrangements
  alleleClusterTable = allele_threshold_table # thresholds (shipped with package)
)

head(genotype)
```

## Where to go next

- Understand the ideas: **[Lessons](../lessons/allele-similarity-clusters.md)**.
- Full walkthroughs: **[Task Guides](../guides/build-allele-clusters.md)**.
- Function-by-function docs: **[API Reference](../reference/index.md)**.
