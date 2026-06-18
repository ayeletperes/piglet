# Lesson: Allele Similarity Clusters

*A concept primer. No prior PIgLET knowledge assumed.*

## The problem: ambiguous IGHV allele calls

When AIRR-seq reads are aligned to an IGHV germline reference, a single read can
match several reference alleles equally well. This happens because the IGHV locus
contains alleles that are **duplicated or highly similar across different genes**.
The aligner then reports an ambiguous call (multiple alleles separated by commas),
which propagates uncertainty into every downstream analysis — genotyping,
repertoire comparison, and lineage work.

## The idea: cluster similar alleles

PIgLET's first tool, the **Allele Similarity Cluster (ASC)**, groups alleles
that cannot be reliably distinguished into clusters, and assigns each cluster a
single unambiguous name. Instead of asking "exactly which of these near-identical
alleles is it?", you ask "which cluster is it?" — a question the data can
actually answer.

Clustering happens at two levels:

1. **Family cluster** — a coarse grouping (default similarity threshold 75%).
2. **Allele cluster** — a fine grouping within a family (default 95%).

## The naming scheme

Each clustered allele is renamed following:

```
IGHVF1-G1*01
    |   |   |
    |   |   +-- allele number (by clustering order)
    |   +------ allele cluster (G) number
    +---------- family cluster (F) number
```

`IGH` is the chain, `V` the region. So `IGHVF1-G1*01` reads as: heavy chain, V
region, family cluster 1, allele cluster 1, allele 1.

## How PIgLET computes it

The workflow behind `inferAlleleClusters()` is three steps:

1. **`ighvDistance()`** — compute pairwise distances between the aligned germline
   sequences.
2. **`ighvClust()`** — hierarchically cluster that distance matrix using the two
   thresholds.
3. **`generateReferenceSet()`** — emit the renamed, de-duplicated germline set.

You normally call the wrapper `inferAlleleClusters()` and let it orchestrate all
three.

## Next steps

- Do it: **[Build allele similarity clusters](../guides/build-allele-clusters.md)**.
- The object it produces: **[`GermlineCluster`](../concepts/germline-cluster.md)**.
- The companion tool: **[Allele-based genotyping](allele-based-genotype.md)**.

PIgLET is described in Peres et al. (2023), *Nucleic Acids Research*,
[doi:10.1093/nar/gkad603](https://doi.org/10.1093/nar/gkad603).
