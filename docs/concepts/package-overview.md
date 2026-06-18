# Package overview

PIgLET (Program for Inferring Immunoglobulin Allele Similarity Clusters and
Genotypes) is a suite of computational tools that improves genotype inference and
downstream AIRR-seq data analysis for the immunoglobulin heavy chain V (IGHV)
locus. It has two main tools.

## Allele Similarity Cluster

Reduces ambiguity among IGHV alleles by grouping duplicated or highly similar
alleles into clusters with an unambiguous naming scheme.

| Function | Role |
| --- | --- |
| [`inferAlleleClusters()`](../reference/index.md#allele-similarity-cluster) | Main entry point: build clusters from a germline set |
| [`ighvDistance()`](../reference/index.md#allele-similarity-cluster) | Pairwise distance between aligned germline sequences |
| [`ighvClust()`](../reference/index.md#allele-similarity-cluster) | Hierarchical clustering of the distance matrix |
| [`generateReferenceSet()`](../reference/index.md#allele-similarity-cluster) | Emit the clustered reference set |
| [`artificialFRW1Germline()`](../reference/index.md#allele-similarity-cluster) | Build an IGHV reference with FWR1 primers |

The general (non-IGHV-specific) primitives `igDistance()` and `igClust()` are
also exported.

## Allele based genotype

Determines allele presence using thresholds derived from a naive population and
the allele cluster thresholds.

| Function | Role |
| --- | --- |
| [`inferGenotypeAllele()`](../reference/index.md#allele-based-genotype) | Main entry point: infer the IGHV genotype |
| [`inferGenotypeAllele_asc()`](../reference/index.md#allele-based-genotype) | Genotype inference from ASC annotations |
| [`assignAlleleClusters()`](../reference/index.md#allele-based-genotype) | Rename V calls to the cluster scheme |
| [`germlineASC()`](../reference/index.md#allele-based-genotype) | Convert an IGHV germline set to an ASC germline set |
| [`recentAlleleClusters()`](../reference/index.md#allele-based-genotype) | Download the latest cluster table from Zenodo |
| [`extractASCTable()`](../reference/index.md#allele-based-genotype) | Extract the cluster table from a Zenodo archive |
| [`zenodoArchive`](../reference/index.md#allele-based-genotype) | R6 client for the Zenodo API |

## Required input

- **For ASC clustering:** a reference set of IUIS/IMGT alleles, gapped by the
  IMGT scheme for best results.
- **For genotype inference:** an AIRR-seq format dataset.

## Reference

Peres et al. (2023), "Immunoglobulin Allele Similarity Clusters", *Nucleic Acids
Research*, [doi:10.1093/nar/gkad603](https://doi.org/10.1093/nar/gkad603).
