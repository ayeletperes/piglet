# API Reference

Function-by-function documentation, generated from the package's roxygen
comments. Grouped by tool.

## Allele Similarity Cluster

- [`inferAlleleClusters`](../topics/inferAlleleClusters.md) — main entry point:
  build clusters from a germline set.
- [`ighvDistance`](../topics/ighvDistance.md) — distance between aligned IGHV
  germline sequences.
- [`ighvClust`](../topics/ighvClust.md) — hierarchical clustering of the distance
  matrix.
- [`generateReferenceSet`](../topics/generateReferenceSet.md) — generate the
  clustered reference set.
- [`plotAlleleCluster`](../topics/plotAlleleCluster.md) — plot the hierarchical
  clustering.
- [`artificialFRW1Germline`](../topics/artificialFRW1Germline.md) — build an IGHV
  reference with FWR1 primers.
- [`alleleClusterNames`](../topics/alleleClusterNames.md) — work with allele
  cluster names.
- [`GermlineCluster-class`](../topics/GermlineCluster-class.md) — the S4 output
  object.

## Allele based genotype

- [`inferGenotypeAllele`](../topics/inferGenotypeAllele.md) — main entry point:
  infer the IGHV genotype.
- [`assignAlleleClusters`](../topics/assignAlleleClusters.md) — rename V calls to
  the cluster scheme.
- [`germlineASC`](../topics/germlineASC.md) — convert an IGHV germline set to an
  ASC germline set.
- [`recentAlleleClusters`](../topics/recentAlleleClusters.md) — download the
  latest cluster table from Zenodo.
- [`extractASCTable`](../topics/extractASCTable.md) — extract the cluster table
  from a Zenodo archive.
- [`zenodoArchive`](../topics/zenodoArchive.md) — R6 client for the Zenodo API.

## Data

- [`HVGERM`](../topics/HVGERM.md) — human IGHV germline reference.
- [`hv_functionality`](../topics/hv_functionality.md) — allele functionality
  annotations.
