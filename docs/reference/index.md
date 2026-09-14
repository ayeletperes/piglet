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
- [`regroupASCByLabel`](../topics/regroupASCByLabel.md) — re-group an existing
  assignment by a per-allele label, such as the IUIS gene.
- [`ascIUISVocabulary`](../topics/ascIUISVocabulary.md) — build the ASC-to-IUIS
  display labels, and the allele set each one stands for.
- [`annotateRepertoireIUIS`](../topics/annotateRepertoireIUIS.md) — apply those
  labels to a repertoire.

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

## Gene-usage QTL

Does a germline variant change how much a gene group is used? The phenotype is the
group's share of its own segment in one subject.

- [`runGeneUsageQTL`](../topics/runGeneUsageQTL.md) — main entry point: phenotype,
  scan, threshold and leads for one locus.
- [`ascUsagePhenotype`](../topics/ascUsagePhenotype.md) — per-subject usage of each
  gene group.
- [`qtlScanUnivariate`](../topics/qtlScanUnivariate.md) — ordinary least squares of
  every phenotype on every variant.
- [`qtlScanMultivariate`](../topics/qtlScanMultivariate.md) — the same against a
  multivariate response, by Pillai's trace.
- [`qtlLDGroups`](../topics/qtlLDGroups.md) — collapse variants carrying identical
  genotype information, which is what a threshold should be set from.
- [`qtlClump`](../topics/qtlClump.md) — independent lead variants.
- [`qtlSmallestGenotypeClass`](../topics/qtlSmallestGenotypeClass.md) — the smallest
  observed genotype class, which a p-value should never be read without.

## Conditional gene pairing

Does a variant change *which* partner a gene pairs with, rather than only how often
it is used? Both marginals are divided out, so the phenotype is pairing preference
rather than usage.

- [`runPairingQTL`](../topics/runPairingQTL.md) — main entry point, one anchoring
  per call.
- [`pairingTable`](../topics/pairingTable.md) — observed against expected pairing,
  per subject.
- [`pairingPhenotype`](../topics/pairingPhenotype.md) — the response matrix one
  anchor is scanned on.
- [`pairingScan`](../topics/pairingScan.md) — variants against an anchor's partner
  distribution.
- [`pairingCellTests`](../topics/pairingCellTests.md) — which partners moved, under
  a significant hit.
- [`pairingLeadCharacter`](../topics/pairingLeadCharacter.md) — whether a hit is a
  reallocation, a rigid shift, one subject, or marginal usage in disguise.

## Data

- [`HVGERM`](../topics/HVGERM.md) — human IGHV germline reference.
- [`hv_functionality`](../topics/hv_functionality.md) — allele functionality
  annotations.

## SNP and sequence utilities

- [`allele_diff`](../topics/allele_diff.md) — differences between aligned alleles.
- [`allele_diff_strings`](../topics/allele_diff_strings.md) — as strings.
- [`allele_diff_indices`](../topics/allele_diff_indices.md) — as positions.
- [`allele_diff_paired`](../topics/allele_diff_paired.md) — count or locate SNPs
  between paired germline and input sequences.
- [`insert_gaps`](../topics/insert_gaps.md) — apply a reference gap pattern to an
  ungapped sequence.

### Deprecated

Still exported and working; each delegates to its replacement.

- [`allele_diff_indices_parallel2`](../topics/allele_diff_indices_parallel2.md) —
  use `allele_diff_paired`. Results unchanged.
- [`allele_diff_indices_parallel`](../topics/allele_diff_indices_parallel.md) —
  use `allele_diff_paired`. **Its counts change**: the old version ignored gaps only
  on the germline side and read past the end of shorter inputs.
- [`insert_gaps2_vec`](../topics/insert_gaps2_vec.md) — use `insert_gaps`. Results
  unchanged.
