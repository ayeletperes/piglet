# Lesson: Allele-based genotyping

*A concept primer. Read the [ASC lesson](allele-similarity-clusters.md) first if
allele ambiguity is new to you.*

## The problem: which alleles does a subject actually carry?

A genotype is the set of germline alleles an individual carries. Inferring it
from AIRR-seq data is hard because low-frequency allele calls can be **artefacts**
— sequencing error, misalignment, or somatic hypermutation mimicking another
allele — rather than evidence that the allele is genuinely present.

A naive rule ("call it present if it appears at all") over-calls. A fixed
frequency cutoff under-calls rare-but-real alleles and over-calls error-prone
ones, because every allele has a different error profile.

## The idea: per-allele thresholds from a naive population

PIgLET's second tool, the **allele-based genotype**, sets an **allele-specific
threshold** learned from a naive population. For each allele, the absolute
frequency observed in a subject is checked against that allele's own threshold.
If it clears the bar, the allele is called present.

Because the thresholds are tied to the allele similarity clusters, this approach
also benefits from the disambiguation described in the ASC lesson.

## How PIgLET computes it

The entry point is **`inferGenotypeAllele()`**. You give it:

- `data` — AIRR-format rearrangements for one subject, with a `v_call` column.
- `alleleClusterTable` — the table of per-allele thresholds.

Useful options:

- `single_assignment = TRUE` restricts inference to reads with a single
  (unambiguous) V call.
- `find_unmutated = TRUE` (with `germline_db`) restricts to unmutated reads.

An ASC-aware variant, `inferGenotypeAllele_asc()`, works directly from ASC
annotations.

## Getting the threshold table

The threshold and cluster tables can be downloaded from Zenodo via
`recentAlleleClusters()` / `extractASCTable()` (using the `zenodoArchive` R6
client), or you can use the versions shipped with the package.

## Next steps

- Do it: **[Infer an allele-based genotype](../guides/infer-genotype.md)**.
- Function reference:
  **[`inferGenotypeAllele`](../reference/index.md#allele-based-genotype)**.
