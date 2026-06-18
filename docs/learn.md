# Choose your path

PIgLET (Program for Inferring Immunoglobulin Allele Similarity Clusters and
Genotypes) provides two complementary tools for AIRR-seq analysis of the
immunoglobulin heavy chain V (IGHV) locus. Pick the entry point that matches
what you are trying to do.

## I am new here

Start with **[Getting Started](getting-started/installation.md)**. Install the
package, then run the **[Quick start](getting-started/quick-start.md)** to build
your first allele similarity cluster set and infer a genotype in a few lines.

## I want to understand the ideas first

Read the **Lessons**. These are short concept primers written for newcomers:

- **[Allele Similarity Clusters](lessons/allele-similarity-clusters.md)** — why
  IGHV allele assignments are ambiguous, and how clustering resolves it.
- **[Allele-based genotyping](lessons/allele-based-genotype.md)** — how PIgLET
  decides whether an allele is truly present using population-derived thresholds.

## I want the precise definitions

Go to **Core Concepts**. This is the reference layer: the
[package overview](concepts/package-overview.md), the
[`GermlineCluster` object](concepts/germline-cluster.md) that holds clustering
results, and the [built-in data](concepts/built-in-data.md) shipped with the
package.

## I have a specific task

Use the **Task Guides** (how-tos):

- **[Build allele similarity clusters](guides/build-allele-clusters.md)** — the
  full ASC workflow from a germline reference set.
- **[Infer an allele-based genotype](guides/infer-genotype.md)** — genotyping
  from AIRR-seq data.
- **[Generating ASC for a BIOMED IGHV reference](vignettes/PIgLET_biomed_asc.md)** — adapting the
  ASC pipeline to a BIOMED/OGRDB-style reference.

## I want the function signatures

Everything exported by the package is documented in the
**[API Reference](reference/index.md)**.
