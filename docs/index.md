# PIgLET

**P**rogram for **I**nferring Immuno**g**lobulin Allele Similarity C**l**usters
and G**e**no**t**ypes.

PIgLET is an R/Rcpp package that improves genotype inference and downstream
AIRR-seq (Adaptive Immune Receptor Repertoire) analysis of the immunoglobulin
heavy chain V (IGHV) locus. It provides two main tools:

- **Allele Similarity Clusters** — reduce ambiguity among IGHV alleles caused by
  duplicated or highly similar alleles shared across genes.
- **Allele-based genotype** — determine the presence of an allele using a
  threshold derived from a naive population.

## Start here

| If you want to… | Go to |
| --- | --- |
| Get up and running fast | **[Getting Started](getting-started/installation.md)** |
| Understand the ideas | **[Lessons](lessons/allele-similarity-clusters.md)** |
| Look up precise definitions | **[Core Concepts](concepts/package-overview.md)** |
| Accomplish a specific task | **[Task Guides](guides/build-allele-clusters.md)** |
| Find a function signature | **[API Reference](reference/index.md)** |

Not sure? See **[Choose your path](learn.md)**.

## Citation

Peres et al. (2023), "Immunoglobulin Allele Similarity Clusters", *Nucleic Acids
Research*, [doi:10.1093/nar/gkad603](https://doi.org/10.1093/nar/gkad603).

## Contact

- [Ayelet Peres](mailto:ayelet.peres@yale.edu)
- [Gur Yaari](mailto:gur.yaari@yale.edu)
- [Issue tracker](https://github.com/ayeletperes/piglet/issues)

## License

PIgLET is released under
[CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/).
