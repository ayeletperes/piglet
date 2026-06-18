# Infer an allele-based genotype

This guide infers an individual's IGHV genotype from AIRR-seq data using
per-allele thresholds. For the reasoning, read the
[allele-based genotyping lesson](../lessons/allele-based-genotype.md).

## 1. Prepare AIRR-format input

`inferGenotypeAllele()` expects a single subject's rearrangements in AIRR format,
with at least a V allele call column (`v_call`) and, if you plan to filter to
unmutated sequences, the IMGT-gapped `sequence_alignment` column.

## 2. Get the threshold table

The per-allele thresholds come from an allele cluster table. Use the version
shipped with the package, or download the most recent one from Zenodo:

```r
library(piglet)
# Download latest from Zenodo:
archive <- recentAlleleClusters()
threshold_table <- extractASCTable(archive)
```

## 3. Run genotype inference

```r
genotype <- inferGenotypeAllele(
  data               = airr_data,
  alleleClusterTable = threshold_table,
  v_call             = "v_call",
  single_assignment  = FALSE,   # TRUE = only unambiguous calls
  find_unmutated     = FALSE,   # TRUE needs germline_db + seq column
  seq                = "sequence_alignment"
)

head(genotype)
```

Set `single_assignment = TRUE` to restrict to reads with a single V call, or
`find_unmutated = TRUE` (with `germline_db`) to restrict to unmutated reads.

## 4. ASC-aware variant

If you already have ASC annotations, `inferGenotypeAllele_asc()` performs the
same inference directly from those annotations.

## Full walkthrough

See the **[PIgLET vignette](../vignettes/PIgLET-vignette.md)** for a complete
example with real data.
