# Generate ASC for a custom reference

This guide adapts the ASC pipeline to a BIOMED-2-style V reference using the
OGRDB IGHV reference, where reads do not cover the full FWR1 region. Read the
[ASC lesson](../lessons/allele-similarity-clusters.md) first for background.

## Goal

- Infer ASC clusters for an OGRDB reference with a masked FWR1 region.
- Merge the result with PIgLET's per-allele thresholds.
- Produce an ASC-formatted reference and threshold table.

## 1. Load and mask the reference

Load the OGRDB reference and mask the FWR1 region so it matches a BIOMED-2-style
library, using `artificialFRW1Germline()`.

```r
library(data.table)
library(tigger)
library(piglet)

ref_ogrdb      <- readIgFasta("HVGERM_OGRDB.fasta")
ref_ogrdb_frw1 <- artificialFRW1Germline(ref_ogrdb)
```

## 2. Infer the allele similarity clusters

```r
asc_frw1 <- inferAlleleClusters(ref_ogrdb_frw1)

allele_table_frw1 <- setDT(asc_frw1@alleleClusterTable)[, .(iuis_allele, new_allele)]
setnames(allele_table_frw1, c("allele", "asc_allele"))
```

## 3. Merge with PIgLET thresholds

Combine the cluster table with PIgLET's threshold table, summing thresholds for
alleles that fall into the same cluster.

```r
allele_table_piglet <- fread("allele_threshold_table.tsv")
# ... assign each ASC allele the summed threshold of its member alleles ...
```

## Full tutorial

The complete, runnable tutorial is the `PIgLET_biomed_asc` vignette
(`vignettes/PIgLET_biomed_asc.Rmd` in the source tree), which also shows how to
load precomputed versions.

```r
vignette("PIgLET_biomed_asc", package = "piglet")
```
