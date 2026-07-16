**genotypeToTigger** - *Convert a PIgLET genotype to a TIgGER/VDJbase genotype table*

Description
--------------------

`genotypeToTigger` reshapes the output of [inferGenotypeAllele](inferGenotypeAllele.md) or
[inferGenotypeAllele_asc](inferGenotypeAllele_asc.md) into the per-gene table layout used by TIgGER
and VDJbase. The set of genotyped alleles is determined from the allele
z-score and a z-score threshold.


Usage
--------------------
```
genotypeToTigger(
genotype,
level = c("gene", "asc"),
z_threshold = 0,
file = NULL
)
```

Arguments
-------------------

genotype
:   a genotype `data.frame` produced by
[inferGenotypeAllele](inferGenotypeAllele.md) (one row per allele, with a
`z_score` column) or by [inferGenotypeAllele_asc](inferGenotypeAllele_asc.md)
(one row per gene, with a `genotype_confidence`
column). The input type is detected automatically.

level
:   the row key for the output table. `"gene"` (default)
keys by V gene; `"asc"` keys by the allele similarity
cluster. `"asc"` requires the genotype to carry ASC
information (i.e. [inferGenotypeAllele_asc](inferGenotypeAllele_asc.md) output, or
[inferGenotypeAllele](inferGenotypeAllele.md) run with `asc_annotation`
or `translate_to_asc`).

z_threshold
:   the z-score threshold for calling an allele present in the
genotype. Alleles with `z_score >= z_threshold` are
listed in `genotyped_alleles`. Default is 0.

file
:   optional path. When supplied, the table is also written as
a tab-separated file with `data.table::fwrite`.




Value
-------------------

A `data.table` with one row per gene (or ASC cluster), and columns:

+  `gene` - the V gene or ASC cluster.
+  `alleles` - the candidate allele numbers, comma-separated and
ordered by descending count.
+  `counts` - the read counts, comma-separated, matching `alleles`.
+  `total` - the total read count for the gene.
+  `depth` - the per-locus repertoire depth used as the denominator
of the z-score.
+  `threshold` - the per-allele presence thresholds, comma-separated,
matching `alleles`.
+  `z_score` - the per-allele z-scores, comma-separated, matching `alleles`.
+  `genotyped_alleles` - the allele numbers with `z_score >= z_threshold`,
comma-separated.




Examples
-------------------

```R
# loading TIgGER AIRR-seq b cell data
data <- tigger::AIRRDb

data(allele_threshold_table)
data(HVGERM)

genotype <- inferGenotypeAllele(
data = data,
allele_threshold_table = allele_threshold_table,
germline_db = HVGERM, find_unmutated = TRUE)

# convert to the TIgGER/VDJbase table layout
geno_table <- genotypeToTigger(genotype)

```

**Error in genotypeToTigger(genotype)**: could not find function "genotypeToTigger"
```R
head(geno_table)

```

**Error**: object 'geno_table' not found

See also
-------------------

[inferGenotypeAllele](inferGenotypeAllele.md) and [inferGenotypeAllele_asc](inferGenotypeAllele_asc.md) for
producing the input genotype.






