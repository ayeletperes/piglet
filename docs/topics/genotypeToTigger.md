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
head(geno_table)

```


```
       gene alleles  counts total    depth   threshold         z_score
     <char>  <char>  <char> <num>    <num>      <char>          <char>
1: IGHV1-18      01    1005  1005 4738.933       0.001         459.717
2:  IGHV1-2   02,04 664,302   966 4738.933 0.001,1e-04 302.994,438.033
3: IGHV1-24      01     105   105 4738.933       1e-04         151.847
4:  IGHV1-3      01     226   226 4738.933       1e-05        1037.957
5: IGHV1-46      01     624   624 4738.933       0.001          284.61
6: IGHV1-58   01,02   23,18    41 4738.933 1e-04,1e-04    32.724,25.46
   genotyped_alleles
              <char>
1:                01
2:             02,04
3:                01
4:                01
5:                01
6:             01,02

```



See also
-------------------

[inferGenotypeAllele](inferGenotypeAllele.md) and [inferGenotypeAllele_asc](inferGenotypeAllele_asc.md) for
producing the input genotype.






