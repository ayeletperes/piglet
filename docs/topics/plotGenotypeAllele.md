**plotGenotypeAllele** - *Plot a PIgLET genotype with a per-allele confidence panel*

Description
--------------------

`plotGenotypeAllele` draws a genotype bar plot in the style of
`tigger::plotGenotype`, with an aligned confidence panel. Unlike TIgGER,
which carries a single confidence value per gene, PIgLET assigns a z-score to
each allele; the confidence panel therefore mirrors the allele bars and colors
each allele segment by its own z-score.


Usage
--------------------
```
plotGenotypeAllele(
genotype,
level = c("gene", "asc"),
z_threshold = 0,
gene_sort = c("name", "position"),
text_size = 12,
silent = FALSE,
...
)
```

Arguments
-------------------

genotype
:   a genotype table from [genotypeToTigger](genotypeToTigger.md), or a raw
genotype from [inferGenotypeAllele](inferGenotypeAllele.md) /
[inferGenotypeAllele_asc](inferGenotypeAllele_asc.md) (converted automatically via
[genotypeToTigger](genotypeToTigger.md)).

level
:   row key when converting a raw genotype: `"gene"`
(default) or `"asc"`. Ignored when `genotype` is
already a converted table.

z_threshold
:   z-score threshold used to select the genotyped alleles when
converting a raw genotype. Default 0.

gene_sort
:   gene ordering, passed to `alakazam::sortGenes`:
`"name"` (default) or `"position"`.

text_size
:   point size of the plotted text.

silent
:   if `TRUE` return the grob without drawing it.

...
:   additional arguments passed to `ggplot2::theme`.




Value
-------------------

A `gridExtra` grob combining the allele panel and the per-allele
z-score confidence panel (returned invisibly).


Details
-------------------

Only the genotyped alleles (those with `z_score >= z_threshold`)
are shown. Each gene/ASC row is split into equal segments, one per genotyped
allele; the left panel colors the segments by allele identity and the right
panel colors the same segments by their z-score on a continuous blue scale.



Examples
-------------------

```R
### Not run:
data <- tigger::AIRRDb
# data(allele_threshold_table)
# data(HVGERM)
# genotype <- inferGenotypeAllele(data, allele_threshold_table = allele_threshold_table,
# germline_db = HVGERM, find_unmutated = TRUE)
# plotGenotypeAllele(genotype)

```



See also
-------------------

[genotypeToTigger](genotypeToTigger.md), [inferGenotypeAllele](inferGenotypeAllele.md),
[inferGenotypeAllele_asc](inferGenotypeAllele_asc.md)






