**runGeneUsageQTL** - *Gene-usage QTL scan for one locus*

Description
--------------------

Runs the whole usage analysis in memory: builds the phenotype, scans every variant
against every gene group, sets the significance threshold from the number of
independent variants, and clumps the significant hits into leads.


Usage
--------------------
```
runGeneUsageQTL(
data,
dosage,
variants,
segments,
positions = NULL,
locus = NA_character_,
subject = "subject",
min_subject_fraction = 0.1,
pseudocount = 0.5,
min_subjects = 60L,
alpha = 0.05,
r2_clump = 0.8,
min_genotype_group = 5L,
cis_window = 50000
)
```

Arguments
-------------------

data
:   A `data.frame` of rearrangements, one row per sequence.

dosage
:   Numeric matrix of genotype dosages, variants in rows, subjects in
columns. Subjects are matched by name against the phenotype; filter on MAF and
missingness before calling.

variants
:   A `data.frame` describing the rows of `dosage`:
`variant`, `contig`, `pos`, `maf`. Not derived from the variant
names, which are a convention of one genotype matrix rather than a fact.

segments
:   Named character vector of the columns holding each segment's gene
call, e.g. `c(V = "v_gene_iuis", D = "d_gene_iuis", J = "j_gene_iuis")`. Names
become the `segment` column.

positions
:   Optional gene coordinates for the locus, with columns `segment`,
`short` and `mid`, used for `asc_position`, `asc_span`,
`n_member`, `distance_to_asc`, `is_cis` and `nearest_gene`. Those
columns are omitted when it is `NULL`.

locus
:   Locus label written onto the output. Optional.

subject
:   Column naming the subject. Default `"subject"`.

min_subject_fraction
:   Minimum fraction of subjects a group must appear in to be
retained. Default 0.1.

pseudocount
:   Added to the count, and `pseudocount * n_asc` to the total.
Default 0.5.

min_subjects
:   Minimum complete-case subjects for a missingness group to be
scanned. Default 60.

alpha
:   Family-wise significance before the independent-variant correction.
Default 0.05.

r2_clump
:   Squared correlation above which a variant is absorbed into a lead.
Default 0.8.

min_genotype_group
:   Smallest genotype class at or above which a lead is called
well powered. Default 5.

cis_window
:   Distance within which a lead is called cis. Default 50000.




Value
-------------------

A list of `data.table`s: `phenotype`, `associations`,
`leads`, `per_asc` and `thresholds`.


Details
-------------------

The threshold is `alpha` over the number of *independent* variants after
exact-LD collapse (`[qtlLDGroups](qtlLDGroups.md)`), not over the number tested. Those differ
by about half in this data, and using the number tested is the more conservative but
wrong denominator.

Power columns are attached to the leads, where `igqtl.R` put them, and the counts
behind every step are messaged: a filter that drops variants without saying so is a bug
even when it returns cleanly.


Ceiling
-------------------


The scan is unweighted OLS on usage fractions, so its standard errors assume constant
variance while binomial sampling variance across subjects spans roughly twentyfold. Any
R-squared read off it is diluted and is an underestimate. Do not report an effect size
from this as exact.



Examples
-------------------

```R
set.seed(1)
n <- 70L
subjects <- sprintf("s%02d", seq_len(n))
dosage <- rbind(v_hit  = rep(c(0, 1, 2), length.out = n),
v_null = rep(c(0, 0, 1, 2), length.out = n))
colnames(dosage) <- subjects

genes <- c("V1-2", "V3-23", "V4-34")
rep_dt <- do.call(rbind, lapply(seq_len(n), function(i) {
w <- c(1 + dosage["v_hit", i], 1, 1)   # dosage raises V1-2 usage
data.frame(subject = subjects[i],
v_gene = sample(genes, 40, replace = TRUE, prob = w / sum(w)))
}))

variants <- data.frame(variant = rownames(dosage), contig = "igh",
pos = c(1000L, 50000L), maf = rowMeans(dosage) / 2)

res <- runGeneUsageQTL(rep_dt, dosage, variants,
segments = c(V = "v_gene"), min_subjects = 60)

```

*ascUsagePhenotype [V]: 3 of 3 groups seen in >= 10% of 70 subjects.**runGeneUsageQTL: 70 subjects, 3 groups scanned, 2 variants.**runGeneUsageQTL: 2 variants collapse to 2 independent; threshold 0.025.**runGeneUsageQTL: 6 associations, 3 significant, 3 leads over 3 groups.*
```R
res$associations

```


```
      asc variant     n         beta         se      t_stat      p_value contig
   <char>  <char> <int>        <num>      <num>       <num>        <num> <char>
1:   V1-2   v_hit    70  0.658063503 0.04280516 15.37346238 8.094034e-24    igh
2:  V4-34   v_hit    70 -0.485647308 0.05398202 -8.99646425 3.425512e-13    igh
3:  V3-23   v_hit    70 -0.372330356 0.04785119 -7.78100489 5.475203e-11    igh
4:  V3-23  v_null    70 -0.020050922 0.06515330 -0.30774992 7.592133e-01    igh
5:  V4-34  v_null    70  0.012542001 0.07915699  0.15844464 8.745760e-01    igh
6:   V1-2  v_null    70  0.006406279 0.08973909  0.07138783 9.432986e-01    igh
     pos       maf segment  locus threshold significant
   <int>     <num>  <char> <char>     <num>      <lgcl>
1:  1000 0.4928571       V   <NA>     0.025        TRUE
2:  1000 0.4928571       V   <NA>     0.025        TRUE
3:  1000 0.4928571       V   <NA>     0.025        TRUE
4: 50000 0.3642857       V   <NA>     0.025       FALSE
5: 50000 0.3642857       V   <NA>     0.025       FALSE
6: 50000 0.3642857       V   <NA>     0.025       FALSE

```


```R
res$thresholds

```


```
    locus analysis n_subjects n_variants n_independent n_asc threshold
   <char>   <char>      <int>      <int>         <int> <int>     <num>
1:   <NA>    usage         70          2             2     3     0.025
   n_significant_variants n_independent_significant
                    <int>                     <int>
1:                      1                         1

```



See also
-------------------

`[qtlScanUnivariate](qtlScanUnivariate.md)`, `[runPairingQTL](runPairingQTL.md)`






