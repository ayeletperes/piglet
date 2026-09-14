**qtlScanMultivariate** - *Multivariate QTL scan by Pillai's trace*

Description
--------------------

Scores every variant against a multivariate response, one MANOVA per variant. For a
one degree of freedom predictor Pillai's trace has a closed form,
<pre class = 'eq'>V = (u' S^{-1} u) / (g'g), \quad u = Y_c' g_c, \quad S = Y_c' Y_c</pre>
and the F approximation is exact, so a whole missingness group is scored by one matrix
product instead of a `stats::manova()` call per variant. The package tests check
this against `stats::manova()`.


Usage
--------------------
```
qtlScanMultivariate(pheno, dosage, min_n = 60L)
```

Arguments
-------------------

pheno
:   Numeric matrix of the multivariate response, subjects in rows, response
variables in columns. Complete cases only.

dosage
:   Numeric matrix of genotype dosages, variants in rows (row names are
variant identifiers) and subjects in columns, in the same subject order as
`pheno`. `NA` marks no call.

min_n
:   Minimum number of complete-case subjects a missingness group needs before
it is scanned. Default 60.




Value
-------------------

A `data.table` with `variant`, `n`, `n_response` (response
columns retained), `pillai`, `f_stat` and `p_value`.


Details
-------------------

Response columns with zero variance are dropped within each missingness group, and a
group is skipped when fewer than two response columns survive or when `n` is too
small for the covariance to be estimable (`n < p + 3`).



Examples
-------------------

```R
set.seed(1)
subjects <- sprintf("s%02d", 1:80)
dosage <- rbind(v_hit  = rep(c(0, 1, 2), length.out = 80),
v_null = rep(c(0, 0, 1, 2), length.out = 80))
colnames(dosage) <- subjects
pheno <- matrix(rnorm(80 * 4), 80, 4,
dimnames = list(subjects, paste0("partner", 1:4)))
pheno[, 1] <- pheno[, 1] + 0.8 * dosage["v_hit", ]

qtlScanMultivariate(pheno, dosage, min_n = 60)

```


```
   variant     n n_response     pillai     f_stat      p_value
    <char> <int>      <int>      <num>      <num>        <num>
1:   v_hit    80          4 0.38201126 11.5903587 2.224249e-07
2:  v_null    80          4 0.01712965  0.3267786 8.591670e-01

```



See also
-------------------

`[pairingScan](pairingScan.md)`, which uses this for the conditional pairing scan






