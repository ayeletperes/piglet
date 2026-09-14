**qtlScanUnivariate** - *Univariate QTL scan of every phenotype on every variant*

Description
--------------------

Ordinary least squares of each phenotype column on each variant, for the whole matrix
at once. The slope of a simple regression is a ratio of cross products, so the entire
scan is two matrix products per missingness group rather than one `lm()` call per
variant-phenotype pair.


Usage
--------------------
```
qtlScanUnivariate(pheno, dosage, min_n = 60L)
```

Arguments
-------------------

pheno
:   Numeric matrix of phenotypes, subjects in rows (row names are subject
identifiers) and phenotypes in columns (column names name the phenotype). No missing
values; drop or impute before calling.

dosage
:   Numeric matrix of genotype dosages, variants in rows (row names are
variant identifiers) and subjects in columns, in the same subject order as
`pheno`. `NA` marks no call.

min_n
:   Minimum number of complete-case subjects a missingness group needs before
it is scanned. Default 60.




Value
-------------------

A `data.table` with one row per variant-phenotype pair that produced a
finite p-value: `variant`, `phenotype`, `n` (subjects in the fit),
`beta`, `se`, `t_stat`, `p_value`. Effect direction is relative
to the coding of `dosage`, which is not necessarily the minor allele, so
`beta` must not be described in terms of a named allele without resolving
polarity first.


Details
-------------------

Variants are grouped by which subjects have no call, so every fit inside a group uses
the same complete-case subject set and `n` is exact for each row rather than an
upper bound. Groups smaller than `min_n` are skipped entirely.

A p-value from this scan must not be read without the smallest genotype class behind
the fit: the extreme tail is anti-conservative, and a non-significant result means
either no effect or no power. Use `[qtlSmallestGenotypeClass](qtlSmallestGenotypeClass.md)` on the same
dosage matrix, as `[runGeneUsageQTL](runGeneUsageQTL.md)` does. Filtering associations on the
p-value alone preferentially removes the underpowered rows, which leaves what survives
looking better powered than it is.



Examples
-------------------

```R
set.seed(1)
subjects <- sprintf("s%02d", 1:80)
dosage <- rbind(v_hit  = rep(c(0, 1, 2), length.out = 80),
v_null = rep(c(0, 0, 1, 2), length.out = 80))
colnames(dosage) <- subjects
pheno <- cbind(gene_a = 0.7 * dosage["v_hit", ] + rnorm(80),
gene_b = rnorm(80))
rownames(pheno) <- subjects

qtlScanUnivariate(pheno, dosage, min_n = 60)

```


```
   variant phenotype     n         beta        se      t_stat      p_value
    <char>    <char> <int>        <num>     <num>       <num>        <num>
1:   v_hit    gene_a    80  0.569260712 0.1236591  4.60346742 1.587254e-05
2:  v_null    gene_a    80  0.085695713 0.1365301  0.62766915 5.320540e-01
3:   v_hit    gene_b    80  0.160276866 0.1264641  1.26737086 2.087934e-01
4:  v_null    gene_b    80 -0.007009597 0.1253978 -0.05589887 9.555653e-01

```


```R

# Read a p-value with the smallest genotype class beside it, never alone.
qtlSmallestGenotypeClass(dosage)

```


```
 v_hit v_null 
    26     20 

```



See also
-------------------

`[qtlScanMultivariate](qtlScanMultivariate.md)` for a multivariate response,
`[qtlSmallestGenotypeClass](qtlSmallestGenotypeClass.md)`, `[runGeneUsageQTL](runGeneUsageQTL.md)`






