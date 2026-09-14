**annotateRepertoireIUIS** - *Add IUIS gene-group columns to a genotype-corrected repertoire*

Description
--------------------

Translates the ASC gene calls on a repertoire into the IUIS labels built by
`[ascIUISVocabulary](ascIUISVocabulary.md)`, adding `v_gene_iuis`, `j_gene_iuis` and,
where the chain has one, `d_gene_iuis`.


Usage
--------------------
```
annotateRepertoireIUIS(
data,
vocabulary,
chain,
v_gene = "v_gene",
j_gene = "j_gene",
d_gene = "d_gene",
d_call = "d_call_new"
)
```

Arguments
-------------------

data
:   A `data.frame` of rearrangements carrying the ASC gene calls.

vocabulary
:   An `ascIUISVocabulary` from `[ascIUISVocabulary](ascIUISVocabulary.md)`.

chain
:   The locus of `data`, one of `"IGH"`, `"IGK"`, `"IGL"`.

v_gene, j_gene, d_gene
:   Column names holding the ASC gene calls. Defaults
`"v_gene"`, `"j_gene"`, `"d_gene"`.

d_call
:   Column name holding the ASC *allele* call for D, used for IGH only.
Default `"d_call_new"`.




Value
-------------------

`data` as a `data.table` with the `*_gene_iuis` columns added.


Details
-------------------

For IGH the D label is per allele, so it is mapped from the allele-level call rather
than from the collapsed gene, and a call spanning more than one ASC group names no
single gene and stays `NA`. V and J map straight from their ASC gene.

Rows in and rows out are reported: a call that maps to nothing is left `NA` rather
than dropped, and the count is messaged so a labelling gap cannot pass unnoticed.



Examples
-------------------

```R
husa <- data.frame(
allele      = c("IGHV1-2*01", "IGHV1-2*02", "IGHD1-CO5H*01"),
chain       = "IGH",
gene_type   = c("IGHV", "IGHV", "IGHD"),
seq         = c("ACGTACGT", "ACGTACGA", "GGTTGGTT"),
iuis_allele = c("IGHV1-2*01", "IGHV1-2*02", "IGHD1-20*01"))
voc <- ascIUISVocabulary(husa)

rep_dt <- data.frame(
v_gene     = c("IGHV1-2", "IGHV1-2"),
d_gene     = c("IGHD1-CO5H", "IGHD1-CO5H"),
j_gene     = c("IGHJ4", "IGHJ4"),
d_call_new = c("IGHD1-CO5H*01", "IGHD1-CO5H*01"))

annotateRepertoireIUIS(rep_dt, voc, chain = "IGH")

```

*annotateRepertoireIUIS [IGH]: v_gene_iuis labelled on 2 of 2 rows (0 unlabelled).**annotateRepertoireIUIS [IGH]: j_gene_iuis labelled on 0 of 2 rows (2 unlabelled).**annotateRepertoireIUIS [IGH]: d_gene_iuis labelled on 2 of 2 rows (0 unlabelled).*
```
    v_gene     d_gene j_gene    d_call_new v_gene_iuis j_gene_iuis d_gene_iuis
    <char>     <char> <char>        <char>      <char>      <char>      <char>
1: IGHV1-2 IGHD1-CO5H  IGHJ4 IGHD1-CO5H*01        V1-2        <NA>    IGHD1-20
2: IGHV1-2 IGHD1-CO5H  IGHJ4 IGHD1-CO5H*01        V1-2        <NA>    IGHD1-20

```



See also
-------------------

`[ascIUISVocabulary](ascIUISVocabulary.md)`






