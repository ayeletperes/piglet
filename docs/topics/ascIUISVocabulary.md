**ascIUISVocabulary** - *Build the ASC subgroup to IUIS gene-group vocabulary*

Description
--------------------

Derives the display labels that translate an Allele Similarity Cluster assignment into
IUIS gene names, and returns the allele set behind each label.


Usage
--------------------
```
ascIUISVocabulary(husa)
```

Arguments
-------------------

husa
:   A `data.frame` of the ASC reference (the `husa.tsv` produced by
the ASC pipeline), with columns `allele` (the ASC allele name), `chain`,
`gene_type`, `seq` and `iuis_allele`. `iuis_allele` may list
several IUIS alleles separated by `,` or `/`.




Value
-------------------

A list of class `ascIUISVocabulary`:

`labels`
:   one row per ASC subgroup: `chain`, `subgroup`,
`iuis_group`, `asc_and_iuis_group`, `iuis_group_size`.
`asc_to_iuis`
:   named character vector, ASC subgroup to IUIS label, for V
and J.
`d_allele_to_iuis`
:   named character vector, ASC *allele* to IUIS
label, for D.
`identity`
:   one row per `chain` and `iuis_label`:
`asc_genes`, `n_alleles`, `member_alleles` and
`member_seq_md5`. This is what the label means; the label itself is not.



Details
-------------------

Two ASC subgroups that share any IUIS gene must carry the same label, so subgroups and
IUIS parts are treated as a bipartite graph and every connected component gets one
merged label. Where a component holds more than one subgroup the label alone does not
identify the subgroup, so `asc_and_iuis_group` suffixes them A, B, C for display.

D is handled differently from V and J, and that asymmetry is the reason ASC labels look
inconsistent downstream. D is relabelled at IUIS gene resolution by
`[regroupASCByLabel](regroupASCByLabel.md)`: an ASC group holding alleles of more than one IUIS gene
is broken apart and groups naming the same gene are joined, so the D label is a property
of the allele rather than of the ASC gene. V and J keep the ASC grouping. This is why
the resulting labels carry the locus for D (`IGHD5-12`) but not for V and J
(`V1-18`, `J4`).

**The label is a label.** It is a display string, and it moves when the reference
set behind it moves. What identifies a group is the set of member allele sequences,
returned in `identity`; join on that, or compare it to detect that an upstream
reference changed underneath a downstream table.



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
voc$labels

```


```
    chain   subgroup iuis_group iuis_group_size asc_and_iuis_group
   <char>     <char>     <char>           <int>             <char>
1:    IGH    IGHV1-2       V1-2               1               V1-2
2:    IGH IGHD1-CO5H      D1-20               1              D1-20

```


```R

# What each label stands for. The label is a display string that moves when the
# reference set moves; these member alleles are what identify the group.
voc$identity

```


```
    chain gene_type iuis_label  asc_genes n_alleles        member_alleles
   <char>    <char>     <char>     <char>     <int>                <char>
1:    IGH      IGHV       V1-2    IGHV1-2         2 IGHV1-2*01,IGHV1-2*02
2:    IGH      IGHD   IGHD1-20 IGHD1-CO5H         1         IGHD1-CO5H*01
                     member_seq_md5
                             <char>
1: c6e558b53b4c9a94ff29b7e9c283c25f
2: 64b463743a2757689315e2263abc7353

```



See also
-------------------

`[annotateRepertoireIUIS](annotateRepertoireIUIS.md)`, `[regroupASCByLabel](regroupASCByLabel.md)`






