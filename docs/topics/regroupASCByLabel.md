**regroupASCByLabel** - *Regroup Allele Similarity Clusters by a per-allele label*

Description
--------------------

Post-processes an ASC assignment (e.g. the `alleleClusterTable` produced by
`[inferAlleleClusters](inferAlleleClusters.md)`) by re-grouping alleles according to a per-allele
label -- typically the IMGT gene held in `imgt_allele`. It corrects the two
ways an ASC grouping and a gene-level view disagree: an ASC group that merges
several distinct genes, and a gene whose alleles are split across ASC groups. The
clustering itself is not re-run; only the output is relabelled, so this drops onto
the end of the pipeline without changing anything upstream.


Usage
--------------------
```
regroupASCByLabel(
x,
group_col = "allele_cluster",
label_col = "imgt_allele",
duplicated_col = "duplicated",
action = c("both", "split", "merge")
)
```

Arguments
-------------------

x
:   A `data.frame` / `AlleleClusterTable` carrying at least the ASC
group column and a label column.

group_col
:   Name of the ASC group column. Default `"allele_cluster"`.
An ASC-allele column such as `"asc_allele"` works too: an allele suffix
(`*NN`) is stripped, so `IGHD1-CO5H*01` collapses to the ASC gene
`IGHD1-CO5H`.

label_col
:   Name of the label column (e.g. `"imgt_allele"` or
`"iuis_allele"`). A single value names one gene; a value listing several
genes separated by `,` or `/` (a shared / ambiguous allele) names them
all and will force those genes together.

duplicated_col
:   Optional column naming *other* genes that share the
same allele, separated by `,` or `/`. Default `"duplicated"`;
ignored if the column is absent or is a logical flag (which carries no gene
names). Not needed when the shared genes are already listed in `label_col`.

action
:   One of:

`"both"`
:   (default) re-group so each label becomes one bin; a shared
allele (a label present under two genes) forces those genes together.
`"split"`
:   only break ASC groups that mix labels; a label split
across two ASC groups stays split.
`"merge"`
:   only merge ASC groups that share a label or a shared
allele; groups that mix labels are left intact.





Value
-------------------

`x` with an added character column `regrouped`. A group made of
several genes (a merge, e.g. a shared allele) is tagged with all of them joined
by `/`, e.g. `"IGHD2-2/IGHD2-8"`, so no gene name is dropped.


Details
-------------------

The label is not produced by the clustering; supply it on the table (the IMGT gene
in `imgt_allele`, and any other genes that share the allele in
`duplicated`). A gene is taken as the part of a label before `*`, so
allele suffixes and mutation tags (`IGHD1-20*01_a7g`) are ignored.

Genes are grouped within a locus. A gene name is unique only inside its locus, so a
label stripped of its prefix (`"J1"` rather than `"IGHJ1"`) would otherwise
put IGHJ1, IGKJ1 and IGLJ1 in one bin. The locus is read from the label where it still
carries one and from `group_col` otherwise; if neither supplies it, the grouping
falls back to plain gene names and warns. `action = "split"` is unaffected, since
its output is prefixed by the ASC group.



Examples
-------------------

```R
tb <- data.frame(
allele         = paste0("a", 1:5),
allele_cluster = c("G1", "G1", "G2", "G3", "G3"),
imgt_allele    = c("IGHD1-1*01", "IGHD1-2*01", "IGHD1-1*02",
"IGHD2-2*01", "IGHD2-8*01"),
duplicated     = c("", "", "", "IGHD2-8*01", ""),   # a4 shared: 2-2 and 2-8
stringsAsFactors = FALSE)
regroupASCByLabel(tb, action = "both")$regrouped

```


```
[1] "IGHD1-1"         "IGHD1-2"         "IGHD1-1"         "IGHD2-2/IGHD2-8"
[5] "IGHD2-2/IGHD2-8"

```



See also
-------------------

`[inferAlleleClusters](inferAlleleClusters.md)`






