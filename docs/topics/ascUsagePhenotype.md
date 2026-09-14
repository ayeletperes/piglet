**ascUsagePhenotype** - *Per-subject usage of each gene group*

Description
--------------------

A group's share of its own segment's repertoire in each subject, with a pseudocount so
the logit is finite. This is the phenotype the usage scan explains.


Usage
--------------------
```
ascUsagePhenotype(
data,
segments,
subject = "subject",
min_subject_fraction = 0.1,
pseudocount = 0.5
)
```

Arguments
-------------------

data
:   A `data.frame` of rearrangements, one row per sequence.

segments
:   Named character vector of the columns holding each segment's gene
call, e.g. `c(V = "v_gene_iuis", D = "d_gene_iuis", J = "j_gene_iuis")`. Names
become the `segment` column.

subject
:   Column naming the subject. Default `"subject"`.

min_subject_fraction
:   Minimum fraction of subjects a group must appear in to be
retained. Default 0.1.

pseudocount
:   Added to the count, and `pseudocount * n_asc` to the total.
Default 0.5.




Value
-------------------

A `data.table` with `subject`, `asc`, `count`,
`total`, `segment`, `n_asc`, `usage` and `logit_usage`.


Details
-------------------

Zero counts are kept deliberately. A variant that deletes a gene drives its usage to
zero, and that is the canonical signal in this family, so dropping zeros would drop the
effect. Every subject gets a row for every retained group, whether or not they used it.

Groups seen in fewer than `min_subject_fraction` of subjects are dropped before
the shares are taken, and the count kept is reported. `n_asc` is how many groups
the share was taken over and is part of the pseudocount, so `usage` cannot be
recomputed from `count` and `total` without it.



Examples
-------------------

```R
set.seed(1)
rep_dt <- data.frame(
subject = rep(sprintf("s%02d", 1:70), each = 20),
v_gene  = sample(c("V1-2", "V3-23", "V4-34"), 1400, replace = TRUE),
j_gene  = sample(c("J4", "J6"), 1400, replace = TRUE))

usage <- ascUsagePhenotype(rep_dt, segments = c(V = "v_gene", J = "j_gene"))

```

*ascUsagePhenotype [V]: 3 of 3 groups seen in >= 10% of 70 subjects.**ascUsagePhenotype [J]: 2 of 2 groups seen in >= 10% of 70 subjects.*
```R
head(usage)

```


```
   subject    asc count total segment n_asc     usage logit_usage
    <char> <char> <int> <int>  <char> <int>     <num>       <num>
1:     s01   V1-2     7    20       V     3 0.3488372  -0.6241543
2:     s01  V3-23     7    20       V     3 0.3488372  -0.6241543
3:     s01  V4-34     6    20       V     3 0.3023256  -0.8362480
4:     s02   V1-2     8    20       V     3 0.3953488  -0.4248832
5:     s02  V3-23     8    20       V     3 0.3953488  -0.4248832
6:     s02  V4-34     4    20       V     3 0.2093023  -1.3291359

```


```R

# Zero counts are kept on purpose: a variant that deletes a gene drives its
# usage to zero, and that is the signal rather than missing data.
sum(usage$count == 0)

```


```
[1] 1

```



See also
-------------------

`[runGeneUsageQTL](runGeneUsageQTL.md)`






