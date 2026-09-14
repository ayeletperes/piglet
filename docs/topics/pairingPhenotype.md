**pairingPhenotype** - *The response matrix for one anchor gene*

Description
--------------------

Subjects by partner gene, of enrichment: the multivariate response
`[pairingScan](pairingScan.md)` tests each variant against. Exported because the follow-up
analyses that characterise a hit need the same matrix the scan saw, and rebuilding it
slightly differently is how a follow-up stops describing the thing it followed up.


Usage
--------------------
```
pairingPhenotype(
pairs,
anchor,
subjects,
min_complete_fraction = 0.9,
min_subjects = 60L
)
```

Arguments
-------------------

pairs
:   A `data.table` from `[pairingTable](pairingTable.md)`.

anchor
:   The anchor gene to build the response for.

subjects
:   Subjects to keep, normally `colnames(dosage)`.

min_complete_fraction
:   Minimum fraction of subjects in which a partner gene must
be observed. Default 0.9.

min_subjects
:   Minimum complete-case subjects. Default 60.




Value
-------------------

A numeric matrix, subjects in rows (named) and partner genes in columns, or
`NULL`.


Details
-------------------

Partner genes observed in too few subjects are dropped *before* complete cases are
taken, so one rare partner cannot cost the anchor its whole subject set. Returns
`NULL` when fewer than two partners or fewer than `min_subjects` subjects
survive, which is the signal to skip this anchor rather than an error.



Examples
-------------------

```R
set.seed(1)
subjects <- sprintf("s%02d", 1:70)
rep_dt <- data.frame(
subject = rep(subjects, each = 20),
d_gene  = sample(c("D1", "D2", "D3"), 1400, replace = TRUE),
j_gene  = sample(c("J1", "J2"), 1400, replace = TRUE))
pairs <- pairingTable(rep_dt, anchor = "j_gene", partner = "d_gene")

```

*pairingTable: 1400 of 1400 rearrangements carry both a j_gene and a d_gene call.*
```R

# Subjects by partner gene, of enrichment: what one anchor is scanned on.
pheno <- pairingPhenotype(pairs, "J1", subjects, min_subjects = 60)
dim(pheno)

```


```
[1] 69  3

```



See also
-------------------

`[pairingScan](pairingScan.md)`






