**pairingCellTests** - *Per-cell follow-up under a significant pairing hit*

Description
--------------------

Regresses each single anchor-partner cell on the variant, so a significant omnibus can
be read as which partners moved and in which direction. Returns the fitted slope
alongside the raw phenotype means in the lowest and highest observed genotype class,
and how many subjects sat either side, which is what says whether a difference in
means is worth anything.


Usage
--------------------
```
pairingCellTests(
pairs,
dosage,
omnibus,
conditional,
alpha = 0.05,
min_subjects = 60L,
min_complete_fraction = 0.9
)
```

Arguments
-------------------

pairs
:   A `data.table` from `[pairingTable](pairingTable.md)`.

dosage
:   Numeric matrix of genotype dosages, variants in rows, subjects in
columns.

omnibus
:   A `[pairingScan](pairingScan.md)` result carrying a logical
`significant` column. It may hold both anchorings stacked together, and
normally should: a variant that cleared the omnibus on either side is worth reading
cells for on both, and the per-row `omnibus_p_value` and
`omnibus_significant` are matched on `variant` and `anchor_gene`,
so the other side's rows simply do not match this side's anchors.

conditional
:   A string naming what this scan estimates, e.g. `"P(J|D)"`.
Required.

alpha
:   Nominal significance for a cell inside a cleared row. Default 0.05.

min_subjects
:   Minimum complete-case subjects for an anchor to be scanned, and
for a missingness group inside it. Default 60.

min_complete_fraction
:   Minimum fraction of subjects in which a partner gene must
be observed for it to stay in the response. Default 0.9.




Value
-------------------

A `data.table` with `conditional`, `variant`,
`anchor_gene`, `partner_gene`, `n`, `beta`, `p_value`,
`n_low`, `n_high`, `mean_low`, `mean_high`, `delta_mean`,
`omnibus_p_value`, `omnibus_significant`, `min_genotype_group`,
`marked` and `marked_strict`.


Details
-------------------

A cell is only meaningful inside a row the omnibus already called significant, so the
scan is restricted to those variants and `omnibus_significant` travels with every
row. `marked` applies nominal alpha inside a cleared row; `marked_strict`
divides alpha by the row width, the family being that row's partner genes. Both are
returned so a downstream table and a figure make the same call rather than each
re-deriving it.



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
dosage <- matrix(rep(c(0, 1, 2), length.out = 70), nrow = 1,
dimnames = list("v1", subjects))

omnibus <- pairingScan(pairs, dosage, conditional = "P(J|D)", min_subjects = 60)
omnibus$significant <- omnibus$p_value < 0.05

pairingCellTests(pairs, dosage, omnibus, conditional = "P(J|D)", min_subjects = 60)

```

*pairingCellTests: 0 variants with at least one significant P(J|D) row.*
```
Null data.table (0 rows and 0 cols)

```



See also
-------------------

`[pairingScan](pairingScan.md)`






