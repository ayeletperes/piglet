**pairingScan** - *Scan variants against a gene's partner distribution*

Description
--------------------

One anchor gene at a time, tests every variant against the whole vector of partner
enrichments by MANOVA (Pillai's trace, via `[qtlScanMultivariate](qtlScanMultivariate.md)`). A
variant that merely raises or lowers how much the anchor is used shifts every partner
in the same direction and produces no signal; one that repartitions the anchor's
partners does.


Usage
--------------------
```
pairingScan(
pairs,
dosage,
conditional,
anchors = NULL,
min_subjects = 60L,
min_complete_fraction = 0.9,
min_genotype_group = 5L
)
```

Arguments
-------------------

pairs
:   A `data.table` from `[pairingTable](pairingTable.md)`.

dosage
:   Numeric matrix of genotype dosages, variants in rows, subjects in
columns.

conditional
:   A string naming what this scan estimates, e.g. `"P(J|D)"`.
Required.

anchors
:   Anchor genes to scan. Default all of them.

min_subjects
:   Minimum complete-case subjects for an anchor to be scanned, and
for a missingness group inside it. Default 60.

min_complete_fraction
:   Minimum fraction of subjects in which a partner gene must
be observed for it to stay in the response. Default 0.9.

min_genotype_group
:   Smallest genotype class at or above which a fit is called
well powered. Default 5.




Value
-------------------

A `data.table` with `conditional`, `variant`,
`anchor_gene`, `n`, `n_response`, `pillai`, `f_stat`,
`p_value`, `min_genotype_group` and `well_powered`. No significance
threshold is applied: set it from the number of independent variants
(`[qtlLDGroups](qtlLDGroups.md)`), not the number tested.


Details
-------------------

`conditional` has no default and must be stated. The two anchorings are separate
scans over the same variants, so a result that does not say which one it is cannot be
pooled, counted or plotted without double-reporting. The direction is easy to invert
by mistake and has been inverted once: at fixed partner,
<pre class = 'eq'>enrichment = \log_2 P(\mathrm{partner} \mid \mathrm{anchor}) - \log_2 P(\mathrm{partner})</pre>
and the second term is constant across the anchors, so the anchored response vector is
<code class = 'eq'>\log_2 P(\mathrm{partner} \mid \mathrm{anchor})</code> shifted by a constant. Anchoring
on J therefore gives the P(J|D) scan, and anchoring on D gives P(D|J).

`min_genotype_group` is counted on each anchor's own subject set, since complete
cases differ between anchors and a variant can be well powered for one anchor and not
another. It is returned on every row, not only on the leads: a p-value from this scan
cannot be read without it.



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

# `conditional` has no default: anchoring on J tests P(J|D), anchoring on D
# tests P(D|J), and the two are separate scans over the same variants.
pairingScan(pairs, dosage, conditional = "P(J|D)", min_subjects = 60)

```


```
   conditional variant anchor_gene     n n_response     pillai   f_stat
        <char>  <char>      <char> <int>      <int>      <num>    <num>
1:      P(J|D)      v1          J2    69          3 0.07124716 1.662109
2:      P(J|D)      v1          J1    69          3 0.05485708 1.257556
     p_value min_genotype_group well_powered
       <num>              <num>       <lgcl>
1: 0.1838316                 22         TRUE
2: 0.2963209                 22         TRUE

```



See also
-------------------

`[pairingTable](pairingTable.md)`, `[pairingCellTests](pairingCellTests.md)`






