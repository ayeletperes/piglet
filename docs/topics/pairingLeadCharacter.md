**pairingLeadCharacter** - *Characterise pairing leads: what kind of signal is this?*

Description
--------------------

A significant omnibus says a variant is associated with an anchor gene's partner
distribution. It does not say the partners were reallocated: the same p-value is
produced by a rigid shift of the whole profile, by one unusual subject, and by a
variant that simply changes how much the anchor or a partner is used. This runs the
follow-ups that separate those, per lead, and returns a verdict.


Usage
--------------------
```
pairingLeadCharacter(
leads,
pairs,
dosage,
positions = NULL,
ancestry = NULL,
replication_ancestry = "EUR",
min_replication_n = 40L,
min_genotype_group = 5L,
alpha = 0.05,
min_sep_d = 0.8,
max_overlap = 0.7,
min_mahalanobis = 1,
min_subjects = 60L,
min_complete_fraction = 0.9
)
```

Arguments
-------------------

leads
:   Significant leads for one anchoring, from `[runPairingQTL](runPairingQTL.md)`,
carrying `variant`, `anchor_gene`, `p_value` and `threshold`.

pairs
:   The `[pairingTable](pairingTable.md)` the scan used.

dosage
:   Numeric matrix of genotype dosages, variants in rows, subjects in columns.

positions
:   Optional `data.frame` of `partner_gene` and
`partner_position`. Without it `rho_position` is `NA` and the
per-partner table is left in column order.

ancestry
:   Optional `data.frame` of `subject` and `ancestry`.
Without it `ancestry_p` and `p_replication` are `NA`.

replication_ancestry
:   Ancestry group to re-test within. Default `"EUR"`.

min_replication_n
:   Minimum subjects in that group before it is attempted.
Default 40.

min_genotype_group
:   Smallest genotype class at or above which a lead is called
well powered. Default 5.

alpha
:   Nominal significance for the per-partner and margin tests. Default 0.05.

min_sep_d, max_overlap
:   A partner counts as separated at or above this Cohen's d
and at or below this overlap. Defaults 0.8 and 0.70.

min_mahalanobis
:   Joint separation a lead needs before it is called separated at
all. Default 1.

min_subjects
:   Minimum complete-case subjects for an anchor to be scanned, and
for a missingness group inside it. Default 60.

min_complete_fraction
:   Minimum fraction of subjects in which a partner gene must
be observed for it to stay in the response. Default 0.9.




Value
-------------------

A list of two `data.table`s: `leads`, the input with the summary
columns and `verdict` merged on, and `per_partner`, one row per lead and
partner with slopes and effect sizes.


What is tested
-------------------



Rigid shift against reallocation
:   The mean across partners regressed on dosage
is the shift every partner shares (`p_common`). Removing that mean and testing
what is left asks whether partners moved by *different* amounts
(`p_heterogeneous`). The second is the reallocation claim; the first is not.
One deviation column is dropped before the test because deviations from a row mean
sum to zero and the covariance would be singular.
Separation
:   `mahalanobis`, the distance between the highest and lowest
genotype class centroids in units of the pooled within-genotype covariance, after
residualising on the per-partner slopes. With `overlap_joint` beside it.
One subject
:   `p_worst_loo` refits the omnibus dropping each subject in turn
and keeps the worst p-value. A signal that rests on one person shows up here and
nowhere else, and a small genotype class is exactly when it happens.
Both margins
:   The enrichment divides both marginals out, so a lead is not
pairing-specific until neither margin moved. Logit marginal usage is regressed on
dosage for the anchor (`anchor_usage_p`) and for every partner
(`partner_usage_p_min`, Bonferroni across partners).
Ancestry
:   `ancestry_p`, whether the dosage itself is stratified, and
`p_replication`, whether the omnibus holds inside one ancestry group. A cis
variant is expected to track ancestry, so neither is a rejection on its own.
Position
:   `rho_position`, Spearman of per-partner slope against partner
position. A monotone trend along the locus is a different claim from a scattered
one. Needs `positions`.



The verdict
-------------------


`verdict` is assigned in order: `single_subject_driven`,
`overlapping_small_effect`, `pairing_plus_both_margins`,
`pairing_plus_anchor_usage`, `pairing_plus_partner_usage`, and otherwise
`pairing_only`. Its cutoffs are arguments rather than constants, and their defaults
are the ones the heavy-chain D/J analysis used. They are that analysis's policy, not
facts about the method: a different cohort or a different pair of segments should set
them deliberately rather than inherit them.



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

leads <- pairingScan(pairs, dosage, conditional = "P(J|D)", min_subjects = 60)
leads$threshold <- 0.05

out <- pairingLeadCharacter(leads, pairs, dosage, min_subjects = 60)

```

*pairingLeadCharacter: single_subject_driven 2*
```R
out$leads$verdict

```


```
[1] "single_subject_driven" "single_subject_driven"

```


```R
out$per_partner

```


```
   anchor_gene variant partner_gene partner_position beta_per_allele         se
        <char>  <char>       <char>            <num>           <num>      <num>
1:          J2      v1           D1               NA     -0.13529153 0.07376043
2:          J2      v1           D2               NA     -0.04162854 0.09274142
3:          J2      v1           D3               NA      0.06883234 0.07652916
4:          J1      v1           D1               NA      0.03381925 0.07330616
5:          J1      v1           D2               NA      0.01673175 0.07427336
6:          J1      v1           D3               NA     -0.12007040 0.07212981
       t_stat    p_value    mean_low    mean_high       sep_d       auc
        <num>      <num>       <num>        <num>       <num>     <num>
1: -1.8342019 0.07106489  0.04587312 -0.225514742 -0.49383099 0.3750000
2: -0.4488668 0.65497626 -0.17035581 -0.255356135 -0.11853933 0.4465580
3:  0.8994262 0.37164697 -0.06692589  0.071002693  0.24301314 0.6240942
4:  0.4613426 0.64604713 -0.06644226  0.001183488  0.11875510 0.5978261
5:  0.2252726 0.82245293 -0.07550623 -0.043070772  0.05617601 0.5760870
6: -1.6646432 0.10065417 -0.06844695 -0.308950695 -0.47385212 0.3505435
     overlap
       <num>
1: 0.8049736
2: 0.9527373
3: 0.9032898
4: 0.9526514
5: 0.9775920
6: 0.8127142

```



See also
-------------------

`[runPairingQTL](runPairingQTL.md)`, `[pairingCellTests](pairingCellTests.md)`






