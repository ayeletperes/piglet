**runPairingQTL** - *Conditional gene-pairing QTL scan for one anchoring*

Description
--------------------

Runs one side of the pairing analysis: builds the enrichment table, scans every variant
against each anchor's partner distribution, sets the threshold from the number of
independent variants, and clumps the hits into leads.


Usage
--------------------
```
runPairingQTL(
data,
dosage,
variants,
anchor,
partner,
conditional,
subject = "subject",
locus = NA_character_,
pseudocount = 0.5,
min_subjects = 60L,
min_complete_fraction = 0.9,
alpha = 0.05,
r2_clump = 0.8,
min_genotype_group = 5L
)
```

Arguments
-------------------

data
:   A `data.frame` of rearrangements, one row per sequence.

dosage
:   Numeric matrix of genotype dosages, variants in rows, subjects in
columns.

variants
:   A `data.frame` describing the rows of `dosage`:
`variant`, `contig`, `pos`, `maf`. Not derived from the variant
names, which are a convention of one genotype matrix rather than a fact.

anchor, partner
:   Columns holding the gene to condition on and the partner gene.

conditional
:   A string naming what this scan estimates, e.g. `"P(J|D)"`.
Required.

subject
:   Column naming the subject. Default `"subject"`.

locus
:   Locus label written onto the output. Optional.

pseudocount
:   Added to the count, and `pseudocount * n_asc` to the total.
Default 0.5.

min_subjects
:   Minimum complete-case subjects for an anchor to be scanned, and
for a missingness group inside it. Default 60.

min_complete_fraction
:   Minimum fraction of subjects in which a partner gene must
be observed for it to stay in the response. Default 0.9.

alpha
:   Family-wise significance before the independent-variant correction.
Default 0.05.

r2_clump
:   Squared correlation above which a variant is absorbed into a lead.
Default 0.8.

min_genotype_group
:   Smallest genotype class at or above which a fit is called
well powered. Default 5.




Value
-------------------

A list of `data.table`s: `pairs`, `associations`, `leads`
and `thresholds`.


Details
-------------------

One anchoring per call, on purpose. Anchoring on each of the two segments gives two
different multivariate tests over the same variants, and a function that ran both and
returned them stacked would invite exactly the pooling that double-reports. Call it
twice, with the `conditional` each side estimates.



Examples
-------------------

```R
set.seed(1)
subjects <- sprintf("s%02d", 1:70)
rep_dt <- data.frame(
subject = rep(subjects, each = 20),
d_gene  = sample(c("D1", "D2", "D3"), 1400, replace = TRUE),
j_gene  = sample(c("J1", "J2"), 1400, replace = TRUE))

dosage <- matrix(rep(c(0, 1, 2), length.out = 70), nrow = 1,
dimnames = list("v1", subjects))
variants <- data.frame(variant = "v1", contig = "igh", pos = 1000L, maf = 0.33)

# One anchoring per call. Running both and stacking them would double-report,
# so the caller states which conditional each side estimates.
jd <- runPairingQTL(rep_dt, dosage, variants,
anchor = "j_gene", partner = "d_gene",
conditional = "P(J|D)", min_subjects = 60)

```

*pairingTable: 1400 of 1400 rearrangements carry both a j_gene and a d_gene call.**runPairingQTL [P(J|D)]: 1 variants collapse to 1 independent; threshold 0.05.**runPairingQTL [P(J|D)]: 2 anchors, 2 associations, 0 significant over 0 independent groups, 0 leads.*
```R
jd$associations

```


```
   variant conditional anchor_gene     n n_response     pillai   f_stat
    <char>      <char>      <char> <int>      <int>      <num>    <num>
1:      v1      P(J|D)          J2    69          3 0.07124716 1.662109
2:      v1      P(J|D)          J1    69          3 0.05485708 1.257556
     p_value min_genotype_group well_powered contig   pos   maf threshold
       <num>              <num>       <lgcl> <char> <int> <num>     <num>
1: 0.1838316                 22         TRUE    igh  1000  0.33      0.05
2: 0.2963209                 22         TRUE    igh  1000  0.33      0.05
   significant
        <lgcl>
1:       FALSE
2:       FALSE

```



See also
-------------------

`[pairingTable](pairingTable.md)`, `[pairingScan](pairingScan.md)`,
`[pairingCellTests](pairingCellTests.md)`






