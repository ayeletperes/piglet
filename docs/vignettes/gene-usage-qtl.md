# Gene-usage QTL: from repertoire to association



# What this analysis asks

Two questions, over the same variants:

**Usage.** Does a germline variant change *how much* a gene group is used? The phenotype
is the group's share of its own segment's repertoire in one subject, and each variant is
regressed on it.

**Conditional pairing.** Does a variant change *which partner* a gene pairs with, rather
than only how often it is used? The phenotype is the anchor-by-partner enrichment, with
both marginals divided out, and each variant is tested against the whole partner vector
at once.

They are different claims. A variant that raises a D segment's usage shifts every J
partner in the same direction and produces no pairing signal; one that repartitions the
segment's partners does.

# Input

Three things go in, and each is a pinned version rather than a file:

| Input | What it is | Pin |
|---|---|---|
| Repertoire | One row per rearrangement, with the gene call for each segment, genotype-corrected so a call the subject does not carry is already gone | the run that produced it |
| Genotype dosages | Variants by subjects, 0/1/2, `NA` for no call | the VCF release and its filters |
| Gene coordinates | Start and end per gene, used for distance and *cis* calls | the annotation release, by its release name, never a `current` symlink |

The package takes these as matrices and tables, not paths. Reading them, and recording
which release each came from, belongs with whatever runs the analysis, next to the output
it produces. What matters here is that the pin travels: an annotation release that gains
or loses genes moves every `nearest_gene`, `distance_to_asc` and `is_cis` without moving
a single coordinate.

## Labels are not identities

The gene groups the scan uses are labels. `ascIUISVocabulary()` builds them from an ASC
reference and, alongside them, returns what each one stands for:


``` r
voc <- ascIUISVocabulary(husa)           # the ASC reference table
rep_dt <- annotateRepertoireIUIS(rep_dt, voc, chain = "IGH")

head(voc$identity)                       # iuis_label, member_alleles, member_seq_md5
```

`identity` is the part to keep beside any downstream table. The label is a display string
that moves when the reference set moves; the member alleles and the hash of their
sequences are what say whether two tables are talking about the same thing. An allele can
also be redefined upstream while keeping its name, which is the change a list of names
cannot show and the hash can.

Note that D is labelled at IUIS gene resolution while V and J keep the ASC grouping, so D
labels carry the locus (`IGHD5-12`) and V and J do not (`V1-18`, `J4`). That asymmetry is
a property of the labelling, not a naming slip.

# Analysis

A small synthetic example, so this vignette runs anywhere. Real cohorts are hundreds of
subjects and thousands of variants; the shapes are the same.


``` r
set.seed(1)
n_sub <- 80L
subjects <- sprintf("S%03d", seq_len(n_sub))

# one variant that shifts IGHV1-2 usage, one that does nothing
dosage <- rbind(v_hit  = rep(c(0, 1, 2), length.out = n_sub),
                v_null = sample(c(0, 1, 2), n_sub, replace = TRUE))
colnames(dosage) <- subjects

genes <- c("V1-2", "V3-23", "V4-34")
rep_dt <- rbindlist(lapply(seq_len(n_sub), function(i) {
  w <- c(1 + dosage["v_hit", i], 1, 1)          # the effect: dosage raises V1-2
  data.table(subject = subjects[i],
             v_gene_iuis = sample(genes, 200, replace = TRUE, prob = w / sum(w)),
             j_gene_iuis = sample(c("J4", "J6"), 200, replace = TRUE))
}))

variants <- data.table(variant = rownames(dosage), contig = "igh",
                       pos = c(1000L, 50000L), maf = rowMeans(dosage) / 2)
```

## The usage scan


``` r
res <- runGeneUsageQTL(
  rep_dt, dosage, variants,
  segments = c(V = "v_gene_iuis", J = "j_gene_iuis"),
  locus = "IGH", min_subjects = 60L)
#> ascUsagePhenotype [V]: 3 of 3 groups seen in >= 10% of 80 subjects.
#> ascUsagePhenotype [J]: 2 of 2 groups seen in >= 10% of 80 subjects.
#> runGeneUsageQTL: 80 subjects, 5 groups scanned, 2 variants.
#> runGeneUsageQTL: 2 variants collapse to 2 independent; threshold 0.025.
#> runGeneUsageQTL: 10 associations, 3 significant, 3 leads over 3 groups.

res$associations[, .(variant, asc, n, beta, p_value, significant)]
#>     variant    asc     n          beta      p_value significant
#>      <char> <char> <int>         <num>        <num>      <lgcl>
#>  1:   v_hit   V1-2    80  0.5394273290 1.691251e-37        TRUE
#>  2:   v_hit  V3-23    80 -0.3547902229 2.224555e-27        TRUE
#>  3:   v_hit  V4-34    80 -0.3257096560 4.635064e-24        TRUE
#>  4:  v_null  V3-23    80 -0.0613359689 1.951152e-01       FALSE
#>  5:  v_null   V1-2    80  0.0439928978 5.182383e-01       FALSE
#>  6:  v_null     J4    80 -0.0046017949 8.477643e-01       FALSE
#>  7:  v_null     J6    80  0.0046017949 8.477643e-01       FALSE
#>  8:   v_hit     J6    80 -0.0040747718 8.592487e-01       FALSE
#>  9:   v_hit     J4    80  0.0040747718 8.592487e-01       FALSE
#> 10:  v_null  V4-34    80  0.0004763526 9.915705e-01       FALSE
```

Every step says how many rows went in and how many came out. A filter that drops variants
without saying so is a bug even when it returns cleanly.

The threshold is `alpha` over the number of **independent** variants, not the number
tested. `qtlLDGroups()` collapses variants whose standardised dosage vectors match up to
sign, which in real IGH data is about half of them:


``` r
res$thresholds
#>     locus analysis n_subjects n_variants n_independent n_asc threshold
#>    <char>   <char>      <int>      <int>         <int> <int>     <num>
#> 1:    IGH    usage         80          2             2     5     0.025
#>    n_significant_variants n_independent_significant
#>                     <int>                     <int>
#> 1:                      1                         1
```

## The conditional pairing scan

`conditional` has no default and must be stated. The two anchorings are separate tests
over the same variants, so a result that does not say which one it is gets double-counted
the moment two tables are stacked. Anchoring on J gives the P(J|D) scan; anchoring on D
gives P(D|J).


``` r
pairs <- pairingTable(rep_dt, anchor = "j_gene_iuis", partner = "v_gene_iuis")
#> pairingTable: 16000 of 16000 rearrangements carry both a j_gene_iuis and a v_gene_iuis call.
pairs[1:3, .(subject, anchor_gene, partner_gene, count, expected, enrichment)]
#> Key: <subject, anchor_gene, partner_gene>
#>    subject anchor_gene partner_gene count expected  enrichment
#>     <char>      <char>       <char> <int>    <num>       <num>
#> 1:    S001          J4         V1-2    32    30.25  0.07985331
#> 2:    S001          J4        V3-23    40    40.70 -0.02472243
#> 3:    S001          J4        V4-34    38    39.05 -0.03881925
```

`enrichment` is `log2((count + pc) / (expected + pc))` where `expected` is the product of
the two margins over the subject's depth. **Both marginals are divided out**, which is
what separates a pairing preference from plain usage. The raw conditionals
(`p_partner_given_anchor`) are returned too, but they should not be scanned: in the heavy
chain the individual-level signal in the raw conditional turned out to be almost entirely
marginal usage and went to roughly zero once both margins were removed.


``` r
jd <- runPairingQTL(rep_dt, dosage, variants,
                    anchor = "j_gene_iuis", partner = "d_gene_iuis",
                    conditional = "P(J|D)")

dj <- runPairingQTL(rep_dt, dosage, variants,
                    anchor = "d_gene_iuis", partner = "j_gene_iuis",
                    conditional = "P(D|J)")

cells <- pairingCellTests(jd$pairs, dosage,
                          omnibus = rbind(jd$associations, dj$associations),
                          conditional = "P(J|D)")
```

`pairingCellTests()` is the follow-up inside a row the omnibus already cleared: which
partners moved, in which direction, and how many subjects sat either side of the split.
`marked` and `marked_strict` are computed once here and carried, so a table and a figure
make the same call rather than each re-deriving it.

### Which combinations are established

IGH D-anchored and J-anchored pairing is the analysis this was developed for. IGK and IGL
V--J, and IGH V--D and V--J, run through the same code and are **exploratory**: nobody has
shown that pairing preference in those segments survives removing the marginals, and the
heavy-chain result above is a reason to expect trouble there rather than reassurance.
Treat their output as a hypothesis, not a finding.

# Output

`runGeneUsageQTL()` returns five tables:

| Table | One row per | Holds |
|---|---|---|
| `phenotype` | subject, group | `count`, `total`, `n_asc`, `usage`, `logit_usage` |
| `associations` | variant, group | `beta`, `se`, `t_stat`, `p_value`, `significant` |
| `leads` | independent signal | the association plus `min_genotype_group`, `well_powered`, `is_cis`, `nearest_gene` |
| `per_asc` | group | `min_p`, `n_significant`, `n_significant_groups`, `median_usage` |
| `thresholds` | analysis | `n_subjects`, `n_variants`, `n_independent`, `threshold`, the significant counts |

`runPairingQTL()` returns `pairs`, `associations`, `leads` and `thresholds` in the same
shape, with `conditional` on every row.

## Reading a p-value

A p-value from either scan must be read with `min_genotype_group`, the smallest genotype
class behind the fit. The extreme tail is anti-conservative, and a non-significant result
means either no effect or no power. That also makes filtering associations on the p-value
alone a mistake: it preferentially removes the underpowered rows, so what survives looks
better powered than it is.

`n_asc` is part of the pseudocount, so `usage` cannot be recomputed from `count` and
`total` without it. And `beta` is signed relative to the coding of the dosage matrix,
which is not necessarily the minor allele, so it must not be described in terms of a named
allele without resolving polarity first.

## One ceiling worth stating

The scan is unweighted OLS on usage fractions. Its standard errors assume constant
variance, while binomial sampling variance across subjects spans roughly twentyfold, so
any R-squared read off it is diluted and is an underestimate. Effect sizes from this
should not be reported as exact. A weighted or beta-binomial fit is the upgrade if that
ever matters.
