**pairingTable** - *Observed against expected gene pairing, per subject*

Description
--------------------

Builds the anchor-by-partner contingency of one subject's repertoire and compares it
to what independence of the two marginals would give. This is the phenotype the
conditional pairing scan explains.


Usage
--------------------
```
pairingTable(data, anchor, partner, subject = "subject", pseudocount = 0.5)
```

Arguments
-------------------

data
:   A `data.frame` of rearrangements, one row per sequence.

anchor
:   Column holding the gene to condition on.

partner
:   Column holding the partner gene whose distribution is the response.

subject
:   Column naming the subject. Default `"subject"`.

pseudocount
:   Added to both count and expectation before the log ratio, so a zero
cell gives a finite enrichment. Default 0.5.




Value
-------------------

A `data.table` with one row per subject, anchor gene and partner gene:
`subject`, `anchor_gene`, `partner_gene`, `count`, `depth`,
`anchor_total`, `partner_total`, `expected`, `enrichment`,
`p_anchor`, `p_partner`, `p_partner_given_anchor`,
`p_anchor_given_partner`. The `anchor` and `partner` column names are
kept on the result as attributes.


Details
-------------------

**Both marginals are divided out.** `expected` is the product of the two
margins over the subject's depth, so `enrichment` measures pairing preference
rather than how much either gene is used. That distinction is the whole analysis: in
the heavy chain the individual-level signal in the raw conditional
`p_partner_given_anchor` was found to be almost entirely marginal usage, and it
went to roughly zero once both marginals were removed. The raw conditionals are
returned for inspection, but `enrichment` is what should be scanned.

Zero cells are created only where both genes occur somewhere in that subject, so a
gene the subject does not carry stays a structural absence rather than a measured
zero. This matters: a deletion is not the same observation as a pairing the subject
never makes.

The enrichment cell is symmetric -- swapping which margin conditions leaves the number
unchanged -- so to scan the other side, call this again with `anchor` and
`partner` swapped. The omnibus test is *not* symmetric, because grouping the
cells by one segment gives each anchor a response vector of a different width, which
is a different multivariate test on a different covariance structure. Both sides are
therefore separate scans and their results must not be pooled.


Which combinations are established
-------------------


IGH D-anchored and J-anchored pairing is the analysis this was developed for and the
one that has been checked. IGK/IGL V--J, and IGH V--D and V--J, run through the same
code and are **exploratory**: nobody has shown that pairing preference in those
segments survives removing the marginals, and the heavy-chain result above is a reason
to expect trouble rather than reassurance. Treat their output as a hypothesis.



Examples
-------------------

```R
set.seed(1)
rep_dt <- data.frame(
subject = rep(sprintf("s%02d", 1:70), each = 20),
d_gene  = sample(c("D1", "D2", "D3"), 1400, replace = TRUE),
j_gene  = sample(c("J1", "J2"), 1400, replace = TRUE))

pairs <- pairingTable(rep_dt, anchor = "j_gene", partner = "d_gene")

```

*pairingTable: 1400 of 1400 rearrangements carry both a j_gene and a d_gene call.*
```R
head(pairs)

```


```
Key: <subject, anchor_gene, partner_gene>
   subject anchor_gene partner_gene count depth anchor_total partner_total
    <char>      <char>       <char> <int> <int>        <int>         <int>
1:     s01          J1           D1     4    20           11             7
2:     s01          J1           D2     3    20           11             7
3:     s01          J1           D3     4    20           11             6
4:     s01          J2           D1     3    20            9             7
5:     s01          J2           D2     4    20            9             7
6:     s01          J2           D3     2    20            9             6
   expected  enrichment p_anchor p_partner p_partner_given_anchor
      <num>       <num>    <num>     <num>                  <num>
1:     3.85  0.04890960     0.55      0.35              0.3636364
2:     3.85 -0.31366048     0.55      0.35              0.2727273
3:     3.30  0.24392558     0.55      0.30              0.3636364
4:     3.15 -0.06054154     0.45      0.35              0.3333333
5:     3.15  0.30202854     0.45      0.35              0.4444444
6:     2.70 -0.35614381     0.45      0.30              0.2222222
   p_anchor_given_partner
                    <num>
1:              0.5714286
2:              0.4285714
3:              0.6666667
4:              0.4285714
5:              0.5714286
6:              0.3333333

```


```R

# The enrichment is symmetric, so swapping the roles gives the mirrored table.
# The omnibus that follows is not, which is why the two are separate scans.
mirrored <- pairingTable(rep_dt, anchor = "d_gene", partner = "j_gene")

```

*pairingTable: 1400 of 1400 rearrangements carry both a d_gene and a j_gene call.*

See also
-------------------

`[pairingScan](pairingScan.md)`, `[pairingCellTests](pairingCellTests.md)`






