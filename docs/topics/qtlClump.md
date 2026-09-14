**qtlClump** - *Independent lead variants by greedy clumping*

Description
--------------------

Takes the strongest signal, absorbs everything correlated with it above `r2`,
and repeats. What comes back is one row per independent signal rather than one row per
significant variant, which is the number to quote.


Usage
--------------------
```
qtlClump(assoc, dosage, r2 = 0.8)
```

Arguments
-------------------

assoc
:   A `data.table` of significant associations carrying at least
`variant` and `p_value`. Reordered in place.

dosage
:   Numeric matrix of genotype dosages, variants in rows, subjects in
columns, covering every variant in `assoc`. Missing calls are mean-imputed
within the variant for the correlation only.

r2
:   Squared-correlation threshold above which a variant is absorbed into an
already-taken lead.




Value
-------------------

The subset of `assoc` that are leads, in increasing p-value order.


Details
-------------------

**Ties are broken on the variant name, not on the p-value alone.** Variants in
perfect linkage disequilibrium give the same fit to the last bit that the arithmetic
carries, and their p-values then differ by around 1e-13 purely from the order of
operations in the matrix product. Sorting on the p-value alone lets that noise decide
which variant is named as the lead, and the lead variant name is what reaches a figure
and a database label. So the sort is on the p-value and then the name, which is still
an arbitrary choice between variants no data distinguishes, but the same arbitrary
choice every run. A lead should be read as one representative of its linkage group.



Examples
-------------------

```R
dosage <- rbind(lead   = c(0, 1, 2, 0, 1, 2, 1, 0),
linked = c(0, 1, 2, 0, 1, 2, 1, 0),
other  = c(0, 0, 1, 2, 2, 1, 0, 1))
colnames(dosage) <- sprintf("s%d", 1:8)

assoc <- data.table::data.table(
variant = c("lead", "linked", "other"),
p_value = c(1e-8, 2e-8, 1e-4))

# `linked` is absorbed into `lead`; two independent signals remain.
qtlClump(assoc, dosage, r2 = 0.8)

```


```
   variant p_value
    <char>   <num>
1:    lead   1e-08
2:   other   1e-04

```








