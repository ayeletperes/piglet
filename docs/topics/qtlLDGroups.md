**qtlLDGroups** - *Collapse variants carrying identical genotype information*

Description
--------------------

Two variants carry the same information when their standardised dosage vectors match
up to sign, so the mirrored case is folded in. This is exact linkage disequilibrium
rather than a correlation threshold; use `[qtlClump](qtlClump.md)` for the latter.


Usage
--------------------
```
qtlLDGroups(dosage, contig)
```

Arguments
-------------------

dosage
:   Numeric matrix of genotype dosages, variants in rows, subjects in
columns. Missing calls are mean-imputed within the variant for the comparison only.

contig
:   Character vector, one contig per row of `dosage`.




Value
-------------------

A character vector named by variant, one group key per row of `dosage`.
The key is opaque: use it to group and count, not to read. Variants with no variance
each get their own group rather than collapsing together.


Details
-------------------

The number of distinct groups is the number of independent tests, which is what a
Bonferroni threshold should be set from -- the number of variants tested overstates it,
often by half. The contig prefix stops identical patterns at different loci from
merging, so positions from different loci are never pooled.



Examples
-------------------

```R
dosage <- rbind(a      = c(0, 1, 2, 0, 1, 2),
copy   = c(0, 1, 2, 0, 1, 2),   # identical
mirror = c(2, 1, 0, 2, 1, 0),   # identical up to sign
other  = c(0, 0, 1, 2, 2, 1))
colnames(dosage) <- sprintf("s%d", 1:6)

groups <- qtlLDGroups(dosage, rep("chr2", 4))
length(unique(groups))   # three independent tests, not four

```


```
[1] 2

```








