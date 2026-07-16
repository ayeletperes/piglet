**allelePalette** - *Allele color palette*

Description
--------------------

`allelePalette` returns a named vector of colors for a set of alleles,
reusing the VDJbase allele color scheme (see the `vdjbasevis` package) for
the curated allele numbers and generating additional, visually distinct colors
when there are more alleles than the curated palette covers. This lets genotype
plots scale to an arbitrary number of alleles while staying consistent with
VDJbase for the common ones.


Usage
--------------------
```
allelePalette(alleles)
```

Arguments
-------------------

alleles
:   a vector of allele identifiers, typically allele numbers such as
`"01"`, `"02"`. A novel-allele suffix after `"_"` shares the
color of its base allele.




Value
-------------------

a named character vector mapping each unique allele to a hex color.


Details
-------------------

The curated colors are taken from `vdjbasevis::allelePalette`
(alleles `01`-`19`). Alleles beyond that first draw from a pool of
additional distinct colors and then, if still short, from evenly spaced hues
generated in HCL space.



Examples
-------------------

```R
allelePalette(c("01", "02", "06"))

```

**Error in allelePalette(c("01", "02", "06"))**: could not find function "allelePalette"
```R
# also covers many alleles without running out of colors
length(allelePalette(sprintf("%02d", 1:40)))

```

**Error in allelePalette(sprintf("%02d", 1:40))**: could not find function "allelePalette"






