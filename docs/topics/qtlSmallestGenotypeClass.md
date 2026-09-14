**qtlSmallestGenotypeClass** - *Smallest observed genotype class per variant*

Description
--------------------

How many subjects sit in the least populated genotype class at each variant. This is
what says whether a fit is worth reading: an association driven by two homozygotes is
not the same claim as one driven by forty, and the scan's extreme tail is
anti-conservative precisely where this number is small.


Usage
--------------------
```
qtlSmallestGenotypeClass(dosage)
```

Arguments
-------------------

dosage
:   Numeric matrix of genotype dosages, variants in rows, subjects in
columns. Dosages are rounded to 0, 1 or 2.




Value
-------------------

A named numeric vector, one entry per row of `dosage`, `NA` where no
class had a subject.


Details
-------------------

Classes with no subject are not counted, so a variant seen only as 0 and 1 reports the
smaller of those two rather than zero. Computed for the whole matrix at once; index
the result by variant name rather than calling this per variant.



Examples
-------------------

```R
dosage <- rbind(all_three = c(0, 0, 0, 1, 1, 2),
two_only  = c(0, 0, 0, 0, 1, 1))
colnames(dosage) <- sprintf("s%d", 1:6)

# The empty class is not counted as zero: two_only reports 2, not 0.
qtlSmallestGenotypeClass(dosage)

```


```
all_three  two_only 
        1         2 

```








