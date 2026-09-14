**insert_gaps** - *Insert gaps into an ungapped sequence based on a gapped reference sequence.*

Description
--------------------

This function inserts gaps (e.g., `.` or `-`) into an ungapped sequence
(`ungapped`) to match the positions of gaps in a reference sequence
(`gapped`). It ensures that the aligned sequence has the same gap structure as
the reference. Vectorised over both arguments, which are paired element by element.


Usage
--------------------
```
insert_gaps(gapped, ungapped)
```

Arguments
-------------------

gapped
:   A vector of strings representing the reference sequences with gaps.

ungapped
:   A vector of strings representing the sequences without gaps.




Value
-------------------

A vector of strings with gaps inserted to match the gapped reference.



Examples
-------------------

```R
gapped <- c("caggtc..aact", "caggtc---aact")
ungapped <- c("caggtcaact", "caggtcaact")

insert_gaps(gapped, ungapped)   # "caggtc..aact", "caggtc---aact"

```


```
[1] "caggtc..aact"  "caggtc---aact"

```








