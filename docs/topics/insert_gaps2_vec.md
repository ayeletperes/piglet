**insert_gaps2_vec** - *Insert gaps into an ungapped sequence based on a gapped reference sequence.*

Description
--------------------

This function inserts gaps (e.g., `.` or `-`) into an ungapped sequence (`ungapped`)
to match the positions of gaps in a reference sequence (`gapped`). It ensures that
the aligned sequence has the same gap structure as the reference.


Usage
--------------------
```
insert_gaps2_vec(gapped, ungapped, parallel = FALSE)
```

Arguments
-------------------

gapped
:   A vector of strings representing the reference sequences with gaps.

ungapped
:   A vector of strings representing the sequences without gaps.

parallel
:   A boolean flag to enable parallel processing (default: FALSE).




Value
-------------------

A vector of strings with gaps inserted to match the gapped reference.



Examples
-------------------

```R
# Example usage
gapped <- c("caggtc..aact", "caggtc---aact")
ungapped <- c("caggtcaact", "caggtcaact")

# Sequential execution
result <- insert_gaps2_vec(gapped, ungapped, parallel = FALSE)
print(result)  # "caggtc..aact", "caggtc---aact"

```


```
[1] "caggtc..aact"  "caggtc---aact"

```


```R

# Parallel execution
result_parallel <- insert_gaps2_vec(gapped, ungapped, parallel = TRUE)
print(result_parallel)

```


```
[1] "caggtc..aact"  "caggtc---aact"

```








