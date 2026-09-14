**insert_gaps2_vec** - *Deprecated: use `[insert_gaps](insert_gaps.md)`*

Description
--------------------

Renamed. The `2` marked a second attempt rather than anything about behaviour,
and `_vec` restated that it is vectorised. Results are unchanged.


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
:   Ignored. Accepted so existing calls keep working.




Value
-------------------

As `[insert_gaps](insert_gaps.md)`.




See also
-------------------

`[insert_gaps](insert_gaps.md)`






