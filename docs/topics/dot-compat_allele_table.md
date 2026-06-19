**.compat_allele_table** - *Backward-compatible allele cluster table helper*

Description
--------------------

Renames the deprecated `imgt_allele` column to `iuis_allele` if present.
Called at the boundary of exported functions that accept an `alleleClusterTable`.


Usage
--------------------
```
.compat_allele_table(tbl)
```

Arguments
-------------------

tbl
:   A data.frame that may contain an `imgt_allele` column.




Value
-------------------

The table with `imgt_allele` renamed to `iuis_allele` if needed.









