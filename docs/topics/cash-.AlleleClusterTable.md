**$.AlleleClusterTable** - *Dollar accessor for AlleleClusterTable (backward-compatible deprecation)*

Description
--------------------

Intercepts access to the deprecated `imgt_allele` column name and
redirects it to `iuis_allele` with a deprecation warning.


Usage
--------------------
```
"$"(x, name)
```

Arguments
-------------------

x
:   An AlleleClusterTable object.

name
:   Column name.











