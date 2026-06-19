**[[.AlleleClusterTable** - *Double-bracket accessor for AlleleClusterTable (backward-compatible deprecation)*

Description
--------------------

Intercepts access to the deprecated `imgt_allele` column name and
redirects it to `iuis_allele` with a deprecation warning.


Usage
--------------------
```
"[["(x, i, ...)
```

Arguments
-------------------

x
:   An AlleleClusterTable object.

i
:   Index or column name.

...
:   Additional arguments.











