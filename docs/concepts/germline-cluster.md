# The `GermlineCluster` object

`inferAlleleClusters()` returns an S4 object of class `GermlineCluster`. It is the
structured container for everything the ASC workflow produces:

- the **allele cluster table** (the renamed alleles and their cluster
  assignments),
- the **hierarchical clustering** of the germline set,
- the **threshold parameters** used (family and allele cluster thresholds), and
- the family cut data.

It provides `print`, `summary`, and `plot` methods. Calling `plot()` on the
object draws the allele cluster dendrogram.

```r
clusters <- inferAlleleClusters(HVGERM)
summary(clusters)
plot(clusters)
```

See the full slot-by-slot reference:

- **[`GermlineCluster-class`](../topics/GermlineCluster-class.md)** — all slots
  and methods.
- **[`new_germline_cluster`](../reference/index.md#allele-similarity-cluster)** —
  low-level constructor.
