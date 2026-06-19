**plotClusterComparison** - *Compare hierarchical and Leiden clustering*

Description
--------------------

Creates a comparison visualization showing cluster assignments from both methods.


Usage
--------------------
```
plotClusterComparison(hierarchical_result, leiden_result, ...)
```

Arguments
-------------------

hierarchical_result
:   GermlineCluster object from hierarchical clustering

leiden_result
:   GermlineCluster object from Leiden clustering

...
:   Additional arguments




Value
-------------------

A ggplot object showing cluster agreement




See also
-------------------

`[inferAlleleClusters](inferAlleleClusters.md)`






