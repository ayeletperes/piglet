**assignAlleleClusters** - *Assign allele similarity clusters*

Description
--------------------

`assignAlleleClusters` uses the allele clusters annotation to change the preliminary allele
assignments to the new annotations before inferring a genotype.


Usage
--------------------
```
assignAlleleClusters(
data,
alleleClusterTable,
v_call = "v_call",
from_col = "iuis_allele",
to_col = "new_allele"
)
```

Arguments
-------------------

data
:   data.frame in AIRR format, containing V allele calls from a single subject and the sample IMGT-gapped V(D)J sequences under seq.

alleleClusterTable
:   A data.frame of the allele clusters new annotations relative to the original reference set. See details.

v_call
:   name of the V allele call column. Default is `v_call`

from_col
:   name of the column in alleleClusterTable to use as the source for the dictionary. Default is `iuis_allele`

to_col
:   name of the column in alleleClusterTable to use as the target for the dictionary. Default is `new_allele`




Value
-------------------

A modified input `data.frame` with the new assigned



Examples
-------------------

```R
# preferably obtain the latest ASC cluster table
# asc_archive <- recentAlleleClusters(doi="10.5281/zenodo.7429773", get_file = TRUE)

# allele_cluster_table <- extractASCTable(archive_file = asc_archive)

# example allele similarity cluster table
data(allele_cluster_table)

# loading TIgGER AIRR-seq b cell data
data <- tigger::AIRRDb

asc_data <- assignAlleleClusters(data, allele_cluster_table)

```

*Warning*:Column 'imgt_allele' has been renamed to 'iuis_allele'. Please update your alleleClusterTable.






