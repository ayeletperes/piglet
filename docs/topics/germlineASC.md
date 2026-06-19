**germlineASC** - *Converts IGHV germline set to ASC germline set.*

Description
--------------------

Converts IGHV germline set to ASC germline set.


Usage
--------------------
```
germlineASC(allele_cluster_table, germline)
```

Arguments
-------------------

allele_cluster_table
:   The allele cluster table.

germline
:   An IGHV germline set with matching names to the "iuis_allele" column in the allele_cluster_table.




Value
-------------------

Returns the IGHV germline set with the ASC allele names.



Examples
-------------------

```R
# preferably obtain the latest ASC cluster table
# asc_archive <- recentAlleleClusters(doi="10.5281/zenodo.7429773", get_file = TRUE)

# allele_cluster_table <- extractASCTable(archive_file = asc_archive)

data(HVGERM)

# example allele similarity cluster table
data(allele_cluster_table)

asc_germline <- germlineASC(allele_cluster_table, germline = HVGERM)

```

*Warning*:Column 'imgt_allele' has been renamed to 'iuis_allele'. Please update your alleleClusterTable.






