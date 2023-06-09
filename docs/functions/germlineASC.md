# `germlineASC`

Converts IGHV germline set to ASC germline set.


## Description

Converts IGHV germline set to ASC germline set.


## Usage

```r
germlineASC(allele_cluster_table, germline)
```


## Arguments

Argument      |Description
------------- |----------------
`allele_cluster_table`     |     The allele cluster table.
`germline`     |     An IGHV germline set with matching names to the "imgt_allele" column in the allele_cluster_table.


## Value

Returns the IGHV germline set with the ASC allele names.


## Examples

```r
asc_archive <- recentAlleleClusters(doi="10.5281/zenodo.7401239", get_file = TRUE)

allele_cluster_table <- extractASCTable(archive_file = asc_archive)

data(HVGERM)

asc_germline <- germlineASC(allele_cluster_table, germline = HVGERM)
```


