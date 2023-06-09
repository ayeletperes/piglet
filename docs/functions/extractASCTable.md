# `extractASCTable`

Extracts the allele cluster table from the archive file.


## Description

Extracts the allele cluster table from the archive file.


## Usage

```r
extractASCTable(archive_file = NULL)
```


## Arguments

Argument      |Description
------------- |----------------
`archive_file`     |     A path to the asc archive file. Default is null. (see details)


## Details

For downloading the latest archive file with the updated allele cluster table, use the function `recentAlleleClusters` .


## Value

Returns the allele cluster table.
 
 The table columns:
 `new_allele` - the ASC given allele name
 `func_group` - the ASC cluster number
 `imgt_allele` - the original IUIS/IMGT allele name
 `thresh` - the allele threshold for ASC-based genotype inference
 `amplicon_length` - is the original length of the reference set.


## Examples

```r
asc_archive <- recentAlleleClusters(doi="10.5281/zenodo.7401239", get_file = TRUE)

allele_cluster_table <- extractASCTable(archive_file = asc_archive)
```


