# `zenodoArchive`

zenodoArchive


## Description

zenodoArchive
 
 zenodoArchive


## Format

[`R6Class`](#r6class) object.


## Value

Object of [`R6Class`](#r6class) for modelling an zenodoArchive for ASC cluster files


## Examples

```r
zenodo_archive <- zenodoArchive$new(
doi = "10.5281/zenodo.7401189"
)

# view available version ins the archive
archive_versions <- zenodo_archive$get_versions()

# Getting the available files in the latest zenodo archive version
files <- zenodo_archive$get_version_files()

# downloading the first file from the latest archive version
zenodo_archive$download_zenodo_files()
```


