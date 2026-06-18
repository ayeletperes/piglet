# Installation

PIgLET is an R package. It depends on several Bioconductor and AIRR-seq
ecosystem packages, so install those first.

## Bioconductor dependencies

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Biostrings", "DECIPHER"))
```

## Install PIgLET from GitHub

```r
install.packages("devtools")
devtools::install_github("ayeletperes/piglet")
```

## Build from source

To build from a local clone (this also generates the documentation):

```r
library(devtools)
install_deps()
document()
build()
install()
```

The build dependencies are:

```r
install.packages(c("devtools", "roxygen2", "testthat", "knitr", "rmarkdown"))
```

## Verify the installation

```r
library(piglet)
# The built-in human IGHV germline reference should be available:
length(HVGERM)
```

Once installed, continue to the **[Quick start](quick-start.md)**.
