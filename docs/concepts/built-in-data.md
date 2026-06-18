# Built-in data

PIgLET ships reference datasets so you can run both tools without downloading
anything first. Load the package and they are available by name.

| Dataset | Description |
| --- | --- |
| [`HVGERM`](../topics/HVGERM.md) | Human IGHV germline reference alleles (IMGT-gapped) |
| [`hv_functionality`](../topics/hv_functionality.md) | Functionality annotations for the IGHV alleles |

The allele cluster and threshold tables used by the genotype tool can also be
fetched in their most recent versions from Zenodo via
[`recentAlleleClusters()`](../reference/index.md#allele-based-genotype) and
[`extractASCTable()`](../reference/index.md#allele-based-genotype).

```r
library(piglet)
length(HVGERM)
head(hv_functionality)
```
