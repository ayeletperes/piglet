**allele_threshold_table** - *Allele thresholds table*

Description
--------------------

A `data.table` of the allele thresholds table. The V alleles are based on the
`HVGERM` and `hv_functionality` germline reference set. The D, and the J are based on
the AIRR-C reference set (https://zenodo.org/records/10489725). The table contains these columns: allele - the IUIS allele name,
asc_allele - the allele name based on allele similarity clusters (only for V), threshold = the genotype threshold for the alleles.


Usage
--------------------
```
allele_threshold_table
```




Format
-------------------

An object of class `data.table` (inherits from `data.frame`) with 262 rows and 4 columns.


References
-------------------

Peres, et al (2022) [doi:10.1101/2022.12.26.521922](doi:10.1101/2022.12.26.521922)









