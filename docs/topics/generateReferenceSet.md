**generateReferenceSet** - *Generate allele similarity reference set*

Description
--------------------

Generates the allele clusters reference set based on the clustering from [ighvClust](ighvClust.md). The function collapse
similar alleles and assign them into their respective allele clusters and family clusters. See details for naming scheme


Usage
--------------------
```
generateReferenceSet(
germline_distance,
germline_set,
alleleClusterTable,
trim_3prime_side = NULL,
family_prefix = TRUE
)
```

Arguments
-------------------

germline_distance
:   A germline set distance matrix created by [ighvDistance](ighvDistance.md).

germline_set
:   A character list of the IMGT aligned IGHV allele sequences. See details for curating options.

alleleClusterTable
:   A data.frame of the alleles and their clusters created by [ighvClust](ighvClust.md).

trim_3prime_side
:   If a 3' position trim is supplied, duplicated sequences will be checked for differential positions past the trim position. Default NULL; NULL will not activate the check. see @details

family_prefix
:   Logical. If TRUE (default), prepend "F" to the family number in ASC names (e.g. IGHVF1-G1*01). If FALSE, omit the "F" (e.g. IGHV1-G1*01).




Value
-------------------

A `list` with the re-named germline set, and a table of the allele clusters and thresholds.


Details
-------------------

Each allele is named by this scheme:
IGHVF1-G1*01 - IGH = chain, V = region, F1 = family cluster numbering (the "F" prefix can be omitted by setting family_prefix = FALSE),
G1 - allele cluster numbering, and 01 = allele numbering (given by clustering order, no connection to the expression)

In case there are alleles that are differentiated in a nucleotide position past the trimming position used for the clustering,
then the alleles are separated and are annotated with the differentiating position as so:
Say A1*01 and A1*02 are similar up to position 318, and thus collapsed in the clusters to G1*01.
Upon checking the sequences past the trim position (318), a differentiating nucleotide was seen in position 319,
A1*01 has a G, and A1*02 has a T.
Then the alleles will be separated, and the new annotation will be as so:
A1*01 = G1*01, and A1*02 = G1*01_G319T.
Where the first nucleotide indicate the base, the following number the position, and the last nucleotide the one the base changed into.









