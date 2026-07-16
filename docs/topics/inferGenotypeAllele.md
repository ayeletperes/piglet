**inferGenotypeAllele** - *Allele based genotype inference*

Description
--------------------

`inferGenotypeAllele` infer an individual's genotype based on the allele-base method.
The method utilize the allele specific threshold to determine the presence of an allele in the genotype.
More specifically, based on the allele frequency, repertoire depth, and the specific allele threshold, a confidence level (Z score) is calculated
for the presence of the allele in the genotype. The user can select the confidence level for the genotype inference.


Usage
--------------------
```
inferGenotypeAllele(
data,
allele_threshold_table = NULL,
call = "v_call",
asc_annotation = FALSE,
single_assignment = FALSE,
translate_to_asc = FALSE,
germline_db = NA,
novel = NA,
find_unmutated = FALSE,
seq = "sequence_alignment",
default_allele_threshold = 1e-04,
depth_adjusted_threshold = FALSE,
z_score_threshold = 0,
quiet = TRUE
)
```

Arguments
-------------------

data
:   data.frame in AIRR format, containing allele calls from a single subject and the sample IMGT-gapped V(D)J sequences under seq.

allele_threshold_table
:   A data.frame of the alleles and their thresholds.

call
:   name of the V,D, or J allele call column, i.e v_call, d_call, j_call. Default is `v_call`

asc_annotation
:   Logical (FALSE by default). Are the allele calls annotated with the allele similarity clusters.

single_assignment
:   if TRUE, the method only considers sequence with single assignment for the genotype inference.

translate_to_asc
:   For V allele calls, collapse identical allele for the genotype inference. Default is FALSE.

germline_db
:   named vector of sequences containing the germline sequences named in V allele calls and the alleleClusterTable. Only required if find_unmutated is TRUE.

novel
:   an optional `data.frame` of novel alleles, as returned by `tigger::findNovelAlleles` (columns `germline_call`, `polymorphism_call`, `novel_imgt`). When supplied, the novel germline sequences are added to `germline_db` and the novel alleles become genotype candidates. A novel allele inherits the threshold of its base allele (or the default when the base is absent). Default `NA` (no novel alleles).

find_unmutated
:   if TRUE, use germline_db to find which samples are unmutated. Not needed if V allele calls only represent unmutated samples. Only meaningful for V calls: D and J segments are heavily trimmed, so unmutated-call detection is biased and is skipped (with a warning) when the call column is a D or J segment.

seq
:   name of the column in data with the aligned, IMGT-numbered, V(D)J nucleotide sequence. Default is sequence_alignment.

default_allele_threshold
:   The default allele threshold for the genotype inference, in case the allele threshold is not in the `allele_threshold_table`. Default is 1e-04.

depth_adjusted_threshold
:   Logical (FALSE by default). If TRUE, each allele's presence threshold is raised to a depth-aware floor `max(Tai, default_allele_threshold, 1/N)`, where `Tai` is the allele threshold and `N` is the per-locus repertoire depth. This prevents shallow repertoires from clearing an unrealistically low threshold.

z_score_threshold
:   Numeric (0 by default). Alleles with `z_score >= z_score_threshold` are flagged as present in the returned `in_genotype` column. This is a flag only; no rows are dropped.

quiet
:   Logical (TRUE by default). Do you want to suppress informative messages




Value
-------------------

A a data.frame with the inferred V genotype. The table contains the following columns:

+  allele: The alleles in the `allele_threshold_table`.
+  counts: The number of reads for each alleles.
+  depth: The total number of reads in the genotype (Sum of counts).
+  threshold: The population driven allele thresholds for genotype presence.
+  z_score: The confidence level for the presence of the allele in the genotype.
+  observed: Logical, whether the allele was seen in the data (in single or multiple assignment). Alleles that were never seen carry a smoothing pseudo-count and are marked FALSE.
+  in_genotype: Logical, whether the allele passes `z_score_threshold` (`z_score >= z_score_threshold`).
+  asc_allele: If `translate_to_asc` is true, the asc allele value from allele_threshold_table.



Details
-------------------

In naive repertoires, allele calls where more than one assignment is assigned is rare. Hence, in case the data represents the naive repertoire of a subject
it is recommended to use the `find_unmutated=TRUE` option, to remove mutated sequences. For non-naive population, the allele calls in cases of multiple assignment
are treated as belonging to all groups.



Examples
-------------------

```R
# loading TIgGER AIRR-seq b cell data
data <- tigger::AIRRDb

# allele threshold table
data(allele_threshold_table)

data(HVGERM)

# inferring the genotype
genotype <- inferGenotypeAllele(
data = data,
allele_threshold_table = allele_threshold_table,
germline_db = HVGERM, find_unmutated=TRUE)

# filter alleles with z_score >= 0 

head(genotype[genotype$z_score >= 0,])

```


```
Key: <allele>
       gene      allele count    depth threshold   z_score
     <char>      <char> <num>    <num>     <num>     <num>
1: IGHV1-18 IGHV1-18*01  1005 4738.942     1e-03  459.7163
2:  IGHV1-2  IGHV1-2*02   664 4738.942     1e-03  302.9940
3:  IGHV1-2  IGHV1-2*04   302 4738.942     1e-04  438.0321
4: IGHV1-24 IGHV1-24*01   105 4738.942     1e-04  151.8468
5:  IGHV1-3  IGHV1-3*01   226 4738.942     1e-05 1037.9557
6: IGHV1-46 IGHV1-46*01   624 4738.942     1e-03  284.6101

```



See also
-------------------

[inferAlleleClusters](inferAlleleClusters.md) will infer the allele clusters based on a supplied V reference set and set the default allele threshold of 1e-04.
See [recentAlleleClusters](recentAlleleClusters.md) to obtain the latest version of the IGHV allele clusters and the naive population based allele threshold.






