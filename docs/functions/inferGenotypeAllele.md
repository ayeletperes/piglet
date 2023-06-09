# `inferGenotypeAllele`

`inferGenotypeAllele` infer an individual's genotype based on the allele-base method.
 The method utilize the allele specific threshold to determine the presence of an allele in the genotype.
 More specifically, the absolute frequency of each allele is calculated and checked against the threshold.


## Description

`inferGenotypeAllele` infer an individual's genotype based on the allele-base method.
 The method utilize the allele specific threshold to determine the presence of an allele in the genotype.
 More specifically, the absolute frequency of each allele is calculated and checked against the threshold.


## Usage

```r
inferGenotypeAllele(
  data,
  alleleClusterTable,
  v_call = "v_call",
  single_assignment = FALSE,
  germline_db = NA,
  find_unmutated = FALSE,
  seq = "sequence_alignment"
)
```


## Arguments

Argument      |Description
------------- |----------------
`data`     |     data.frame in AIRR format, containing V allele calls from a single subject and the sample IMGT-gapped V(D)J sequences under seq.
`alleleClusterTable`     |     A data.frame of the allele clusters new annotations relative to the original reference set.
`v_call`     |     name of the V allele call column. Default is `v_call`
`single_assignment`     |     if TRUE, the method only considers sequence with single assignment for the genotype inference.
`germline_db`     |     named vector of sequences containing the germline sequences named in V allele calls and the alleleClusterTable. Only required if find_unmutated is TRUE.
`find_unmutated`     |     if TRUE, use germline_db to find which samples are unmutated. Not needed if V allele calls only represent unmutated samples.
`seq`     |     name of the column in data with the aligned, IMGT-numbered, V(D)J nucleotide sequence. Default is sequence_alignment.


## Details

In naive repertoires, allele calls where more than one assignment is assigned is rare. Hence, in case the data represents the naive repertoire of a subject
 it is recommended to use the `find_unmutated=TRUE` option, to remove mutated sequences. For non-naive population, the allele calls in cases of multiple assignment
 are treated as belonging to all groups.


## Value

A a data.frame with the inferred V genotype. The table contains the following columns: list(list("llllllll"), list("\n", "   gene ", list(), " alleles ", list(), " imgt_alleles ", list(), " counts ", list(), " absolute_fraction ", list(), " absolute_threshold ", list(), " genotyped_alleles ", list(), " genotype_imgt_alleles ", list(), "\n", "   allele cluster ", list(), " the present alleles ", list(), " the imgt nomenclature ", list(), " the number of reads ", list(), " the absolute fraction ", list(), " the population driven allele ", list(), " the alleles which ", list(), " the imgt nomenclature ", 
    list(), "\n", "    ", list(), " in the repertoire ", list(), " of the alleles ", list(), " for each alleles ", list(), " of the alleles ", list(), " thresholds for genotype presence ", list(), " entered the genotype ", list(), " of the alleles ", list(), "\n"))


## Seealso

[inferAlleleClusters](#inferalleleclusters) will infer the allele clusters based on a supplied V reference set and set the default allele threshold of 1e-04.
 See [recentAlleleClusters](#recentalleleclusters) to obtain the latest version of the IGHV allele clusters and the naive population based allele threshold.


## Examples

```r
# loading TIgGER AIRR-seq b cell data
data <- tigger::AIRRDb

# getting the archive
asc_archive <- recentAlleleClusters(doi="10.5281/zenodo.7401239", get_file = TRUE)

# extracting the allele cluster table
allele_cluster_table <- extractASCTable(archive_file = asc_archive)

data(HVGERM)

# reforming the germline set
asc_germline <- germlineASC(allele_cluster_table, germline = HVGERM)

# assigning the ASC alleles
asc_data <- assignAlleleClusters(data, allele_cluster_table)

# inferring the genotype
asc_genotype <- inferGenotypeAllele(asc_data,
alleleClusterTable = allele_cluster_table,
germline_db = asc_germline, find_unmutated=T)
```


