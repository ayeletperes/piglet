# README #

The code in the repository is for creating an R package for the allele clusters (ACs) and the allele-based genotype.

The Scripts folder will hold the function of the package. The functions are divided into several sections.

## Functions

1. Allele Clusters creation
	* Naming description: The names for the new allele cluster will be as follows: IGHVF1-G1\*01 - IGH = chain, V = region, F1 = family cluster [:number:], G1 - allele cluster [:number:], and 01 = allele numbering (given by clustering order, no connection to the expression)
	* A function to calculate the distances given an IMGT aligned allele reference set
		* Input: An IMGT aligned germline set
		* Output: Distance matrix - based on levenshtein
	* A function to cluster the alleles based on two thresholds - 75% for family and 95% for allele-cluster
		* Input: A distance matrix from the function above
		* Output: A data.frame with columns for the alleles, the 75% threshold, and the 95% threshold
	* A function that re-names alleles base on the new cluster and gives the default absolute fraction threshold
		* Input: An IMGT aligned germline set, and the data.frame from the function above
		* Output: A data.frame with columns for the original allele name, the new allele annotation, the family cluster, the allele cluster, the default threshold.
	* A wrapper function `inferAlleleClusters` which runs all the function above by order to produce the new germline set and thresholds
		* Input: An IMGT aligned germline set
		* Output: An IMGT aligned germline set with the new names, and a data.frame of the allele cluster annotations with the threhsolds.
		* Note: Maybe use S4 object, this will hold the data.frame, germline set, and the hclust object for ploting?	
	* A function to plot the clustering with the clustering thresholds?
		* Input: The S4 object from the function above
		* Output: An Hierarchical Clustering plot with dotted lines for the thresholds and colored boxes around the clustered alleles.
2. Allele-base genotype
	* This part relies on the user running by himself the repertoire through igblast with the new reference set and obtaining an AIRR format file
	* A function to scrape the most updated allele cluster thresholds for the chains from the docker bitbucker repository
		* Input: Ig chain. Default is `IGH`. currently supports only IGH
		* Output: a data.frame of the allele cluster thresholds and annotations
	* A function that re-names the V column annotation based on a new reference set. This function is time consuming, it will havely relay on regex. Not recommended.
		* Input: An AIRR complient data.frame, the V column name, and the allele cluster table. Default for the V column is: `v_call`; If the allele table is not supplied the previous step is activated and re-names the alleles based on default germline set. 
		* Output: a data.frame with the re-named annoation in a new column - `v_call_new`
		* Note: we can use the code from tigger with updates
	* A function `inferGenotypeAllele` to genotype the repertoire based on the allele-based method
		* Input: An AIRR complient data.frame, the V column name, the allele cluster table. Default for the V column is: `v_call`; If the allele table is not supplied the previous step is activated and re-names the alleles based on default germline set. 
		* Output: a data.frame with the V genotype based on the thresholds.
		* Note: Should we add unmutated option as in tigger? we do address in the pipeline if the sequences are naive or not. If they are then we take unmutated sequences and those that have one assignment.
		* Value: V genotype table:
		
			   |gene            | alleles             | imgt_alleles          | counts              | absolute_fraction     | absolute_threshold               | genotyped_alleles    | genotype_imgt_alleles|               
			   |----------------|---------------------|-----------------------|---------------------|-----------------------|----------------------------------|----------------------|----------------------|
			   | allele cluster | the present alleles | the imgt nomenclature | the number of reads | the absolute fraction | the population driven allele     | the alleles which    | the imgt nomenclature| 
			   |                | in the repertoire   | of the alleles        | for each alleles    | of the alleles        | thresholds for genotype presence | entered the genotype | of the alleles 		 |
				 
	* A genotype plot? we need to think what is suitable here
3. Allele Cluster reporting
	* A function to plot the frequency of the top N combination of alleles from each allele cluster compared to the overall fraction. This function is not sutiable for one sample
		* Note: Here we can maybe utilize P11, and to plot the new sample against the results from there. 
		
## General notes

* We can wrap some of the functions from tigger to produce the same outputs or we can send the user there to use. Such as personal reference set in fasta format.


	
	
