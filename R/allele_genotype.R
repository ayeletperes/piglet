# ------------------------------------------------------------------------------
# Allele clusters functions
# The functions in this scripts are for generating the allele clusters from a given reference set

#' @include alleleclustergenotype.R
NULL


## TODO: decide if a class is needed

# ------------------------------------------------------------------------------

#' Download the most recent allele clusters reference set and thresholds based on \href{https://yaarilab.github.io/IGHV_reference_book}{https://yaarilab.github.io/IGHV_reference_book} 
#' At the moment onlu available for human IGHV reference set.
#'
#' @param version    The version of the allele cluster to download. Default is `latest`.
#' 
#' @return
#' An object of type \code{GermlineCluster} that includes the allele cluster table \code{alleleClusterTable} with the new names and the default thresholds, 
#' the renamed germline set \code{alleleClusterSet}, the germline set hierarchical clustering \code{hclustAlleleCluster},
#' and the threshold parameters \code{threshold}.
#'
#'

# TODO: the reference will need to be pulled from zenodo
recentAlleleClusters <- function(version = "latest"){
  
}

# TODO: a function to get the version number from zenodo.

# ------------------------------------------------------------------------------

#' \code{assignAlleleClusters} uses the allele clusters annotation to change the preliminary allele 
#' assignments to the new annotations before inferring a genotype.
#'
#' @param data                  data.frame in AIRR format, containing V allele calls from a single subject and the sample IMGT-gapped V(D)J sequences under seq.
#' @param v_call                name of the V allele call column. Default is `v_call`
#' @param alleleClusterTable    A data.frame of the allele clusters new annotations relative to the original reference set. See details.
#' 
#' @return
#' A modified input \code{data.frame} with the new assigned 
#'
#' @export
assignAlleleClusters <- function(data, alleleClusterTable, v_call = "v_call"){
  
  ## set the dictionary
  germline_set <- setNames(alleleClusterTable$new_allele, alleleClusterTable$original_allele)
  
  # switch the assignments
  data[[v_call]] <- sapply(data[[v_call]], function(x){
    calls <- unlist(strsplit(x, ","))
    calls <- germline_set[calls]
    calls <- calls[!duplicated(calls)]
    paste0(calls, collapse = ",")
  }, USE.NAMES = F)
  
  return(data)
}

# ------------------------------------------------------------------------------

#' \code{inferGenotypeAllele} infer an individual's genotype based on the allele-base method.
#' The method utilize the allele specific threshold to determine the presence of an allele in the genotype.
#' More specifically, the absolute frequency of each allele is calculated and checked against the threshold. 
#'
#' @param data                  data.frame in AIRR format, containing V allele calls from a single subject and the sample IMGT-gapped V(D)J sequences under seq.
#' @param alleleClusterTable    A data.frame of the allele clusters new annotations relative to the original reference set. 
#' @param v_call                name of the V allele call column. Default is `v_call`
#' @param single_assignment     if TRUE, the method only considers sequence with single assignment for the genotype inference. 
#' @param germline_db           named vector of sequences containing the germline sequences named in V allele calls and the alleleClusterTable. Only required if find_unmutated is TRUE.
#' @param find_unmutated        if TRUE, use germline_db to find which samples are unmutated. Not needed if V allele calls only represent unmutated samples.
#' @param seq                   name of the column in data with the aligned, IMGT-numbered, V(D)J nucleotide sequence. Default is sequence_alignment.
#' 
#' @return
#' A a data.frame with the inferred V genotype. The table contains the following columns:
#' 	|gene            | alleles             | imgt_alleles          | counts              | absolute_fraction     | absolute_threshold               | genotyped_alleles    | genotype_imgt_alleles|               
#'  |----------------|---------------------|-----------------------|---------------------|-----------------------|----------------------------------|----------------------|----------------------|
#'  | allele cluster | the present alleles | the imgt nomenclature | the number of reads | the absolute fraction | the population driven allele     | the alleles which    | the imgt nomenclature| 
#'  |                | in the repertoire   | of the alleles        | for each alleles    | of the alleles        | thresholds for genotype presence | entered the genotype | of the alleles 		  |
#'
#' @details 
#' 
#' In naive repertoires, allele calls where more than one assignment is assigned is rare. Hence, in case the data represents the naive repertoire of a subject
#' it is recommended to use the `find_unmutated=TRUE` option, to remove mutated sequences. For non-naive population, the allele calls in cases of multiple assinment
#' are treated as belonging to all groups.
#' 
#' @seealso 
#' 
#' \link{inferAlleleClusters} will infer the allele clusters based on a supplied V reference set and set the default allele threshold of 1e-04.
#' See \link{recentAlleleClusters} to obtain the latest version of the IGHV allele clusters and the naive population based allele threshold.
#' 
#' @export
# Parts are adapted from tigger::inferGenotype
inferGenotypeAllele <- function(data, alleleClusterTable, v_call = "v_call", single_assignment = FALSE, germline_db = NA, find_unmutated = FALSE, seq = "sequence_alignment"){
  
  
  . = NULL
  
  allele_calls = alakazam::getAllele(data[[v_call]], first = FALSE, 
                           strip_d = FALSE)
  
  ## check that the allele calls are in the supplied alleleClusterTable
  unique_calls <- unique(unlist(strsplit(allele_calls,",")))
  match <- unique_calls %in% alleleClusterTable$new_allele
  if(!all(match)){
    stop("The are allele calls that are not in the alleleClusterTable. Please check the allele call column.")
  }
  
  ## check unmutated
  if (find_unmutated) {
    if (is.na(germline_db[1])) {
      stop("germline_db needed if find_unmutated is TRUE")
    }
    allele_calls <- tigger::findUnmutatedCalls(allele_calls, as.character(data[[seq]]), 
                                       germline_db)
    if (length(allele_calls) == 0) {
      stop("No unmutated sequences found! Set 'find_unmutated' to 'FALSE'.")
    }
  }
  
  geno_V <- data.table::data.table(v_call = allele_calls)
  
  geno_V <-
    geno_V[, "gene" := alakazam::getGene(get("v_call"),
                                       first = F,
                                       collapse = T,
                                       strip_d = F)]
  
  # removing multiple gene assignments
  geno_V <- geno_V[!grepl(",", get("gene"))]
  
  # clean the allele class
  geno_V[, "v_allele" :=
           gsub("(IG[HKL][VDJADEGMC]|TR[ABDG])[A-Z0-9\\(\\)]+[-/\\w]*[*]",
                "",
                get("v_call"),
                perl = T)]
  
  ## get single assignments
  if(single_assignment){
    geno_V <- geno_V[!(grepl(",", get("v_allele")))]
    geno_V <- geno_V[!is.na(get("v_allele"))]
    geno_V[, "n_row_sub" := .N]
    geno_V[, "frac" := 1]
    geno_V_fraction <-
      geno_V[, .("absolute_fraction" = round(sum(get("frac")) / unique(get("n_row_sub")), 8),
                 "count" = sum(get("frac"))), by = list(get("gene"), get("v_allele"))]
  }else{
  ### distribute the multiple assignments for non naive sequences
    geno_V <- geno_V[!is.na(get("v_allele"))]
    geno_V[, "n_row_sub" := .N]
    geno_V[, "frac" := 1]
    geno_V <-
      geno_V[, .("count" = sum(get("frac"))), by = mget(c("gene", "v_allele", "n_row_sub"))]
    
    n_row_sub <- unique(geno_V$n_row_sub)
    
    geno_V_fraction <- c()
    #### code from TIgGER inferGentoype
    for(g in unique(geno_V$gene)){
      ac <- geno_V[get("gene")==g,get("v_allele")]
      t_ac <- setNames(geno_V[get("gene")==g,get("count")], geno_V[get("gene")==g,get("v_allele")]) # table of allele calls
      potentials <- unique(unlist(strsplit(names(t_ac),","))) # potential alleles
      
      regexpotentials <- paste(gsub("\\*","\\\\*", potentials),"$",sep="")
      regexpotentials <- 
        paste(regexpotentials,gsub("\\$",",",regexpotentials),sep="|")
      tmat <- 
        sapply(regexpotentials, function(x) grepl(x, names(t_ac),fixed=FALSE))
      
      if (length(potentials) == 1 | length(t_ac) == 1){ 
        seqs_expl <- t(as.data.frame(apply(t(as.matrix(tmat)), 2, function(x) x * 
                                             t_ac)))
        rownames(seqs_expl) <- names(t_ac)[1]
      }else{
        seqs_expl <- as.data.frame(apply(tmat, 2, function(x) x * 
                                           t_ac))
      }
      #       seqs_expl = as.data.frame(apply(tmat, 2, function(x) x*t_ac))
      colnames(seqs_expl) <- potentials
      # Add low (fake) counts
      sapply(colnames(seqs_expl), function(x){if(sum(rownames(seqs_expl) %in% paste(x)) == 0){
        seqs_expl <<- rbind(seqs_expl,rep(0,ncol(seqs_expl))); 
        rownames(seqs_expl)[nrow(seqs_expl)] <<- paste(x)
        seqs_expl[rownames(seqs_expl) %in% paste(x),paste(x)] <<- 0.01
        
      }}) 
      
      # Build ratio dependent allele count distribution of multi assigned reads
      seqs_expl_single <- seqs_expl[grep(',',rownames(seqs_expl),invert = T),] 
      
      seqs_expl_multi <- seqs_expl[grep(',',rownames(seqs_expl),invert = F),] 
      if(is.null(nrow(seqs_expl_multi))){
        seqs_expl_multi <- t(as.data.frame(seqs_expl_multi))
        rownames(seqs_expl_multi) <- grep(',',rownames(seqs_expl),invert = F,value = T)
      }
      
      if(!is.null(nrow(seqs_expl_single))  && nrow(seqs_expl_single) !=0 && nrow(seqs_expl_single) != nrow(seqs_expl)){
        if(nrow(seqs_expl_multi)>1){
          seqs_expl_multi <- seqs_expl_multi[order(nchar(row.names(seqs_expl_multi))),]
        }
        sapply(1:nrow(seqs_expl_multi),function(x){
          genes <- unlist(strsplit(row.names(seqs_expl_multi)[x],','));
          counts <- seqs_expl_single[rownames(seqs_expl_single) %in% genes,genes]
          counts <- colSums(counts)
          counts_to_distribute <- seqs_expl_multi[x,genes]
          
          new_counts <- counts+((counts_to_distribute*counts)/sum(counts))
          for(i in 1:length(new_counts)){
            gene_tmp <- names(new_counts)[i] 
            seqs_expl_single[rownames(seqs_expl_single) %in% gene_tmp,gene_tmp] <<- new_counts[i]
          }
        })
      } 
      
      seqs_expl <- if(is.null(nrow(seqs_expl_single)) || nrow(seqs_expl_single) ==0 ){seqs_expl}else{seqs_expl_single}
      seqs_expl <- round(seqs_expl)
      if(sum(rowSums(seqs_expl) == 0 ) != 0){
        seqs_expl <- seqs_expl[rowSums(seqs_expl)!= 0, ]
      }
      
      allele_tot <- sort(apply(seqs_expl, 2, sum),decreasing=TRUE)
      
      gene_table <- data.table::data.table("gene" = g, "v_allele" = names(allele_tot), "count" = allele_tot, "n_row_sub" = n_row_sub)
      gene_table <- gene_table[get("count")!=0]
      
      geno_V_fraction <- bind_rows(geno_V_fraction, gene_table)
    }
    ############
    geno_V_fraction <-
      geno_V_fraction[, .("absolute_fraction" = round(get("count") / unique(get("n_row_sub")), 8),
                          "count" = get("count")), by = mget(c("gene", "v_allele"))]
    
  }
  
  ## add original allele and cut off
  geno_V_fraction[, "v_call" := paste0(get("gene"), "*", get("v_allele"))]
  
  alleles_clusters <-
    setNames(alleleClusterTable$original_allele, alleleClusterTable$new_allele)
  geno_V_fraction[, "v_call_or" := alleles_clusters[get("v_call")]]
  
  
  na_id <- which(is.na(geno_V_fraction$v_call_or))
  if (length(na_id) != 0)
    geno_V_fraction$v_call_or[na_id] <- sapply(na_id, function(i) {
      new_allele <- geno_V_fraction$v_call[i]
      closest <- strsplit(geno_V_fraction$v_call[i], "_")[[1]][1]
      or_allele <- alleles_clusters[closest]
      paste0(or_allele, gsub(closest, "", new_allele, fixed = T))
    })
  
  allele_cluster_threshold <- setNames(alleleClusterTable$thresh, alleleClusterTable$new_allele)
  geno_V_fraction <-
    geno_V_fraction[, "absolute_thresh" := allele_cluster_threshold[get("v_call")]]
  na_id <- which(is.na(geno_V_fraction$absolute_thresh))
  if (length(na_id) != 0)
    geno_V_fraction$absolute_thresh[na_id] <-
    sapply(na_id, function(i) {
      new_allele <- geno_V_fraction$v_call[i]
      closest <- strsplit(geno_V_fraction$v_call[i], "_")[[1]][1]
      allele_cluster_threshold[closest]
    })
  
  ## check if allele is above thresh.
  
  geno_V_fraction[, "above_thresh" := get("absolute_fraction") >= get("absolute_thresh")]
  
  
  sortBy <- c('gene', 'absolute_fraction')
  sortType <- c( 1, -1)
  data.table::setorderv(geno_V_fraction, sortBy, sortType)
  genoV <-
    geno_V_fraction[, .(
      "alleles" = paste0(get("v_allele"), collapse = ","),
      "imgt_alleles" = paste0(get("v_call_or"), collapse = ","),
      "counts" = paste0(get("count"), collapse = ","),
      "absolute_fraction" = paste0(round(get("absolute_fraction"), 7), collapse = ","),
      "absolute_threshold" = paste0(formatC(get("absolute_thresh"), format = "f"), collapse = ","),
      "genotyped_alleles" = paste0(get("v_allele")[get("absolute_fraction") >=
                                            get("absolute_thresh")], collapse = ","),
      "genotyped_imgt_alleles" = paste0(get("v_call_or")[get("absolute_fraction") >=
                                                  get("absolute_thresh")], collapse = ",")
    ), by = mget(c("gene"))]
  
  
  # geno <- as.data.frame(genotype, stringsAsFactors = FALSE)
  # if (find_unmutated == TRUE) {
  #   seqs <- genotypeFasta(geno, germline_db)
  #   dist_mat <- seqs %>% sapply(function(x) sapply((getMutatedPositions(seqs, 
  #                                                                       x)), length)) %>% as.matrix
  #   rownames(dist_mat) <- colnames(dist_mat)
  #   for (i in 1:nrow(dist_mat)) {
  #     dist_mat[i, i] = NA
  #   }
  #   same <- which(dist_mat == 0, arr.ind = TRUE)
  #   if (nrow(same) > 0) {
  #     for (r in 1:nrow(same)) {
  #       inds <- as.vector(same[r, ])
  #       geno[getGene(rownames(dist_mat)[inds][1]), ]$note <- paste(rownames(dist_mat)[inds], 
  #                                                                  collapse = " and ") %>% paste("Cannot distinguish", 
  #                                                                                                .)
  #     }
  #   }
  # }
  # rownames(geno) <- NULL
  
  return(geno_V)
}




