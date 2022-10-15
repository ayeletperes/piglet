# alleleclustergenotype package documentation and import directives

#' The alleleclustergenotype package
#'
#' The \code{alleleclustergenotype} package provides two tools. The first is Allele Clusters,
#' this tool is designed to reduce the ambiguity within the IGHV alleles. The ambiguity
#' is caused by duplicated or similar alleles which are shared among different genes.
#' The second tool is an allele based genotype, that determined the presence of an allele based on
#' a threshold derived from a naive population.
#'
#' @section  Allele Cluster:
#' This section provides the functions that support the main tool of creating the allele cluster form
#' an IGHV germline set.
#'
#' \itemize{
#'   \item  \link{inferAlleleClusters}:      The main function of the section to create the allele clusters based on a germline set.
#'   \item  \link{ighvDistance}:             Calculate the distance between IGHV aligned germline sequences..
#'   \item  \link{ighvClust}:                Hierarchical clustering of the distance matrix from `ighvDistance`.
#'   \item  \link{generateReferenceSet}:     Generate the allele clusters reference set.
#'   \item  \link{plotAlleleCluster}:        Plots the Hierarchical clustering.
#' }
#'
#' @section  Allele based genotype:
#' This section provides the functions to infer the IGHV genotype using
#' the allele based method and the allele clusters thresholds
#'
#' \itemize{
#'   \item  \link{inferGenotypeAllele}:      Infer the IGHV genotype using the allele based method.
#'   \item  \link{updateThresh}:             Download the most recent version of the allele clusters threshold.
#'   \item  \link{assignAlleleClusters}:     Renames the v allele calls based on the new allele clusters.
#' }
#'
#' @name     alleleclustergenotype
#' @docType  package
#' @references
#' \enumerate{
#'   \item
#'  }
#'
#' @import   methods
#' @import   utils
#' @import   dendextend
#' @import   dplyr
#' @import   ggplot2
#' @import   circlize
#' @importFrom  data.table       := rbindlist data.table .N setDT CJ setorderv setkey .SD
#' @importFrom  stats            hclust as.dendrogram as.dist binom.test p.adjust setNames weighted.mean
#' @importFrom  alakazam         getGene
#' @importFrom  rlang            .data
#' @importFrom  tigger           readIgFasta
#' @importFrom  Biostrings       DNAStringSet
#' @importFrom  DECIPHER         DistanceMatrix
#' @importFrom  RColorBrewer     brewer.pal.info brewer.pal
NULL
