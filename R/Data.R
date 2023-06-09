# Documentation and definitions for data and constants

#### Data ####

#' Human IGHV germlines
#'
#' A \code{character} vector of all 498 human IGHV germline gene segment alleles
#' in IMGT Gene-db release July 2022, with an additional 25 undocumented alleles from VDJbase.
#'
#' @name HVGERM
#' @docType data
#' @format Values correspond to IMGT-gaped nuceltoide sequences (with
#' nucleotides capitalized and gaps represented by '.').
#'
#' @references Xochelli \emph{et al}. (2014) Immunoglobulin heavy variable
#' (IGHV) genes and alleles: new entities, new names and implications for
#' research and prognostication in chronic lymphocytic leukaemia.
#' \emph{Immunogenetics}. 67(1):61-6.
#' @keywords data
"HVGERM"

#' Human IGHV germlines functionality description
#'
#' A \code{data.table} of all 498 human IGHV germline gene segment alleles
#' in IMGT Gene-db release July 2022, with an additional 25 undocumented alleles from VDJbase.
#' The first column is the allele name, the second column is the functionality annotation, the 
#' third column is the nt sequence and the last column is the aa sequence.
#'
#' @name hv_functionality
#' @docType data
#'
#' @references Xochelli \emph{et al}. (2014) Immunoglobulin heavy variable
#' (IGHV) genes and alleles: new entities, new names and implications for
#' research and prognostication in chronic lymphocytic leukaemia.
#' \emph{Immunogenetics}. 67(1):61-6.
#' @keywords data
"hv_functionality"
