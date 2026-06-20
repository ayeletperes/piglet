# ------------------------------------------------------------------------------
# Utils functions

#' @include piglet.R
NULL

## Function to clean allele calls
clean_allele_calls <- function(segment_call, sep = ","){
  
  segment_regex <- "(IG[HKL]|TR[ABDG])[VDJADEGMC][A-R0-9()]*[-/\\w]*[-*]*[.\\w]+"
  nl_regex <- "(IG[HKL]|TR[ABDG])[VDJADEGMC][0-9]+-NL[0-9]([-/\\w]*[-*][.\\w]+)*"
  segments <- strsplit(segment_call, sep)
  r <- sapply(segments, function(segment){
    segment <- sub(paste0("[^", sep, "]*(", segment_regex, ")[^", sep, "]*"), "\\1", segment, perl = TRUE)
    segment[!grepl(nl_regex, segment)]
  }
  )
  return(r)
}

## Z score function for the allele based genotype
z_score <- function(Ni, N, Pi) {
  (Ni - Pi * N) / sqrt(Pi * N * (1 - Pi))
}

## Incorporate a novel-allele table (as returned by tigger::findNovelAlleles)
## into a set of comma-joined allele calls, mirroring tigger::inferGenotypeBayesian:
##  - append each novel `polymorphism_call` to the calls that carry its
##    `germline_call`, so the novel allele becomes a genotype candidate;
##  - add the novel IMGT-gapped germline sequences to `germline_db` so the
##    unmutated check can match them.
## Returns the augmented `allele_calls` (character vector), `germline_db`, and the
## cleaned `novel` data.frame (or NULL when no usable novels were supplied).
incorporate_novel <- function(allele_calls, germline_db, novel) {
  if (is.null(nrow(novel)) || all(is.na(novel))) {
    return(list(allele_calls = allele_calls, germline_db = germline_db, novel = NULL))
  }
  novel <- as.data.frame(novel, stringsAsFactors = FALSE)
  req <- c("germline_call", "polymorphism_call", "novel_imgt")
  if (!all(req %in% names(novel))) {
    stop("'novel' must contain the columns: ", paste(req, collapse = ", "))
  }
  novel <- novel[!is.na(novel[["polymorphism_call"]]), req, drop = FALSE]
  if (nrow(novel) == 0) {
    return(list(allele_calls = allele_calls, germline_db = germline_db, novel = NULL))
  }
  novel_gl <- novel[["novel_imgt"]]
  names(novel_gl) <- novel[["polymorphism_call"]]
  germline_db <- c(germline_db, novel_gl)
  for (r in seq_len(nrow(novel))) {
    ind <- grep(novel[["germline_call"]][r], allele_calls, fixed = TRUE)
    if (length(ind) > 0) {
      allele_calls[ind] <- paste(allele_calls[ind], novel[["polymorphism_call"]][r],
                                 sep = ",")
    }
  }
  list(allele_calls = allele_calls, germline_db = germline_db, novel = novel)
}

## Base (germline) allele of a possibly-novel allele name:
## "IGHV1-69*01_A45T" -> "IGHV1-69*01"; non-novel names are returned unchanged.
base_allele <- function(allele) sub("_.*$", "", allele)
