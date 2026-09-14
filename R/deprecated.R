# Deprecated names for the C++ sequence-comparison helpers.
#
# These are kept exported and working indefinitely: they are called from roughly forty
# analysis scripts outside this package, most of which are not under version control.
# Each one delegates to the current function, so there is one implementation to maintain.

#' Deprecated: use \code{\link{allele_diff_paired}}
#'
#' Renamed. The old name said \code{parallel}, which was an argument rather than a
#' behaviour and never had any effect (the package is not built with OpenMP), and
#' \code{indices}, which is wrong whenever \code{return_count = TRUE}.
#'
#' Results are unchanged: this is the same implementation under a clearer name.
#'
#' @inheritParams allele_diff_paired
#' @param parallel Ignored. Accepted so existing calls keep working.
#' @return As \code{\link{allele_diff_paired}}.
#' @seealso \code{\link{allele_diff_paired}}
#' @export
allele_diff_indices_parallel2 <- function(germs, inputs, X = 0L, parallel = FALSE,
                                          return_count = FALSE,
                                          non_mismatch_chars_nullable = NULL) {
  .Deprecated("allele_diff_paired")
  allele_diff_paired(germs, inputs, X = X, return_count = return_count,
                     non_mismatch_chars_nullable = non_mismatch_chars_nullable)
}

#' Deprecated: use \code{\link{allele_diff_paired}}
#'
#' Renamed, and \strong{the counts change}. The old implementation ignored gaps and
#' ambiguous bases only on the germline side, so an aligned position where the germline
#' carried a base and the input carried a gap was counted as a mismatch. It also looped
#' to the length of the germline while indexing the input, reading past the end of the
#' input whenever the input was shorter -- undefined behaviour, and on real IGHV data it
#' returned counts derived partly from adjacent memory.
#'
#' Calls now route to \code{\link{allele_diff_paired}}, which ignores such positions on
#' either side and pads the shorter sequence. Against the old behaviour, 194 of 400 real
#' germline pairs differ, by a mean of 8 mismatches. Anything that depended on the old
#' numbers needs re-checking rather than re-running.
#'
#' @inheritParams allele_diff_paired
#' @param parallel Ignored. Accepted so existing calls keep working.
#' @return As \code{\link{allele_diff_paired}}.
#' @seealso \code{\link{allele_diff_paired}}
#' @export
allele_diff_indices_parallel <- function(germs, inputs, X = 0L, parallel = FALSE,
                                         return_count = FALSE) {
  .Deprecated("allele_diff_paired",
              msg = paste("'allele_diff_indices_parallel' is deprecated and its counts",
                          "have changed. It ignored gaps only on the germline side and",
                          "read past the end of shorter inputs. Use",
                          "'allele_diff_paired', and re-check any result that depended",
                          "on the old counts."))
  allele_diff_paired(germs, inputs, X = X, return_count = return_count)
}

#' Deprecated: use \code{\link{insert_gaps}}
#'
#' Renamed. The \code{2} marked a second attempt rather than anything about behaviour,
#' and \code{_vec} restated that it is vectorised. Results are unchanged.
#'
#' @inheritParams insert_gaps
#' @param parallel Ignored. Accepted so existing calls keep working.
#' @return As \code{\link{insert_gaps}}.
#' @seealso \code{\link{insert_gaps}}
#' @export
insert_gaps2_vec <- function(gapped, ungapped, parallel = FALSE) {
  .Deprecated("insert_gaps")
  insert_gaps(gapped, ungapped)
}
