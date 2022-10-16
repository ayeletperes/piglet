#' .onAttach start message
#'
#' @param libname defunct
#' @param pkgname defunct
#'
#' @return invisible()
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("AlleleClusterGenotype version: ",packageVersion(pkgname)))
  invisible()
}
