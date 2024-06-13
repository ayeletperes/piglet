#' .onAttach start message
#'
#' @param libname defunct
#' @param pkgname defunct
#'
#' @return invisible()
.onAttach <- function(libname, pkgname) {
  msg <- paste0("PIgLET version: ",packageVersion(pkgname))
  msg <- paste(msg, 'New feature was added! A confidence level to the genotype inference. Check the news for more details', collapse = "\n\n")
  cite <- citation(pkgname)
  msg <-paste(msg,c(format(cite,"citation")),collapse="\n\n")
  packageStartupMessage(msg)
  invisible()
}
