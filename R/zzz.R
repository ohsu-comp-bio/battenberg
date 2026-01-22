#' Battenberg: Subclonal Copy Number Caller
#'
#' @useDynLib Battenberg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onLoad <- function(libname, pkgname) {
  # Keep your scipen setting
  options(scipen = 999)
}

.datatable.aware <- TRUE
