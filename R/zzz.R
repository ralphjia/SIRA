# =============================================================================
# R/zzz.R
# Package-level documentation and Rcpp shared library loader.
#
# zzz.R is loaded last by R's package loader (alphabetically after all other
# R/ files), making it the conventional place for:
#   - The package-level roxygen2 docstring (creates the ?SIRA help page)
#   - .onLoad() / .onAttach() hooks
#   - Rcpp module registration via loadModule() if needed
# =============================================================================


#' @keywords internal
"_PACKAGE"

#' @importFrom Matrix Matrix sparseMatrix Diagonal rowSums
#' @useDynLib SIRA, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
