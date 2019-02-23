#' @importFrom stats lm residuals rnorm cor sd
#' @importFrom Rcpp evalCpp
#' @useDynLib mnq
NULL

.onUnload <- function(libpath)
{
    print(paste("Unloading", libpath))
    library.dynam.unload("mnq", libpath)
}
