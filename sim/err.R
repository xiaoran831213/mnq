#' mean correlation error
#'
#' @param u correlation estimates
#' @param v correlation truth
mce <- function(u, v, diag=TRUE)
{
    mean(abs(u - v)[upper.tri(v, diag=diag)])
}

#' mean absolute error
#'
#' @param h estimates
#' @param x truth
mae <- function(h, x)
{
    mean(abs(h - x), na.rm=TRUE)
}

#' mean square error
#'
#' @param h estimates
#' @param x truth
mse <- function(h, x)
{
    mean((h - x)^2, na.rm=TRUE)
}
