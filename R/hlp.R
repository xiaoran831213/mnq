#' Force Postive Definite
fpd <- function(x, eps=NULL)
{
    if(is.null(eps))
        eps <- sqrt(.Machine$double.eps)
    . <- svd(x)
    u <- .$u
    v <- .$v
    d <- .$d

    . <- d > d[1] * eps
    u <- u[, .]
    v <- v[, .]
    d <- d[.]

    x <- u %*% (d * t(v))       # U diag(d) V'
    x <- 0.5 * (x + t(x))

    x
}

#' Kernel Scaling
ksc <- function(a)
{
    a - outer(rowMeans(a), colMeans(a), `+`) + mean(a)
}
