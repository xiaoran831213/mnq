#' Generate Multi-Variant-Normal Samples
#'
#' A simplified copy of MASS::mvrnorm.
#'
#' @param N number of samples to be drawn
#' @param M vector of mean, to be rotated to the length of sample
#' features \code{D}
#' @param V square matrix of covariance, with \code{D} dimensions
#' @param drop TRUE to drop matrix of single sample to vector
#' 
#' @return matrix of N row samples and D column features, whose
#' covariance is \code{V}
mvn <- function (N=1, M=0, V=rbind(c(1, .5), c(.5, 1)), drop=TRUE)
{
    ## dimensionality
    D <- nrow(V)

    ## mean vector
    M <- drop(rep(M, length=D))

    ## eigen decomposition
    .e <- eigen(V, symmetric = TRUE)
    ev <- .e$values
    ec <- .e$vectors
    ev <- sqrt(pmax(ev, 0))          # PSD projection

    ## random values
    X <- matrix(rnorm(D * N), D, N)
    y <- t(M + ec %*% (ev * X))

    if(drop)
        y <- drop(y)
    y
}


#' Simulation 1: Cpp MINQUE versus R-language MINQUE
#'
#' @param N sample size
#' @param P number of features
#' @export
sm1 <- function(N=2000, P=2 * N)
{
    Z1 <- as.matrix(scale(matrix(rnorm(N * P), N, P)))
    Z2 <- as.matrix(scale(matrix(rnorm(N * P), N, P)))
    knl <- within(list(),
    {
        EPS <- diag(N)
        KN1 <- tcrossprod(Z1) / P
        KN2 <- tcrossprod(Z2) / P
        KN1 <- KN1 / mean(diag(KN1))
        KN2 <- KN2 / mean(diag(KN2))
    })
    K <- length(knl)
    
    ## allowed kernels
    v <- knl[c('EPS', 'KN1', 'KN2')]
    w <- matrix(0, 0, 0)
    
    ## true covariance
    X <- matrix(1, N, 1)
    V <- with(v, .8 * EPS + .5 * KN1 + 1.5 * KN2)
    y <- as.matrix(mvn(1, -1.3, V))
    
    ## minque
    t1 <- system.time(r1 <- rln_mnq(y, v, X, NULL))
    t2 <- system.time(r2 <- .Call('amo_mnq', y, v, X, w, PACKAGE='mnq'))

    print(list(tm.rln=t1, tm.amo=t2))

    list(r1=r1, r2=r2)
}

#' Simulation 2: MINQUE main function
#'
#' @param N sample size
#' @param P number of features
#' @param ... additional argument for \code{mnq} function
#' @export
sm2 <- function(N=2000, P=2 * N, ...)
{
    Z1 <- as.matrix(scale(matrix(rnorm(N * P), N, P)))
    Z2 <- as.matrix(scale(matrix(rnorm(N * P), N, P)))
    knl <- within(list(),
    {
        EPS <- diag(N)
        KN1 <- tcrossprod(Z1) / P
        KN2 <- tcrossprod(Z2) / P
    })
    K <- length(knl)

    ## allowed kernels
    v <- knl[c('KN1', 'KN2')]
    
    ## true covariance
    X <- cbind(X00=rep(1, N))
    S <- with(knl, .5 * EPS + .8 * KN1 + 1.3 * KN2)
    y <- mvn(1, -1.2, S)
    
    ## minque
    mdl <- mnq(y, v, X, ...)
    mdl
}

#' Simulation 3: Foward selection for MINQUE
#'
#' @param N sample size
#' @param P number of features
#' @param ... additional argument for \code{mnq} function
#' @export
sm3 <- function(N=2000, P=2 * N, ...)
{
    dot <- list(...)
    set.seed(dot$seed)

    ## true vcs
    w <- c(EPS=3.00, LN1=0.67,  LN2=0.36,  LN3=3.57,  LN4=0.66,  LN5=0.46, LN6=1.11)

    ## data
    zs <- replicate(length(w) - 1,
                    as.matrix(scale(matrix(rnorm(N * P), N, P))), simplify=FALSE)
    ks <- lapply(zs, function(z) tcrossprod(z) / P)
    names(ks) <- sprintf("LN%d", seq_along(ks))
    ks <- c(list(EPS=diag(N)), ks)
    
    ## allowed kernels
    v <- ks[-1]
    
    ## true covariance
    x <- cbind(X00=rep(1, N))
    names(w) <- names(ks)
    S <- .vsum(ks, w)
    S <- with(svd(S), u %*% diag(pmax(0, d)) %*% t(v))

    ## response
    y <- mvn(1, -1, S)

    ## minque
    ref <- within(list(),
    {
        par <- c(X00=-1, w)
        rpt <- vpd(y, v, x, par, ...)
    })

    kcv <- fcv(y, v, x, ...)
    
    ful <- within(list(),
    {
        par <- mnq(y, v, x, ...)$par
        par[names(par) %in% c("EPS", names(v)) & par < 0] <- 0
        rpt <- vpd(y, v, x, par, ...)
    })

    sel <- fwd(y, v, x, ...)

    list(sel=sel, ref=ref, ful=ful, kcv=kcv)
}
