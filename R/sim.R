#' Multivariant Normal Samples
#'
#' A simplified copy of MASS::mvrnorm.
#'
#' @param N number of samples to be drawn
#' @param  M vector  of mean,  to be rotated  to the  length of  sample features
#'     \code{D}
#' @param V square matrix of covariance, with \code{D} dimensions
#' @param drop TRUE to drop matrix of single sample to vector
#' 
#' @return matrix of N row samples and D column features, whose
#' covariance is \code{V}
.mvn <- function (N=1, M=0, V=NULL, drop=TRUE)
{
    ## default V is for demonstration
    if(is.null(V))
        V = rbind(c(1, .5, .25), c(.5, 1, .5), c(.25, .5, 1))

    ## dimensionality
    D <- nrow(V)

    ## mean vector
    M <- drop(rep(M, length=D))

    ## eigen decomposition
    e <- eigen(V, symmetric = TRUE)
    d <- e$values
    U <- e$vectors
    s <- sqrt(pmax(d, 0))          # square root of V

    ## random values
    X <- matrix(rnorm(D * N), D, N)
    y <- t(M + U %*% (s * X))

    if(drop)
        y <- drop(y)
    y
}

#' generate a matrix
#'
#' @param N number of row samples
#' @param P number of column features
#' @param f funtion to draw random numbers (i.e., rnorm, rpois)
#' @param ... additional arguments for \code{f}.
#' @return a matrix of N rows and P columns.
gen.mtx <- function(N=500, P=2*N, ...)
{
    
    x <- matrix(rnorm(N * P), N, P)
    rownames(x) <- sprintf("I%04X", seq(1, l=N))
    colnames(x) <- sprintf("X%04X", seq(1, l=P))
    x
}

#' generate PSD matrix
#'
#' @param N number of row samples
#' @param P number of column features
#' @param f funtion to draw random samples (i.e., rnorm, rpois)
#' @param ... additional arguments for \code{f}.
#' @return a square PSD matrix of N rows/columns.
gen.psd <- function(N=500, P=2*N, tau=NULL, scl=FALSE, ...)
{
    a <- gen.mtx(N, P, ...)
    a <- tcrossprod(a) / P
    a <- 0.5 * (a + t(a))
    if(scl)
        a <- a - outer(rowMeans(a), colMeans(a), `+`) + mean(a)
    if(is.null(tau))
        tau <- sqrt(.Machine$double.eps) #  * N
    a + diag(tau, N)
}

#' Generate kernels
#'
#' @param N number of row samples
#' @param P number of column features
#' @param L number of kernels
#' @param ... additional arguments for \code{f}.
#' @return a matrix of N rows and P columns.
gen.knl <- function(N=500, P=N*2, L=2, ...)
{
    K <- replicate(L, gen.psd(N, P), simplify=FALSE)
    names(K) <- sprintf("K%02d", seq_along(K))
    K
}

#' Simulation 1: Cpp MINQUE versus R-language MINQUE
#'
#' No fixed effect
#' 
#' @param N sample size
#' @param P number of features
#' @param vcs variance components, noise included
#' @param fix fixed coefficients
#' @param itc intercept
#' @export
sim <- function(N=500, P=N, vcs=c(1.0, 0.5, 1.5), fix=NULL, itc=NULL, mtd=rln_mnq, seed=NULL)
{
    L <- length(vcs)
    M <- length(fix)
    C <- length(itc)

    set.seed(seed)
    ## fixed effect
    X <- cbind(X0000=rep(1, N * C), gen.mtx(N, M)) # covariants
    m <- 0
    if(length(X) > 0)
        m <- X %*% c(itc, fix)

    ## kernel (random effect)
    K <- c(list(EPS=diag(N)), gen.knl(N, P, L-1))  # kernels

    ## generate outcome
    y <- as.matrix(.mvn(M=m, V=.vsum(K, vcs)))
    set.seed(NULL)
    
    ## minque
    Rprof(interval=1e-3)
    ret <- mtd(y, K, X)[c('vcs', 'fix', 'rtm')]
    Rprof(NULL)
    ret
}


bmk <- function(N=500, P=N, vcs=c(1, 1, 1), fix=NULL, itc=NULL, mtd=rln_mnq, times=100)
{
    L <- length(vcs)
    M <- length(fix)
    C <- length(itc)

    ## fixed effect
    X <- cbind(X0000=rep(1, N * C), gen.mtx(N, M)) # covariants
    m <- 0
    if(length(X) > 0)
        m <- X %*% c(itc, fix)

    ## kernel (random effect)
    K <- c(list(EPS=diag(N)), gen.knl(N, P, L-1))  # kernels

    ## generate outcome
    y <- as.matrix(.mvn(M=m, V=.vsum(K, vcs)))
    
    ## minque
    rtm <- replicate(times,
    {
        tm <- mtd(y, K, X)[['rtm']]
        tm["user.self"]
    })
    mean(rtm)
}
