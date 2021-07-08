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
#' @noRD
mvn <- function (N=1, M=0, V=NULL, drop=TRUE)
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
mtx <- function(N=500, P=2*N, ...)
{
    
    x <- matrix(rnorm(N * P), N, P)
    rownames(x) <- sprintf("I%02X", seq(1, l=N))
    colnames(x) <- sprintf("X%02X", seq(1, l=P))
    x
}

#' generate PSD matrix
#'
#' @param N number of row samples
#' @param P number of column features
#' @param f funtion to draw random samples (i.e., rnorm, rpois)
#' @param ... additional arguments for \code{f}.
#' @return a square PSD matrix of N rows/columns.
psd <- function(N=500, P=2*N, tau=NULL, scl=FALSE, ...)
{
    a <- mtx(N, P, ...)
    a <- tcrossprod(a) / P
    a <- 0.5 * (a + t(a))
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
knl <- function(N=500, P=N*2, L=2, ...)
{
    K <- replicate(L, psd(N, P), simplify=FALSE)
    K <- lapply(K, ksc)
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
#' @param fix fixed effect coefficients, intercept included
#' @param itc intercept
#' @export
sim <- function(N=500, P=N, vcs=c(1, 1), fix=c(1, 1), ...)
{
    arg <- get.arg()
    times <- arg$times %||% 1; arg$times <- NULL
    L <- length(vcs)
    M <- length(fix)
    res <- list()

    set.seed(arg$seed)
    for(i in seq(times))
    {
        X <- cbind(X00=rep(1, N), mtx(N, M - 1))    # covariants
        K <- c(list(EPS=diag(N)), knl(N, P, L - 1)) # kernels
        y <- mvn(1, X %*% fix, vsm(K, vcs))         # response
        r <- lm(y~X-1)
        p1 <- st0(r, K[[2]])
        p2 <- st1(r, K[[2]])
        res[[i]] <- with(r, .d(p1=p1, p2=p2))
    }
    set.seed(NULL)
    res <- .d(arg, do.call(rbind, res))
    res
}

st0 <- function(r, R)
{
    N <- nrow(R)
    rsd <- resid(r)
    mu2 <- mean(rsd^2)
    mu4 <- mean(rsd^4)
    
    Q <- t(rsd) %*% R %*% rsd / mu2

    EQ <- sum(diag(R))
    VQ <- 2 * sum(diag(R)^2) + (mu4/mu2^2 - 3) * sum(diag(R)^2)
    
    zsc <- drop((Q - EQ) / sqrt(VQ))
    pvl <- 2 * (1 - pnorm(abs(zsc)))

    pvl
}

st1 <- function(r, R)
{
    N <- nrow(R)
    rsd <- resid(r)
    sg2 <- var(rsd)

    qt <- t(rsd) %*% R %*% rsd / sg2

    rs <- ksc(R)
    eq <- tr(rs)
    vq <- 2 / (N + 1) * ((N - 1) * sum(diag(rs)^2) - sum(diag(rs))^2)
    zs <- drop((qt - eq) / sqrt(vq))
    
    pvl <- 2 * (1 - pnorm(abs(zs)))
    pvl
}

bmk <- function(N=500, P=N, vcs=c(1, 1, 1), fix=NULL, itc=NULL, mtd=rln_mnq, times=100)
{
    L <- length(vcs)
    M <- length(fix)
    C <- length(itc)

    ## fixed effect
    X <- cbind(X0000=rep(1, N * C), mtx(N, M)) # covariants
    m <- 0
    if(length(X) > 0)
        m <- X %*% c(itc, fix)

    ## kernel (random effect)
    K <- c(list(EPS=diag(N)), knl(N, P, L-1))  # kernels

    ## generate outcome
    y <- as.matrix(mvn(M=m, V=.vsum(K, vcs)))
    
    ## minque
    rtm <- replicate(times,
    {
        tm <- mtd(y, K, X)[['rtm']]
        tm["user.self"]
    })
    mean(rtm)
}
