#' Multivariant Normal Samples
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
#' @param 
#' @export
sm1 <- function(N=500, P=N, vcs=c(1.0, 0.5, 1.5), fix=NULL, itc=NULL)
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
    op <- list(gls=1)
    t0 <- Sys.time()
    r1 <- rln_mnq(y, K, X)
    t1 <- Sys.time()
    r2 <- try(.Call('amo_mnq', y, K, X, matrix(0, 0, 0), op, PACKAGE='mnq'))
    t2 <- Sys.time()
    r3 <- rln_mq0(y, K, X)
    t3 <- Sys.time()
    r4 <- try(.Call('amo_mn0', y, K, X, op, PACKAGE='mnq'))
    t4 <- Sys.time()

    t4 <- t4 - t3; units(t4) <- 'secs'; t4 <- as.numeric(t4)
    t3 <- t3 - t2; units(t3) <- 'secs'; t3 <- as.numeric(t3)
    t2 <- t2 - t1; units(t2) <- 'secs'; t2 <- as.numeric(t2)
    t1 <- t1 - t0; units(t1) <- 'secs'; t1 <- as.numeric(t1)

    cat("Runing Time:\n")
    cat("R MN1:", t1, "\n")
    cat("C MN1:", t2, "\n")
    cat("R MN0:", t3, "\n")
    cat("C MN0:", t4, "\n")
    list(r_mn1=r1, c_mn1=r2, r_mn0=r3, c_mn0=r4)
}


#' Simulation 2: MINQUE main function
#'
#' @param N sample size
#' @param P number of features
#' @param ... additional argument for \code{mnq} function
#' @export
sm2 <- function(N=500, P=2*N, ...)
{
    Z1 <- as.matrix(scale(matrix(rnorm(N * P), N, P)))
    Z2 <- as.matrix(scale(matrix(rnorm(N * P), N, P)))
    knl <- within(list(),
    {
        EPS <- diag(N)
        ## KN1 <- tcrossprod(Z1) / P
        ## KN2 <- tcrossprod(Z2) / P
        KN1 <- gen.psd(N, P)
        KN2 <- gen.psd(N, P)
    })
    K <- length(knl)

    ## allowed kernels
    v <- knl[c('KN1', 'KN2')]
    
    ## true covariance
    ## X <- cbind(X00=rep(1, N))
    S <- with(knl, .5 * EPS + .8 * KN1 + 1.3 * KN2)
    y <- .mvn(1, 1.2, S)
    
    ## minque
    mdl <- mnq(y, v, NULL, ...)
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
    y <- .mvn(1, -1, S)

    ## minque
    ref <- within(list(),
    {
        par <- c(X00=-1, w)
        rpt <- vpd(y, v, x, par, ...)
    })

    ## kcv <- fcv(y, v, x, ...)
    
    ful <- mnq(y, v, x, zcp=1, rpt=1, ...)
    sel <- fwd(y, v, x, ...)
    list(sel=sel, ref=ref, ful=ful) #, kcv=kcv)
}
