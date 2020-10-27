#' Multivariant Normal Samples
#'
#' A simplified MASS::mvrnorm.
#'
#' @param N number of samples to draw
#' @param m vector of means, rotated if necessary
#' @param V matrix of covariance dimensions
#' @param drop TRUE to drop matrix of single sample to vector
#' @return matrix of N row samples and D column features, whose
#' covariance is \code{V}
#' @noRd
mvn <- function (N=1, m=0, V=NULL, drop=TRUE)
{
    ## default V is for demonstration
    if(is.null(V))
        V = rbind(c(1, .5, .25), c(.5, 1, .5), c(.25, .5, 1))

    ## dimensionality
    D <- nrow(V)

    ## mean vector
    m <- drop(rep(m, length=D))

    ## eigen decomposition
    e <- eigen(V, symmetric = TRUE)
    d <- e$values
    U <- e$vectors
    S <- sqrt(pmax(d, 0))          # square root of V

    ## random values
    y <- t(m + U %*% (S * matrix(rnorm(D * N), D, N)))

    if(drop)
        y <- drop(y)
    y
}

#' generate a matrix
#'
#' @param N number of row samples
#' @param P number of column features
#' @param ... additional arguments for \code{f}.
#' @return a matrix of N rows and P columns.
gen.mtx <- function(N=500, P=2*N, ...)
{
    x <- matrix(rnorm(N * P), N, P)
    rownames(x) <- sprintf("I%02X", seq(1, l=N))
    colnames(x) <- sprintf("X%02X", seq(1, l=P))
    x
}

#" standardize a matrix
std.mtx <- function(x)
{
    x - outer(rowMeans(x), colMeans(x), `+`) + mean(x)
}

#' generate PSD matrix
#'
#' @param N number of row samples
#' @param P number of column features
#' @param ... additional arguments for \code{f}.
#' @return a PSD matrix of N rows/columns.
gen.psd <- function(N=500, P=2*N, tau=NULL, ...)
{
    a <- tcrossprod(gen.mtx(N, P, ...)) / P
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
gen.knl <- function(N=500, P=N*2, L=2, ...)
{
    K <- replicate(L, gen.psd(N, P), simplify=FALSE)
    names(K) <- sprintf("K%02d", seq_along(K))
    K
}

#' Simulation 1
#'
#' No fixed effect
#' 
#' @param N sample size
#' @param P number of features
#' @param vcs variance components, noise included
#' @param fix fixed coefficients
#' @param ... additional arguments passed to \code{mtd}
#' @export
sim <- function(N=500, P=N, vcs=c(1.0, 0.5, 1.5), fix=c(1.0, 0.8), times=5e2, rs=NULL, ...)
{
    L <- length(vcs)
    M <- length(fix)

    res <- list()
    set.seed(rs)
    for(i in seq(times))
    {
        ## fixed effect
        X <- cbind(X00=rep(1, N), gen.mtx(N, M-1))
    
        ## random effect
        K <- c(list(EPS=diag(N)), gen.knl(N, P, L-1))  # kernels

        ## generate outcome
        y <- as.matrix(mvn(M=X %*% fix, V=vsm(K, vcs)))
    
        ## minque
        ## Rprof(interval=1e-3)
        r <- mnq(y, K, X, sst=1, zcp=0, itc=0, eps=0, ...)
        r <- with(r, data.frame(rtm=rtm, pvl=stt$P))
        
        res[[i]] <- cbind(N=N, P=P, r)
    }
    set.seed(NULL)
    res <- do.call(rbind, res)
    rownames(res) <- NULL
    res
}


bmk <- function(N=500, P=N, vcs=c(1, 1, 1), fix=NULL, itc=NULL, mtd=rln_mnq, times=100)
{
    L <- length(vcs)
    M <- length(fix)
    C <- length(itc)

    ## fixed effect
    X <- cbind(X00=rep(1, N * C), gen.mtx(N, M)) # covariants
    m <- 0
    if(length(X) > 0)
        m <- X %*% c(itc, fix)

    ## kernel (random effect)
    K <- c(list(EPS=diag(N)), gen.knl(N, P, L-1))  # kernels

    ## generate outcome
    y <- as.matrix(.mvn(M=m, V=vsm(K, vcs)))
    
    ## minque
    rtm <- replicate(times,
    {
        tm <- mtd(y, K, X)[['rtm']]
        tm["user.self"]
    })
    mean(rtm)
}

pow <- function(res, ...)
{
    val <- c('rtm', 'pvl')
    grp <- setdiff(names(res), val)

    res <- by(res, res[, grp], function(.)
    {
        cfg <- .[1, grp]
        rtm <- mean(.$rtm)
        rej <- mean(.$pvl < 0.05)
        data.frame(cfg, rep=nrow(.), rtm=rtm, pow=rej)
    })
    res <- do.call(rbind, res)
    res
}
