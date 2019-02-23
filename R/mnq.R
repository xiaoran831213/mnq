## Modified MINQUE

#' Generalized Inverse
#'
#' A simplified copy of MASS::ginv
#'
#' @param x matrix to be inverted
#' @return generalized inverse of \code{x}
.ginv <- function(x)
{
    with(svd(x),
    {
        i <- d > d[1L] * sqrt(.Machine$double.eps)
        if(all(i))
            v %*% (t(u) / d)
        else
            v[, i, drop=FALSE] %*% (t(u[, i, drop=FALSE]) / d[i])
    })
}

#' Weight Sum of Kernels
#'
#' @param v list of  kernels
#' @param w vector of weights
#' @return sum of \code{v} weighted by \code{w}
.vsum <- function(v, w)
{
    Reduce(`+`, mapply(`*`, v, w, SIMPLIFY=FALSE))
}

#' R-Language MINQUE
#'
#' Given a model specified by a response variable \code{y}, a series
#' of K kernels \code{V_1 ... V_K}, and an optional matrix of
#' covariate, solve the variance components and fixed coefficients.
#'
#' @param y. vector of response variable
#' @param v. list of kernel matrices
#' @param x. matrix of covariate
#' @param w. initial values for variance components
#'
#' @return list of results:
#'   * vcs: vector of variance component estimates
#'   * fix: vector of fixed coeficient estimates
rln_mnq <- function(y., v., x.=NULL, w.=NULL)
{
    K <- length(v.)

    ## initial weights
    if(is.null(w.))
        w. <- rep(1/K, K)

    ## sum of V_i, i = 1 .. K by initial weights
    V <- .vsum(v., w.)
    A <- try(chol2inv(chol(V)), silent=TRUE)
    if(inherits(A, 'try-error'))
       A <- .ginv(V)
    
    ## get P, Q = I - P, and R V_i and R y, where R = V^Q
    Rv <- list()
    if(is.null(x.))
    {
        ## no X, P = 0, Q = I, R V_i = V^ I V_i = V^V_i
        for(i in seq.int(K))
            Rv[[i]] <- A %*% v.[[i]]
        Ry <- A %*% y.
    }
    else
    {
        ## P = X (X'V^X)^ (V^X)'
        B <- solve(V, x.)                # V^X
        C <- .ginv(t(x.) %*% B) %*% t(B) # (X'V^X)^ (V^X)'
        P <- x. %*% C

        ## R V_i = V^ (I - P) V_i = V^ (V_i - P V_i)
        for(i in seq.int(K))
            Rv[[i]] <- A %*% (v.[[i]] - P %*% v.[[i]])
        Ry <- A %*% (y. - P %*% y.)
    }

    ## u_i = e' V_i e = y' R V_i R y
    u<- double(K)
    for(i in seq.int(K))
        u[i] <- sum(Ry * v.[[i]] %*% Ry)

    ## Caculate F: F_ij = Tr(R V_i R V_j)
    F <- matrix(.0, K, K)
    for(i in seq.int(K))
    {
        for(j in seq.int(K)[-1L])
        {
            F[i, j] <- sum(Rv[[i]] * Rv[[j]])
            F[j, i] <- F[i, j]
        }
        F[i, i] <- sum(Rv[[i]] * Rv[[i]])
    }

    ## [s2_1, .., s2_K]' = [yA1y, .., yAKy]' = solve(F, u) = F^u
    w <- .ginv(F) %*% u 
    ##  bypass A1, .. A_K in (Rao. 1971)

    ## GLS for fixed effects
    if(!is.null(x.))
    {
        b <- C %*% y.
        names(b) <- colnames(x.)
    }
    else
        b <- NULL

    ## pack up
    list(vcs=w, fix=b)
}

#' Cpp MINQUE
#'
#' Given a model specified by a response variable \code{y}, a series
#' of K kernels \code{V_1 ... V_K}, and an optional matrix of
#' covariate, solve the variance components and fixed coefficients.
#'
#' The function is implemented by RcppArmadillo, slightly faster then
#' R-language.
#'
#' @param y. vector of response variable
#' @param v. list of kernel matrices
#' @param x. matrix of covariate
#' @param w. initial values for variance components
#' @return list of results, includs:
#'   * vcs: vector of variance component estimates
#'   * fix: vector of fixed coeficient estimates
#'
#' @seealso \code{rln_mnq}
amo_mnq <- function(y., v., x.=NULL, w.=NULL)
{
    y. <- as.matrix(y.)
    x. <- if(is.null(x.)) matrix(0, 0, 0) else as.matrix(x.)
    w. <- if(is.null(w.)) matrix(0, 0, 0) else as.matrix(w.)
    .Call('amo_mnq', y., v., x., w., PACKAGE='mnq')
}


#' MINQUE
#'
#' Main function of MINQUE, allowing iteration.
#' 
#' @param y vector of dependent variable
#' @param v list of covariance kernels
#' @param x matrix of covariat
#' @param w vector of initial paramters, fixed effect by followed
#' by variance components
#' 
#' @param ... additional options
#'   * tol: convergence tolerence (def=1e-5)
#'   * itr: iterations allowed (def=Inf)
#'   * cpp: use C++ function for speed (def=TRUE)
#'   * quiet: TRUE to omit messages (def=FALSE)
#' 
#' @return list of results:
#'   * par: variance component esimates and fixed coefficients
#'   * rtm: running time
#'
#' @export
mnq <- function(y, v=NULL, x=NULL, w=NULL, ...)
{
    dot <- list(...)
    tol <- if(is.null(dot$tol)) 1e-5 else dot$tol
    cpp <- if(is.null(dot$cpp)) TRUE else dot$cpp
    itr <- if(is.null(dot$itr))   20 else dot$itr
    vbs <- if(is.null(dot$quiet)) TRUE else !dot$quiet
    rpt <- if(is.null(dot$rpt)) FALSE else dot$rpt
    
    ## append noise kernel
    N <- length(y)                      # sample size
    v <- c(list(EPS=diag(N)), v)
    K <- length(v)                      # kernel count

    ## prepend intercept if necessary
    if(is.null(x) || !any(grepl('X00', colnames(x)[1], TRUE)))
        x <- cbind(X00=rep(1.0, N), x)
    x <- as.matrix(x)
    y <- as.matrix(y)

    ## print('begin MINQUE')
    t0 <- Sys.time()

    ## initialize variance components
    .m <- lm(y ~ x - 1)
    .s <- sum(residuals(.m)^2) / (N - NCOL(x))
    if(is.null(w))
        vc0 <- rep(.s / K, K)
    else
        vc0 <- w[seq(NCOL(x) + 1, length(w))]
    vc1 <- c(.s, rep(0, K - 1))
    rm(.m, .s)

    ## MINQUE implementation
    .mq <- if(cpp) amo_mnq else rln_mnq
    
    ## messages
    ps0 <- function(...) paste(..., collapse="")
    sp0 <- function(...) ps0(sprintf(...))
    ca0 <- function(...) if(vbs) cat(..., "\n", sep="")
    hdr <- ps0("ITR", sp0("%6s", names(v)), sp0("%7s", "MDF"))

    ## print the header
    ca0(hdr)

    dff <- max(abs(vc1 - vc0))
    ca0(ps0(sp0('%03d', itr), sp0("%6.3f", vc0), sp0("%7.1e", dff)))
    while(itr > 0)
    {
        ## call MINQUE core function, skip kernels with zero weights.
        rt0 <- .mq(y, v, x, vc0)

        vc1 <- rt0$vcs
        dff <- max(abs(vc1 - vc0))

        ## report
        ca0(ps0(sp0('%03d', itr), sp0("%6.3f", vc1), sp0("%7.1e", dff)))
        
        ## convergence check
        if(dff < tol)
            break
        vc0 <- vc1
        itr <- itr - 1
    }

    ## fixed effects
    fx0 <- rt0$fix

    td <- Sys.time() - t0; units(td) <- 'secs'; td <- as.numeric(td)
    ## print('end MINQUE')

    ## pack up and return: estimates
    names(vc0) <- names(v)
    names(fx0) <- colnames(x)
    par <- c(fx0, vc0)

    ## reports & return
    ret=list(par=par, rtm=td)
    if(rpt)
    {
        ret$rpt <- vpd(y, v[-1], x, par, ...)
    }
    ret
}

#' Variance component model prediction
#' 
#' @param y vector of response variable
#' @param v list of covariance kernels
#' @param w vector of parameters (beta, and sigma^2)
#' @param x matrix of covariate
#' @param ... additional arguments
vpd <- function(y, v=NULL, x=NULL, w, ...)
{
    dot <- list(...)
    y <- unname(y)
    N <- NROW(y)
    Q <- length(w)

    ## kernels, prepended with noise
    v <- c(list(EPS=diag(N)), v)
    K <- length(v)

    ## variance components, and fixed coeficients
    vcs <- w[seq(Q - K + 1, Q)]
    fix <- w[seq(1, Q - K)]
    
    ## matrix of fixed effect, prepended with intercept
    M <- if(is.null(x)) 0 else NCOL(x)
    if(length(fix) == M + 1)
    {
        x <- cbind(X00=rep(1, N), x)
        M <- M + 1
    }
    
    xb <- if(M > 0) x %*% fix else 0    # x beta
    e <- y - xb                         # fix effect residual
    SST <- sum(e^2) / (N - M)           # sum square total

    ## heritability
    hsq <- unname(1 - vcs[1] / SST)      # hsq 2
    ## h2c <- unname(1 - vcs[1] / sum(vcs)) # hsq 3
    
    ## predicted random effects
    V <- unname(.vsum(v, vcs))
    A <- try(chol2inv(chol(V)), silent=TRUE)
    if(inherits(A, 'try-error'))
       A <- .ginv(V)
    Ae <- A %*% e
    
    zub <- (V - diag(vcs[1], N)) %*% Ae   # use BLUP
    zul <- e - Ae / diag(A)               # LOO
    yhb <- zub + xb
    yhl <- zul + xb

    zeb <- mean((e - zub)^2)     # e MSE 1, BLUP
    zel <- mean((e - zul)^2)     # e MSE 2, LOOV
    yeb <- mean((y - yhb)^2)     # y MSE 1, BLUP
    yel <- mean((y - yhl)^2)     # y MSE 2, LOOV
    
    ## correlation between truth and prediction
    ycb <- if(abs(sd(yhb)) < 1e-8) 0 else cor(y, yhb)
    ycl <- if(abs(sd(zul)) < 1e-8) 0 else cor(y, yhl)
    zcb <- if(abs(sd(zub)) < 1e-8) 0 else cor(e, zub)
    zcl <- if(abs(sd(zul)) < 1e-8) 0 else cor(e, zul)

    ## negative likelihood
    ldt <- with(determinant(V), modulus * sign) / N
    eae <- sum(Ae * e) / N    # e^T V^{-1} e
    ## nlk <- .5 * (eae + ldt + log(2 * pi))
    nlk <- eae + ldt
    attributes(nlk) <- NULL
    
    ## return
    rpt <- data.frame(N=N, hsq=hsq, SST=SST, nlk=nlk,
                      zeb=zeb, zel=zel, yeb=yeb, yel=yel,
                      ycb=ycb, ycl=ycl, zcb=zcb, zcl=zcl)
    rpt
}
