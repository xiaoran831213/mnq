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
#' @param Y vector of response variable
#' @param K list of kernel matrices
#' @param X matrix of covariate
#' @param a initial guess for variance components
#'
#' @return list of results:
#'   * vcs: vector of variance component estimates
#'   * fix: vector of fixed coeficient estimates
rln_mnq <- function(Y, K, X=NULL, a=NULL)
{
    L <- length(K)
    N <- length(Y)
    
    ## initial weights
    if(is.null(a))
        a <- rep(1, L)

    ## V = a_0 K_0 + a_1 K_1 + ... + a_L K_L
    ## A = V^{-1}
    V <- Reduce(`+`, mapply(`*`, K, a, SIMPLIFY=FALSE))
    A <- chol2inv(chol(V)) # this may fail
    
    ## get P, Q = I - P, and R V_i and R y, where R = V^Q
    
    ## Q = V^{-1} - V^{-1}X (X'V^{-1}X)^{+} X'V^{-1}
    ## get Q K_0, Q K_1, ... Q K_L, and Qy
    Rv <- list()
    if(is.null(X) || length(X) == 0)
    {
        ## no X, Q = V^{-1} = A
        for(i in seq.int(K))
            Rv[[i]] <- A %*% K[[i]]
        Ry <- A %*% Y
    }
    else
    {
        ## Q = A - B (X'B )^{+} B'
        B <- A %*% X
        Q <- A - B %*% .ginv(t(X) %*% B) %*% t(B)
        for(i in seq.int(K))
            Rv[[i]] <- Q %*% K[[i]]
        Ry <- Q %*% Y
    }

    ## get y' Q K_l Q y = (Qy)' K_l (Qy), for l=0 ... L
    u <- double(L)
    for(i in seq(L))
        u[i] <- t(Ry) %*% K[[i]] %*% Ry  # y'Q K_l Qy

    ## get Tr(QK_l QK_m) = | (QK_l)*(QK_m)'|_f
    F <- matrix(.0, L, L)
    for(i in seq(L))
    {
        for(j in seq(i + 1, l= L - i))
        {
            F[i, j] <- sum(Rv[[i]] * t(Rv[[j]]))
            F[j, i] <- F[i, j]
        }
        F[i, i] <- sum(Rv[[i]] * t(Rv[[i]]))
    }

    ## solve:
    ## [Tr(QK_0 QK_0) ... Tr(QK_0 QK_L)]  [s_0^2] = [y'Q K_0 Q y]
    ## [Tr(QK_1 QK_0) ... Tr(QK_1 QK_L)]  [s_1^2] = [y'Q K_1 Q y]
    ## ...            ...            ...    ...          ...
    ## [Tr(QK_L QK_0) ... Tr(QK_L QK_L)]  [s_L^2] = [y'Q K_L Q y]
    ## for s_0^2, s_1^2, ... s_L^2
    w <- drop(.ginv(F) %*% u)

    ## GLS: beta = (X'V^{-1}X)^+ X'V^{-1}y
    if(length(X) > 0)
    {
        V <- Reduce(`+`, mapply(`*`, K, w, SIMPLIFY=FALSE))
        IX <- solve(V, X) # V^{-1}X
        b <- drop(.ginv(t(X) %*% IX) %*% t(IX) %*% Y)
        names(b) <- colnames(X)
    }
    else
        b <- NULL
    
    ## pack up
    list(vcs=w, fix=b, ini=a)
}

#' MINQUE 0
#'
#' A special case  of MINQUE, where the initial guess  for noise ($\alpha_0$) is
#' set to 1, and the rest be 0, thus,
#'
#' $V = a_0 K_0 + a_1... + a_L K_L = I$,
#'
#' which speed up computation.
#'
#' @param Y vector of obseved outcomes
#' @param K list of kernels, the first must be diagonal I(N x N).
#' @param X matrix of covariates
rln_mq0 <- function(Y, K, X=NULL)
{
    L <- length(K)
    N <- length(Y)

    ## sum of V_i, i = 1 .. K by initial weights
    ## V = 1 * K_1 + 0*K_2 ... 0*K_L = I  (since K_0 = I)

    ## Q = V^{-1} - V^{-1}X (X'V^{-1}X)^{+} X'V^{-1}
    ## get Q K_1, Q K_2, ... Q K_L, and Qy
    Rv <- list()
    if(is.null(X) || length(X) == 0)
    {
        ## No X, Q = V^{-1} = I => Q K_l = K_l, l = 1 ... L
        Rv <- K
        Ry <- Y # Q y = y
    }
    else
    {
        ## Q = V^{-1} - V^{-1}X (X'V^{-1}X)^+ X'V^{-1}
        ##   = I - X (X'X)^{+} X',   (V^{-1} = I)
        ##
        ## Q K = K - X (X'X)^{+} (X' K)     (K: K_1 ... K_L)
        ##
        ## Q <- diag(N) - X %*% .ginv(t(X) %*% X) %*% t(X)
        B <- X %*% .ginv(t(X) %*% X)
        for(l in seq(L))
        {
            ## Rv[[l]] <- Q %*% K[[l]] # QK_1 ... QK_L, Q K_0=Q
            Rv[[l]] <- K[[l]] - B %*% crossprod(X, K[[l]])
        }
        Ry <- Y - B %*% t(X) %*% Y 
        ## Ry <- Q %*% Y               # Qy
    }
    
    ## get y' Q K_l Q y = (Qy)' K_l (Qy), for l=0 ... L
    u <- double(L)
    for(l in seq(L))
        u[l] <- t(Ry) %*% K[[l]] %*% Ry  # y'Q K_l Qy

    ## get Tr(QK_l QK_m) = | (QK_l)*(QK_m)'|_f
    F <- matrix(.0, L, L)
    for(i in seq(L))
    {
        for(j in seq.int(i + 1, l = L - i))
        {
            F[i, j] <- sum(Rv[[i]] * t(Rv[[j]]))
            F[j, i] <- F[i, j]
        }
        F[i, i] <- sum(Rv[[i]] * t(Rv[[i]]))
    }

    ## solve:
    ## [Tr(QK_0 QK_0) ... Tr(QK_0 QK_L)]  [s_0^2] = [y'Q K_0 Q y]
    ## [Tr(QK_1 QK_0) ... Tr(QK_1 QK_L)]  [s_1^2] = [y'Q K_1 Q y]
    ## ...            ...            ...    ...          ...
    ## [Tr(QK_L QK_0) ... Tr(QK_L QK_L)]  [s_L^2] = [y'Q K_L Q y]
    ## for s_0^2, s_1^2, ... s_L^2
    w <- .ginv(F) %*% u 

    ## GLS for fixed effects
    ## beta = (X'V^{-1}X)^+ X'V^{-1}y
    ## approximate V with diag(V)
    if(length(X) > 0)
    {
        ## V <- Reduce("+", mapply("*", K, w, SIMPLIFY=FALSE))
        ## D <- solve(V, X) # V^{-1}X
        D <- X / drop(sapply(K, diag) %*% w)
        b <- drop(.ginv(t(X) %*% D) %*% t(D) %*% Y)
        names(b) <- colnames(X)
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
#' @param Y vector of response variable
#' @param K list of kernel matrices
#' @param X matrix of covariate
#' @param W initial values for variance components
#' @return list of results, includs:
#'   * vcs: vector of variance component estimates
#'   * fix: vector of fixed coeficient estimates
#'
#' @seealso \code{rln_mnq}
amo_mnq <- function(Y, K, X=NULL, W=NULL, ...)
{
    opt <- list(...)
    Y <- as.matrix(Y)
    X <- if(is.null(X)) matrix(0, 0, 0) else as.matrix(X)
    W <- if(is.null(W)) matrix(0, 0, 0) else as.matrix(W)
    .Call('amo_mnq', Y, K, X, W, opt, PACKAGE='mnq')
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
#'
#    * zcp enable zero capping? (def=0)
#'   * eht: halt on error (def=TRUE)
#'   * tol: convergence tolerence (def=1e-5)
#'   * itr: iterations allowed (def=Inf)
#'   * cpp: use C++ function for speed (def=TRUE)
#'   * quiet: TRUE to omit messages (def=FALSE)
#'   * zht: halt on negative variance component (def=!zcp)
#'
#'   * itc: prepend a intercept term to fixed terms? (def=TRUE)
#' 
#' @return list of results:
#' 
#'   * par: variance component esimates and fixed coefficients
#'   * rtm: running time
#'
#' @export
mnq <- function(y, v=NULL, x=NULL, w=NULL, ...)
{
    dot <- list(...)
    tol <- if(is.null(dot$tol))   1e-5 else dot$tol
    cpp <- if(is.null(dot$cpp))   TRUE else dot$cpp
    itr <- if(is.null(dot$itr))     20 else dot$itr
    vbs <- if(is.null(dot$quiet)) TRUE else !dot$quiet
    zcp <- if(is.null(dot$zcp))      0 else dot$zcp
    zht <- if(is.null(dot$zht))   !zcp else dot$zht
    eht <- if(is.null(dot$eht))   TRUE else dot$eht
    itc <- if(is.null(dot$itc))   TRUE else dot$itc
    rpt <- if(is.null(dot$rpt))  FALSE else dot$rpt

    ## append noise kernel
    N <- length(y)                      # sample size
    v <- c(list(EPS=diag(N)), v)        # first kernel is random noise
    K <- length(v)                      # kernel count

    ## sanity check: kernel size
    stopifnot(all(sapply(v, dim) == N))
    
    ## prepend intercept if necessary
    if(itc)
        x <- cbind(X00=rep(1.0, N), x)
    x <- as.matrix(x)
    y <- as.matrix(y)
    
    ## print('begin MINQUE')
    t0 <- Sys.time()

    ## initialize variance components
    if(is.null(w))
        w <- rep(1, K)          # equal by default
    vc0 <- c(w, rep(0, K))[1:K] # 0 for unspecified components
    vc1 <- c(1, rep(0, K - 1))

    ## MINQUE implementation
    imp <- if(cpp) amo_mnq else rln_mnq

    ## messages
    ps0 <- function(...) paste(..., collapse="")
    sp0 <- function(...) ps0(sprintf(...))
    ca0 <- function(...) if(vbs) cat(..., "\n", sep="")
    hdr <- ps0("ITR", sp0("%6s", names(v)), sp0("%7s", "MDF"))

    ## print the header
    ca0(hdr)

    dff <- max(abs(vc1 - vc0))
    ca0(ps0(sp0('%03d', itr), sp0("%6.3f", vc0), sp0("%7.1e", dff)))
    neg <- (vc0 < 0) * 1
    while(itr > 0)
    {
        ## call MINQUE core function, skip kernels with zero weights.
        rt0 <- imp(y, v, x, vc0)
        if(rt0$err && eht)
            stop("MNQ: halt on error: ", rt0$err)
        
        vc1 <- rt0$vcs
        dff <- max(abs(vc1 - vc0))

        ## report
        ca0(ps0(sp0('%03d', itr), sp0("%6.3f", vc1), sp0("%7.1e", dff)))
        
        ## convergence check
        if(dff < tol)
            break
        vc0 <- vc1

        ## halt on negative variance component?
        neg <- neg + (vc0 < 0)
        if(zht > 0 & max(neg) >= zht)
            break

        itr <- itr - 1 # next round
    }
    
    ## fixed effects
    fx0 <- drop(rt0$fix)

    ## pack up and return: estimates
    vcs <- drop(vc0)
    ## zero capping?
    if(zcp)
        vcs <- pmax(vcs, 0)
    fxs <- drop(fx0)
    names(vcs) <- names(v)
    names(fxs) <- colnames(x)
    par <- c(fxs, vcs)

    ## reports & return
    if(rpt)
        rpt <- list(rpt=vpd(y, v[-1], x, par, ...))
    else
        rpt <- NULL
    par <- list(par=par)

    td <- Sys.time() - t0; units(td) <- 'secs'; td <- as.numeric(td)
    ## print('end MINQUE')
    rtm <- list(rtm=c(rtm=td))

    c(par, rpt, rtm)
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
    ## 2*nlk <- eae + ldt + N*log(2*pi)
    nlk <- eae + ldt
    attributes(nlk) <- NULL
    aic <- 2 * sum(abs(w) > 1e-6) + nlk * N

    ## return
    rpt <- data.frame(N=N, hsq=hsq, SST=SST, nlk=nlk, aic=aic,
                      zeb=zeb, zel=zel, yeb=yeb, yel=yel,
                      ycb=ycb, ycl=ycl, zcb=zcb, zcl=zcl)
    rpt
}

#' Fixed Effect Getter
#'
#' by convention, fixed effect coefficients preceeds the
#' variance components, and the first variance component
#' is named 'EPS'.
#' 
#' @param x a vector of parameters
#' @return fixed effect coefficients
#' @export
fx <- function(x)
{
    i <- grep('eps', names(x), TRUE)
    if(length(i) > 0)
        x[seq(1, i - 1)]
    else
        x
}

#' Fixed Effect Setter
#'
#' by convention, fixed effect coefficients preceeds the
#' variance components, and the first variance component
#' is named 'EPS'.
#' 
#' @param x a vector of parameters
#' @param value to be assign to fixed effect coefficients
#' @export
`fx<-` <- function(x, value)
{
    i <- grep('eps', names(x), TRUE)
    if(length(i) > 0)
        i <- seq(1, i - 1)
    else
        i <- seq_along(x)
    x[i] <- rep(value, l=length(i))
    x
}

#' Variance Component Getter
#'
#' by convention, fixed effect coefficients preceeds the
#' variance components, and the first variance component
#' is named 'EPS'.
#' 
#' @param x a vector of parameters
#' @return variance components
#' @export
vc <- function(x)
{
    i <- grep('eps', names(x), TRUE)
    if(length(i) > 0)
        x[seq(i, length(x))]
    else
        NULL
}

#' Variance Component Setter
#'
#' by convention, fixed effect coefficients preceeds the
#' variance components, and the first variance component
#' is named 'EPS'.
#' 
#' @param x a vector of parameters
#' @param value to assign to variance components
#' @export
`vc<-` <- function(x, value)
{
    i <- grep('eps', names(x), TRUE)
    if(length(i) > 0)
        i <- seq(i, length(x))
    else
        i <- seq_along(x)
    x[i] <- rep(value, l=length(i))
    x
}
