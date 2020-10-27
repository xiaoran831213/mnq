## Modified MINQUE

#' Generalized Inverse
#'
#' A simplified copy of MASS::ginv
#'
#' @param x positive definite matrix to be inverted
#' @return generalized inverse of \code{x}
#' @noRd
inv <- function(x)
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
#' @noRd
vsm <- function(v, w)  Reduce(`+`, mapply(`*`, v, w, SIMPLIFY=FALSE))

#' MINQUE in R language
#'
#' Fit  an linear  mixed model  (LMM) by  MINQUE, given  the vector  of response
#' \code{\bold{y}}, the list of L kernels \code{\bold{K}_1 ...  \bold{K}_L}, and
#' (optional) covariate(s) in \code{\bold{X}}.
#'
#' @param Y vector of the response variable;
#' @param K list of the L kernel matrices;
#' @param X matrix of the covariate(s);
#' @param a vector of the initial guess of variance components.
#'
#' @return a list of results and context: \itemize{
#'   \item{vcs: }{estimated variance components};
#'   \item{fix: }{estimated fixed effect coeficient};
#'   \item{rtm: }{running time};
#'   \item{err: }{error code}.
#'   \item{the input data: }{Y, K, X, and a}.}
#'
#' If MINQUE  fails because of  the weighted sum  of kernels being  not positive
#' semi definit, an  error code of 2  is assigned to \code{err}  in the returned
#' list, otherwise, \code{err} is set to 0.
#'
#' In most LMM the first kernel  \bold{K}_1 is a noise kernel (identity matrix),
#' MINQUE does not enforce such rule, but MINQUE-0 does.
#' 
#' @seealso \code{link{rln_mq0}} for the special high speed MINQUE-0.
rln_mnq <- function(Y, K, X=NULL, a=NULL, ...)
{
    tm <- Sys.time() # timer start

    L <- length(K)
    N <- length(Y)
    
    ## initial weights
    if(is.null(a))
        a <- rep(1, L)

    ## V = a_0 K_0 + a_1 K_1 + ... + a_L K_L, A = V^{-1}
    V <- Reduce(`+`, mapply(`*`, K, a, SIMPLIFY=FALSE))
    A <- chol2inv(chol(V)) # this may fail
    
    ## get P, Q = I - P, and R V_i and R y, where R = V^Q
    ##
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
        Q <- A - B %*% inv(t(X) %*% B) %*% t(B)
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
    w <- drop(inv(F) %*% u)

    ## GLS: beta = (X'V^{-1}X)^+ X'V^{-1}y
    if(length(X) > 0)
    {
        V <- Reduce(`+`, mapply(`*`, K, w, SIMPLIFY=FALSE))
        IX <- solve(V, X) # V^{-1}X
        b <- drop(inv(t(X) %*% IX) %*% t(IX) %*% Y)
        names(b) <- colnames(X)
    }
    else
        b <- NULL
    
    ## timer stop
    tm <- Sys.time() - tm; units(tm) <- 'secs'; tm <- as.numeric(tm)

    ## pack up
    list(vcs=w, fix=b, rtm=tm, Y=Y, K=K, X=X, a=a)
}

#' MINQUE 0 in R language
#'
#' A special case of MINQUE, where  the foremost kernel \eqn{K_1} is an identity
#' kernel matrix that models noise, and the initial guess of variance components
#' are 1 for  the foremost kernel (\eqn{a_1=1}),  and 0 for the  rest of kernels
#' (\eqn{a_l=0, l=2 \dots L}).
#'
#' In this special case, the initial covariance is simplifed to identity matirx
#'
#' \deqn{V = a_1 K_1 + a_2 K_2 ... + a_L K_L = I_{N \times N}},
#'
#' which speeds up  the computation considerably since  the spectral decompotion
#' of identity matrix \eqn{I_{N \times N}} takes no time.
#'
#' @param Y vector of the dependent variable;
#' @param K list of the kernels, the first must be identity I(N x N);
#' @param X matrix of the covariates.
#'
#' MINQUE-0 requires the first kernel to be identity \eqn{K_1 = I_{N \times N}}.
#' Since the  initial guess of variance  components is now fixed  to \eqn{(1, 0,
#' ... 0)}, the paramter \code{a} in standard MINQUE is dropped.
#' @seealso rln_mnq
rln_mq0 <- function(Y, K, X=NULL)
{
    tic <- proc.time() # timer start

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
        ## Q <- diag(N) - X %*% inv(t(X) %*% X) %*% t(X)
        B <- X %*% inv(crossprod(X))
        for(l in seq(L))
        {
            ## Rv[[l]] <- Q %*% K[[l]] # QK_1 ... QK_L, Q K_0=Q
            Rv[[l]] <- K[[l]] - B %*% crossprod(X, K[[l]])
        }
        Ry <- Y - B %*% crossprod(X, Y)
        ## Ry <- Q %*% Y               # Qy
    }
    
    ## get y' Q K_l Q y = (Qy)' K_l (Qy), for l=0 ... L
    u <- double(L)
    for(l in seq(L))
        u[l] <- crossprod(Ry, K[[l]]) %*% Ry  # y'Q K_l Qy

    ## get Tr(QK_l QK_m) = | (QK_l)*(QK_m)'|_f
    F <- matrix(.0, L, L)
    for(i in seq(L))
    {
        for(j in seq.int(i + 1, l = L - i))
        {
            ## F[i, j] <- sum(Rv[[i]] * t(Rv[[j]]))
            F[i, j] <- sum(Rv[[i]] * Rv[[j]])
            F[j, i] <- F[i, j]
        }
        F[i, i] <- sum(Rv[[i]] * Rv[[i]])
    }

    ## solve:
    ## [Tr(QK_0 QK_0) ... Tr(QK_0 QK_L)]  [s_0^2] = [y'Q K_0 Q y]
    ## [Tr(QK_1 QK_0) ... Tr(QK_1 QK_L)]  [s_1^2] = [y'Q K_1 Q y]
    ## ...            ...            ...    ...          ...
    ## [Tr(QK_L QK_0) ... Tr(QK_L QK_L)]  [s_L^2] = [y'Q K_L Q y]
    ## for s_0^2, s_1^2, ... s_L^2
    w <- inv(F) %*% u 

    ## GLS for fixed effects
    ## beta = (X'V^{-1}X)^+ X'V^{-1}y
    ## approximate V with diag(V)
    if(length(X) > 0)
    {
        ## V <- Reduce("+", mapply("*", K, w, SIMPLIFY=FALSE))
        ## D <- solve(V, X) # V^{-1}X
        D <- X / drop(sapply(K, diag) %*% w)
        b <- drop(inv(t(X) %*% D) %*% t(D) %*% Y)
        names(b) <- colnames(X)
    }
    else
        b <- NULL

    ## timer stop
    tm <- proc.time() - tic

    ## pack up
    list(vcs=w, fix=b, rtm=tm, Y=Y, K=K, X=X)
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
#' @param Y vector of the response variable;
#' @param K list of the kernel matrices;
#' @param X matrix of the covariate(s);
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
#'   * zcp: enable zero capping? (def=0)
#'   * stt: compute test statistics?
#'   * err: halt on error (def=TRUE)
#'   * tol: convergence tolerence (def=1e-5)
#'   * itr: iterations allowed (def=Inf)
#'   * cpp: use C++ function for speed (def=TRUE)
#'   * zht: halt on negative variance component (def=!zcp)
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
    eps <- dot$eps %||% TRUE # need noise kernel?
    itc <- dot$itc %||% TRUE # need intercept?
    tol <- dot$tol %||% 1e-4 # convergence
    cpp <- dot$cpp %||% FALSE
    itr <- dot$itr %||% 0x10
    zcp <- if(is.null(dot$zcp)) TRUE else dot$zcp
    zht <- if(is.null(dot$zht)) !zcp else dot$zht
    stt <- dot$rpt %||% TRUE
    err <- dot$err %||% TRUE
    rpt <- dot$rpt %||% FALSE
    par <- dot$par %||% TRUE

    ## append noise kernel
    N <- length(y)                      # sample size
    if(eps)
        v <- c(list(EPS=diag(N)), v) # first kernel is random noise
    K <- length(v)                      # kernel count

    ## sanity check: kernel size
    stopifnot(all(sapply(v, dim) == N))
    
    ## prepend intercept if necessary
    if(itc)
        x <- cbind(X00=rep(1.0, N), x)
    x <- as.matrix(x)
    y <- as.matrix(y)

    ## print('begin MINQUE')
    tic <- Sys.time()

    ## initialize variance components
    if(is.null(w))
    {
        rt0 <- mq0(y, v, x) # equal by default
        w <- rt0$vcs
    }
    vc0 <- c(w, rep(0, K))[1:K] # 0 for unspecified components
    vc1 <- c(1, rep(0, K - 1))

    ## MINQUE implementation
    imp <- if(cpp) amo_mnq else rln_mnq

    ## messages
    ps0 <- function(...) paste(..., collapse="")
    sp0 <- function(...) ps0(sprintf(...))
    ca0 <- function(...) cat(..., "\n", sep="")
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
        if(!is.null(rt0$err) && err)
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
    fix <- drop(rt0$fix)     # fixed effects
    vcs <- drop(vc0)         # random effects
    if(zcp)                  # zero capping?
        vcs <- pmax(vcs, 0)
    names(vcs) <- names(v)
    names(fix) <- colnames(x)

    ## reports & return
    if(rpt)
        rpt <- list(rpt=vpd(y, v[-1], x, par, ...))
    else
        rpt <- NULL

    if(par)
        par <- list(par=c(fix, vcs))
    else
        par <- NULL

    ## test
    if(stt)
        stt <- list(stt=qst(y, v[-1], x, fix, vcs[-1]))
    else
        stt <- NULL
    toc <- Sys.time()

    rtm <- toc - tic; units(rtm) <- 'secs'; rtm <- as.numeric(rtm)
    ## print('end MINQUE')
    rtm <- list(rtm=c(rtm=rtm))
    
    c(par, rpt, rtm, stt)
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
    V <- unname(vsm(v, vcs))
    A <- try(chol2inv(chol(V)), silent=TRUE)
    if(inherits(A, 'try-error'))
       A <- inv(V)
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
