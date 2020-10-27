#' MINQUE-0
#'
#' A special case of MINQUE, where  the foremost kernel \eqn{K_1} is an identity
#' kernel matrix that models noise, and the initial guess of variance components
#' are 1  for the  foremost noise kernel  (\eqn{a_1=1}), and 0  for the  rest of
#' the kernels (\eqn{a_l=0, l=2 \dots L}).
#'
#' In this case, the initial covariance is simplifed to an identity matirx
#'
#' \deqn{V = a_1 K_1 + a_2 K_2 ... + a_L K_L = I_{N \times N}},
#'
#' which  speeds   up  MINQUE   considerably  since  the   spectral  decompotion
#' (inversion) of identity matrix takes no time.
#'
#' @param Y vector of the dependent variable;
#' @param K list of the kernels, the first must be identity I(N x N);
#' @param X matrix of the covariates.
#'
#' MINQUE-0 requires the first kernel to be identity \eqn{K_1 = I_{N \times N}}.
#' Since the  initial guess of variance  components is now fixed  to \eqn{(1, 0,
#' ... 0)}, the paramter \code{a} in standard MINQUE is dropped.
#' @seealso rln_mnq
mq0 <- function(Y, K, X=NULL)
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
    w <- drop(inv(F) %*% u )

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
    toc <- proc.time()

    ## pack up
    list(vcs=w, fix=b, rtm=toc - tic, Y=Y, K=K, X=X)
}

#' MINQUE-0 by Page
#'
#' Perform MINQUE-0, but  scan the input kernels several lines  at a time, which
#' is meant to process out-of-memory kernels.
#'
#' At this testing stage, the kernel is actually not out-of-memory.
#'
#'
#' @param Y vector of the dependent variable;
#' @param K list of the kernels, the first must be identity I(N x N);
#' @param X matrix of the covariates.
#' @param R the number of pages (def=1).
#'
#' MINQUE-0 requires the first kernel to be identity \eqn{K_1 = I_{N \times N}}.
#' Since the  initial guess of variance  components is now fixed  to \eqn{(1, 0,
#' ... 0)}, the paramter \code{a} in standard MINQUE is dropped.
#' @seealso mq0
mqp <- function(Y, K, X=NULL, R=1, diag=TRUE)
{
    tic <- proc.time()
    L <- length(K) # number of kernels
    N <- length(Y) # number of samples

    u <- double(L)
    T <- matrix(0, L, L)

    B <- tcrossprod(inv(crossprod(X)), X) # B <- (X'X)^{+}X'
    QY <- Y - crossprod(B, crossprod(X, Y))  # Qy = y - X(X'X)^{+}(X'y)
    for(r in seq(1, length=R))
    {
        . <- seq(1 + ceiling((r - 1) * N / R), ceiling(r * N / R))
        for(l in seq(1, L))
        {
            Klr <- K[[l]][., ]           # K_l^{(r)}
            KlQ <- Klr - Klr %*% X %*% B # (K_l Q)^{(r)}
            for(m in seq(l, L))
            {
                if(m > l)
                {
                    Kmr <- K[[m]][., ]           # K_m^{(r)}
                    KmQ <- Kmr - Kmr %*% X %*% B # (K_m Q)^{(r)}
                }
                else
                {
                    Kmr <- Klr
                    KmQ <- KlQ
                }
                T[l, m] <- T[l, m] + sum(KlQ * KmQ)
            }
            u[l] <- u[l] + crossprod(QY[.], Klr %*% QY)
        }
    }
    T[lower.tri(T)] <- t(T)[lower.tri(T)]

    ## solve:
    ## [Tr(QK_0 QK_0) ... Tr(QK_0 QK_L)]  [s_0^2] = [y'Q K_0 Q y]
    ## [Tr(QK_1 QK_0) ... Tr(QK_1 QK_L)]  [s_1^2] = [y'Q K_1 Q y]
    ## ...            ...            ...    ...          ...
    ## [Tr(QK_L QK_0) ... Tr(QK_L QK_L)]  [s_L^2] = [y'Q K_L Q y]
    ## for s_0^2, s_1^2, ... s_L^2
    w <- drop(inv(T) %*% u)

    ## GLS for fixed effects
    ## beta = (X'V^{-1}X)^+ X'V^{-1}y
    ## approximate V with diag(V)
    if(length(X) > 0)
    {
        ## V <- Reduce("+", mapply("*", K, w, SIMPLIFY=FALSE))
        ## D <- solve(V, X) # V^{-1}X
        if(diag)
        {
            b <- X / drop(sapply(K, diag) %*% w)
        }
        else
        {
            V <- Reduce(`+`, mapply(`*`, K, w, SIMPLIFY=FALSE))
            b <- solve(V, X) # V^{-1}X
        }
        b <- drop(inv(t(X) %*% b) %*% t(b) %*% Y)
        names(b) <- colnames(X)
    }
    else
        b <- NULL

    ## timer stop
    toc <- proc.time()

    ## pack up
    list(vcs=w, fix=b, rtm=toc - tic, Y=Y, K=K, X=X)
}

