## model selection functions

#' Forward Selection for MINQUE
#'
#' build a parsimonious model by adding subsequently the most
#' predictive kernel into it, while preserving its stability.
#' 
#' @param y vector of response variable
#' @param v list of covariance kernels
#' @param w vector of baseline parameters (beta and sigma^2),
#' where the VCs will always be preserved.
#' @param x matrix of covariat
#' @param ... additional arguments
#'   * tlr: threshold of likelihood ratio, in p-value.
#'   * msl: maximum number of selections
#'   * sby: select by (NLK, MSE, etc.)
#'   * vld: validation data set - a list of \code{y} and \code{v},
#'   if null, the validation is done on the training set.
#' @return the selected model:
#'   * par: parameter estimates
#'   * rpt: performance report
#' @export
fwd <- function(y, v, x=NULL, w=NULL, ...)
{
    dot <- list(...)
    msl <- if(is.null(dot$msl)) Inf   else dot$msl
    tlr <- if(is.null(dot$tlr)) 0.05  else dot$tlr
    sby <- if(is.null(dot$sby)) 'nlk' else dot$sby
    vld <- if(is.null(dot$vld)) list(y=y, v=v, x=x) else dot$vld

    N <- length(y)                      # sample size
    K <- length(v)                      # kernel count
    L <- if(is.null(x)) 0 else NCOL(x)  # covariate
    TLR <- stats::qchisq(1.0 - tlr, 0.5)

    ## initialization
    ksl <- intersect(names(v), names(w)) # ini selection
    kpl <- setdiff(names(v), ksl)        # ini kernel pool
    par <- NULL                          # ini par
    if(!is.null(w))
    {
        nms <- unique(c('X00', colnames(x), 'EPS', ksl))
        par <- rep(0, length(nms))
        names(par) <- nms
        par[names(w)] <- w
    }

    ## baseline model
    par <- mnq(y, v[ksl], x, par)$par
    vc(par) <- pmax(0, vc(par))
    rpt <- with(vld, vpd(y, v[ksl], x, par, ...))
    
    ## model selection
    while(length(kpl) > 0 && length(ksl) < msl)
    {
        err <- rpt[[sby]]
        nlk <- rpt$nlk
        ## try grow the model by one kernel
        mds <- lapply(kpl, function(n)
        {
            ksl <- c(ksl, n)
            par <- c(par, 0)
            par <- mnq(y, v[ksl], x, par, ...)$par
            vcs <- vc(par)
            rpt <- with(vld, vpd(y, v[ksl], x, par, ...))

            ## reject the growth if the new model is worse
            if(rpt[[sby]] > err || par[n] < 0 || par['EPS'] < 0)
                return(NULL)

            ## if any VC for other than EPS and new kernel is negative
            nps <- which(vcs[-c(1, length(vcs))] < 0)
            for(i in rev(nps))
            {
                ## drop the least predictive one and re-estimate
                par <- mnq(y, v[ksl[-i]], x, ...)$par
                vcs <- vc(par)
                rpt <- with(vld, vpd(y, v[ksl[-i]], x, par, ...))
                if(rpt[[sby]] < err && all(vcs > 0))
                {
                    nps <- integer()
                    ksl <- ksl[-1]
                    break
                }
            }

            ## discard if the adjusted model is STILL worse
            if(length(nps) > 0)
                return(NULL)
            list(par=par, rpt=rpt, ksl=ksl)
        })
        mds <- mds[!sapply(mds, is.null)]

        ## some condidate growth to choose from?
        if(length(mds) < 1)
            break
        ## rank the new models, select the best
        i <- order(sapply(mds, function(.) .$rpt[[sby]]))[1]
        m <- mds[[i]]
        par <- m$par
        rpt <- m$rpt
        ksl <- m$ksl

        ## cross the threshold of likelihood?
        LRT <- N * (nlk - rpt$nlk)
        cat(sprintf("LRT=%7f, TLR=%7f\n", LRT, TLR))
        if(LRT < TLR)
            break

        ## update kernel pool
        kpl <- setdiff(kpl, ksl)
    }

    rpt$nsl <- length(vc(par))
    list(par=par, rpt=rpt)
}


#' K-fold CV for Forward Selection
#'
#' use k-fold CV to decide the proper number of variance
#' component (\code{msl}) for forward selection (\code{fwd}).
#' 
#' @param y vector of response variable
#' @param v list of covariance kernels
#' @param x matrix of covariat
#' @param k number of folds (def=5)
#' @param ... additional arguments to pass to \code{fwd}
#' 
#' @return the selected model:
#'   * par: parameter estimates
#'   * rpt: performance report
#' @export
fcv <- function(y, v=NULL, x=NULL, k=5, ...)
{
    dot <- list(...)
    N <- NROW(y)

    ## k-fold mask
    f <- sample(rep(seq(k), length=N))

    mds <- list()                       # models
    for(i in seq(k))
    {
        cat(sprintf("CV: %02d\n", i))
        m <- f != i                     # training
        yi <- y[m]
        vi <- lapply(v, `[`, m, m)
        xi <- x[m, , drop=FALSE]

        n <- f == i                     # validation
        vld <- within(list(),
        {
            y <- y[n]
            v <- lapply(v, `[`, n, n)
            x <- x[n, , drop=FALSE]
        })

        ## model by foward selection
        mds[[i]] <- fwd(yi, vi, xi, vld=vld, ...)
    }

    ## decide maximum selection
    par <- lapply(mds, `[[`, 'par')
    vcs <- lapply(par, vc)
    rpt <- do.call(rbind, lapply(mds, `[[`, 'rpt'))
    len <- sapply(vcs, length)
    msl <- round(mean(sapply(vcs, length)))
    
    ## selection by complete data
    fwd(y, v, x, msl=msl, ...)
}
