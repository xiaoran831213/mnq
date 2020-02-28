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
#' @param zht halt on negative esimates.
#' @param ... additional arguments
#' 
#'   * tlr: threshold of likelihood ratio, in p-value.
#'   * nks: number of kernels to be selected.
#'   * vld: validation data set - a list of \code{y} and \code{v},
#'   if null, the validation is done on the training set.
#' @return the selected model:
#'   * par: parameter estimates
#'   * rpt: performance report
#' @export
fwd <- function(y, v, x=NULL, w=NULL, zht=2, ...)
{
    dot <- list(...)
    ## nkn <- if(is.null(dot$nkn)) Inf   else dot$nkn
    tlr <- if(is.null(dot$tlr)) 1.0 else dot$tlr
    nks <- if(is.null(dot$nks)) Inf else dot$nks
    vld <- if(is.null(dot$vld)) list(y=y, v=v, x=x) else dot$vld
    
    N <- length(y)                      # sample size
    K <- length(v)                      # kernel count
    L <- if(is.null(x)) 0 else NCOL(x)  # covariate
    TLR <- stats::qchisq(tlr, 0.5, lower.tail=FALSE)
    MLR <- Inf
    
    ## initialization
    ksl <- intersect(names(v), names(w)) # ini selection
    kpl <- setdiff(names(v), ksl)        # ini kernel pool
    par <- NULL                          # ini par
    if(!is.null(w))
    {
        nms <- unique(c('X00', colnames(x), 'EPS', ksl))
        par <- rep(0, length(nms))
        names(par) <- nms
        nms <- intersect(nms, names(w))
        par[nms] <- w[nms]
    }

    ## print('begin MINQUE')
    t0 <- Sys.time()

    ## baseline model
    par <- mnq(y, v[ksl], x, par)$par
    vc(par) <- pmax(0, vc(par))
    rpt <- with(vld, vpd(y, v[ksl], x, par, ...))
    
    ## model selection
    while(length(kpl) > 0  && length(ksl) < nks)
    {
        nlk <- rpt$nlk
        aic <- rpt$aic
        yel <- rpt$yel
        ycl <- rpt$ycl
        ## try grow the model by one kernel
        mds <- lapply(kpl, function(n)
        {
            ksl <- c(ksl, n)
            par <- c(par, 0)
            par <- mnq(y, v[ksl], x, par, zht=zht, ...)$par
            vcs <- vc(par)
            rpt <- with(vld, vpd(y, v[ksl], x, par, ...))
            lrt <- NROW(vld$y) * (nlk - rpt$nlk) # LRT
            dic <- aic - rpt$aic
            cat(sprintf("LRT=%7f, DIC=%7f\n", lrt, dic))

            ## reject the growth if the new model is worse
            ## if(lrt < 0 || par[n] < 0 || par['EPS'] < 0)
            if(dic < 0 || par[n] < 0 || par['EPS'] < 0)
                return(NULL)

            ## if any VC for other than EPS and new kernel is negative
            nps <- which(vcs[-c(1, length(vcs))] < 0)
            for(i in rev(nps))
            {
                ## drop the least predictive one and re-estimate
                par <- mnq(y, v[ksl[-i]], x, ...)$par
                vcs <- vc(par)
                rpt <- with(vld, vpd(y, v[ksl[-i]], x, par, ...))
                lrt <- NROW(vld$y) * (nlk - rpt$nlk) # LRT
                dic <- aic - rpt$aic
                cat(sprintf("LRT=%7f, DIC=%7f\n", lrt, dic))
                ## find a better model
                ## if(lrt > 0 && all(vcs > 0))
                if(dic > 0 && all(vcs > 0))
                {
                    nps <- integer()
                    ksl <- ksl[-i]
                    break
                }
            }

            ## discard if the adjusted model is STILL worse
            if(length(nps) > 0)
                return(NULL)
            list(par=par, rpt=rpt, ksl=ksl)
        })
        mds <- mds[!sapply(mds, is.null)]

        ## any condidate growth to choose from?
        if(length(mds) < 1)
        {
            cat("FWD  HALT\n")
            break
        }
        ## rank the new models, select the best
        i <- order(sapply(mds, function(.) .$rpt$nlk))[1]
        m <- mds[[i]]

        ## cross the threshold of likelihood ratio test?
        lrt <- NROW(vld$y) * (nlk - m$rpt$nlk)
        dic <- aic - m$rpt$aic
        cat(sprintf("SL?: LRT=%7f, DIC=%7f, TLR=%7f\n", lrt, dic, TLR))
        if(lrt < TLR)
            break

        ## remember the minimum LRT encountered
        MLR <- min(lrt, MLR)

        ## update kernel pool
        par <- m$par
        rpt <- m$rpt
        ksl <- m$ksl
        kpl <- setdiff(kpl, ksl)
    }

    td <- Sys.time() - t0; units(td) <- 'secs'; td <- as.numeric(td)
    ## print('end MINQUE')

    rpt$nks <- length(vc(par)) - 1
    ## minimum threshold required to include all selected kernels
    rpt$mlr <- stats::pchisq(MLR, 0.5, lower.tail=FALSE)
    list(par=par, rpt=rpt, rtm=c(rtm=td))
}


#' K-fold CV for Forward Selection
#'
#' use k-fold CV to decide the proper number of variance
#' component (\code{nks}) for forward selection (\code{fwd}).
#' 
#' @param y vector of response variable
#' @param v list of covariance kernels
#' @param x matrix of covariat
#' @param k number of folds (def=5)
#' @param ... additional arguments to pass to \code{fwd}
#'
#' @return the selected model:
#' 
#'   * par: parameter estimates
#'   * rpt: performance report
#' @export
fcv <- function(y, v=NULL, x=NULL, ncv=5, ...)
{
    dot <- list(...)
    set.seed(dot$seed)
    N <- NROW(y)

    ## k-fold mask
    f <- sample(rep(seq(ncv), length=N))

    ## print('begin MINQUE')
    t0 <- Sys.time()

    mds <- list()                       # models
    for(i in seq(ncv))
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

        ## foward selection
        mds[[i]] <- fwd(yi, vi, xi, vld=vld, tlr=1.0, ...)
    }

    ## decide maximum selection
    par <- lapply(mds, `[[`, 'par')
    nms <- Reduce(union, sapply(par, names))
    par <- do.call(rbind, lapply(par, `[`, nms))
    colnames(par) <- nms
    par[is.na(par)] <- 0
    rpt <- as.matrix(do.call(rbind, lapply(mds, `[[`, 'rpt')))
    
    ## weighted estimates
    wgt <- exp(-rpt[, "nlk"])
    par <- colSums(par * wgt) / sum(wgt)
    rpt <- colSums(rpt * wgt) / sum(wgt)

    ## selection by complete data
    mlr <- rpt['mlr']
    nks <- min(1, ceiling(rpt['nks']))
    cat("\nNKS=", nks, "\n", sep="")
    ret <- fwd(y, v, x, nks=nks, ...)
    
    ## pack up  and return
    td <- Sys.time() - t0; units(td) <- 'secs'; td <- as.numeric(td)
    ret$rtm <- c(rtm=td)
    ret
}
