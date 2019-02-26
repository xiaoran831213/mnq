## model selection functions

#' Forward Selection for MINQUE
#'
#' build a parsimonious model by adding subsequently the most
#' predictive kernel into it, while preserving its stability.
#' 
#' @param y vector of response variable
#' @param v list of covariance kernels
#' @param w vector of parameters (beta, and sigma^2)
#' @param x matrix of covariat
#' @param ... additional arguments
#'   * msl: maximum number of selections
#'   * sby: select by (NLK, MSE, etc.)
#' @return the selected model:
#'   * par: parameter estimates
#'   * rpt: performance report
#' @export
fwd <- function(y, v, x=NULL, w=NULL, ...)
{
    dot <- list(...)
    msl <- if(is.null(dot$msl)) Inf   else dot$msl
    sby <- if(is.null(dot$sby)) 'nlk' else dot$sby
    vld <- if(is.null(dot$vld)) list(y=y, v=v, x=x) else dot$vld
    N <- length(y)                      # sample size
    K <- length(v)                      # kernel count

    ## model selection
    kpl <- v                            # kernel pool
    ksl <- list()                       # selected
    par <- mnq(y, NULL, x)$par          # NULL model
    rpt <- vpd(y, NULL, x, par)
    while(length(kpl) > 0 && length(ksl) < msl)
    {
        err <- rpt[[sby]]
        ## grow the model
        mds <- lapply(names(kpl), function(n)
        {
            ksl <- c(ksl, kpl[n])
            par <- mnq(y, ksl, x, ...)$par
            vcs <- par[c("EPS", names(ksl))]
            rpt <- vpd(y, ksl, x, par, ...)

            ## reject the growth if the new model is worse
            if(rpt[[sby]] > err)
                return(NULL)

            ## if any VC for other than EPS and new kernel is negative
            nps <- which(vcs[-c(1, length(vcs))] < 0)
            for(i in rev(nps))
            {
                ## drop the least predictive one and re-estimate
                par <- mnq(y, ksl[-i], x, ...)$par
                vcs <- par[c("EPS", names(ksl[-i]))]
                rpt <- vpd(y, ksl[-i], x, par, ...)
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
        if(length(mds) > 0)
        {
            ## rank the new models, select the best
            idx <- order(sapply(mds, function(.) .$rpt[[sby]]))[1]
            mds <- mds[[idx]]
            par <- mds$par
            rpt <- mds$rpt
            ksl <- mds$ksl

            ## update kernel pool
            kpl <- kpl[setdiff(names(kpl), names(ksl))]
        }
        else
            break
    }
    list(par=par, rpt=rpt)
}
