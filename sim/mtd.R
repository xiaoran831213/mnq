#' SKAT
#'
#' @param y: vector of phenotype
#' @param g: genotype matrix
#' @param int: use intercept (def=no)?
#' @return the P-value
skt <- function(y, g, ...)
{
    obj <- SKAT_Null_Model(y ~ 1)
    skt <- SKAT(g, obj, 'linear', is_check_genotype = TRUE, is_dosage = FALSE)
    list(P=skt$p.value)
}

#' multi-regression analysis
#'
#' @param y: vector of phenotype
#' @param g: genotype matrix
#' @param int: use intercept (def=no)?
#' @return the degree of freedom, the  sum of squre, the F-test statistics, and
#'     the P-value of F test.
mra <- function(y, g, int=1, ...)
{
    g <- imp.avg(g)
    if(int > 0)
        ret <- anova(lm(y ~ g + 1))[1, c(1, 3, 4, 5)]
    else
        ret <- anova(lm(y ~ g - 1))[1, c(1, 3, 4, 5)]
    names(ret) <- c('df', 'ssq', 'F', 'P')
    ret
}

#' test of the sum of squares
tsq <- function(Z, C, d=NULL, eps=NULL, ...)
{
    if(is.null(eps))
        eps <- sqrt(.Machine$double.eps)
    S <- sum(Z^2) # sum square

    ## chi-square mixture
    W <- eigen(C, TRUE, TRUE)$value
    L <- length(Z)
    egv <- sum(W > eps)
    P <- imhof(S, W, delta=rep(0, length(Z)))$Qq
    list(Z=Z, C=C, W=W, Y=S, P=P, L=L)
}

#' False Discovery Rate
#'
#' @param P p-values
fdr <- function(P, ...) list(P=min(p.adjust(P, 'fdr')))

#' bonferroni correction
#'
#' @param P p-values
bon <- function(P, ...) list(P=min(p.adjust(P, 'bon')))
