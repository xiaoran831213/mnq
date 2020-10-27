#' Cached genotype data
#'
#' A piece of genotype from 1000 genome project.
#'
#' Taken from chromosome 17, region q12. The MAF is no less than 0.05.
if(!exists("c17"))
    c17 <- readRDS("17q12.rds")
if(!exists("f17"))
    f17 <- as.vector(table(as.vector(c17), useNA="no") / sum(!is.na(c17)))
n17 <- nrow(c17)
m17 <- ncol(c17)

#' NL is short for NULL
NL <- NULL

#' retain non-colinear variables
ncv <- function(ldm, lo=0, up=1, ...)
{
    r <- abs(ldm)
    r[upper.tri(r, TRUE)] <- lo + (up - lo) / 2
    b <- which(lo <= r & r <= up, arr.ind=TRUE)
    apply(r, 2, function(.) all(lo <= . & . <= up))
}

#' calculate minor allele frequency
#'
#' @param g genotype matrix, N row samples and M column variants
#' @return MAF of M variants 
maf <- function(g) {a <- colMeans(g, na.rm=TRUE) / 2; pmin(a, 1 - a)}

#' set random element to missing
#'
#' @param x true data
#' @param f frequency of missing
#' @param nan the NaN to use (def=NA).
#' @return the incomplete observation of \code{x}
set.nan <- function(x, f=.1, nan=NA)
{
    x[sample.int(length(x), length(x) * f)] <- nan
    x
}

#' Draw a genotype matrix
#'
#' @param N: samples to draw (rows)
#' @param M: variants to craw (columns)
#' @param psd: postive semi-definite threshold
#' @param bat: the genotype data, if NULL, use beta distribution
#'
#' The between variants correlation must not exceed \code{nlr}.
#'
#' The smallest eigenvalue of the LD matrix must not below (-\code{psd}) * (the
#' largest eigenvalue).
get.gmx <- function(N=5e2, M=5, gno=0, NAF=.1, MAF=.04, psd=NULL, ...)
{
    if(is.null(psd))
        psd <- sqrt(.Machine$double.eps)

    ## M variants is in demand, reserve 3 times more
    P <- min(M * 4, ncol(c17))
    while(TRUE)
    {
        i <- sample.int(n17, N, N > n17)            # N
        j <- seq(sample(m17 - P, 1) + 1, l=P)       # P
        gmx <- c17[i, j]                            # N x P
        gmx <- imp.mod(gmx)                         # no NA
        gmx <- gmx[, maf(gmx) > MAF]                # MAF
        ldm <- cor(gmx)                             # P x P
        
        ## enforce non-linearity; if there are less than M remaining, try again
        ## with a bigger reserve.
        kpp <- ncv(ldm, ...)
        gmx <- gmx[, kpp]
        if(NCOL(gmx) < M)
        {
            P <- min(P + M - NCOL(gmx), m17)
            ## cat("P = ", P, "\n", sep="")
            next
        }

        ## select M variants now
        j <- seq(sample(ncol(gmx) - M, 1) + 1, l=M)
        gmx <- gmx[, j, drop=FALSE]
        ldm <- cor(gmx)
        
        ## if min(eigenvalue) < (threshold) * max(eigenvalue), try again
        egv <- eigen(ldm, TRUE, TRUE)$values
        if (egv[M] < psd * egv[1L])
        {
            cat("Non-PSD!\n")
            next
        }
        break
    }

    ## generate genotype?
    if(gno == 1)
        gmx <- mvrnorm(N, rep(0, M), ldm)

    ## introduce missing
    obs <- set.nan(gmx, NAF)

    ## remove rows that are entirely NA.
    i <- rowSums(is.na(obs)) < M
    gmx <- gmx[i, ]
    obs <- obs[i, ]

    list(gmx=gmx, obs=obs, ldm=ldm)
}


#' Generate phenotype given genotype
#'
#' @param frq: fraction of casual variants;
#' @param hsq: h^2, the desired heritability.
get.phe <- function(gmx, frq=.1, hsq=.1)
{
    N <- nrow(gmx)
    M <- ncol(gmx)

    if(hsq > 0 && frq > 0) # genetic effect exists
    {
        ## the number of casual variants
        Q <- max(1, floor(M * frq))
        
        ## effect sizes (0 for non-casual variants)
        b <- rnorm(M, 0, 1) * sample(c(rep(1, Q), rep(0, M - Q)))
        
        ## genetic effect vector
        gvt <- drop(gmx %*% b)
        
        gvr <- var(gvt)            # variance of genetic effect
        nvr <- gvr / hsq - gvr     # variance of noise
        msk <- b == 0              # mask ineffective variants
    }
    else
    {
        gvt <- 0.0          # no genotype effect
        nvr <- 1.0          # pure noise
        msk <- rep(TRUE, M) # mask ineffective variants
    }
    
    ## noise vector
    nvt <- rnorm(N, 0, sqrt(nvr))
    
    ## genetic effect + noise, standardized to mu=0, var=1
    phe <- gvt + nvt + 1
    ## phe <- as.vector(scale(gvt + nvt))
    list(phe=phe, msk=msk)
}

#' GWAS
#' 
#' @param y: vector of phenotype
#' @param g: genotype matrix
#' @param int: use intercept? (def=No)
get.gwa <- function(y, g, int=1, ...)
{
    ## g <- scale(g, TRUE, FALSE)
    if(int == 1)
    {
        r <- t(apply(g, 2, function(x) summary(lm(y ~ x + 1))$coef[2, ]))
    }
    else if(int == 2)
    {
        r <- drop(cor(y, g) * sqrt(nrow(g)))
        r <- .d(BETA=NA, SE=NA, Z=r, P=pnorm(abs(r), low=FALSE) * 2)
    }
    else
    {
        r <- t(apply(g, 2, function(x) summary(lm(y ~ x - 1))$coef[1, ]))
    }
    colnames(r) <- c("BETA", "SE", "Z", "P")
    rownames(r) <- rownames(g)
    .d(r)
}

#' Try to generate one simulation data
#'
#' @param N sample size
#' @param M number of variants
#' @param hsp intended h^2, the heritability
#' @param MAF minor allele frequency threshold
#' @param ... additional parameters
try.gen <- function(N=500, M=10, hsq=.1, frq=.1, msk=0, ...)
{
    gmx <- get.gmx(N, M, ...)
    phe <- get.phe(gmx$gmx, frq, hsq)

    ## present incomplete data to competing methods
    g <- gmx$obs
    l <- gmx$ldm

    ## present non-effective variants
    if(msk)
    {
        g <- g[, phe$msk, drop=FALSE]
        l <- l[phe$msk, phe$msk, drop=FALSE]
    }

    ## GWAS
    y <- phe$phe
    gwa <- get.gwa(y, g, ...)

    list(gmx=g,       # genotype matrix, with NA(s)
         rsp=y,       # response variable
         tld=l,       # true LD
         zsc=gwa$Z, pvl=gwa$P)
}
