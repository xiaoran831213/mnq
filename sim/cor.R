cor.avg <- function(x, ...) cor(imp.avg(x))
cor.cxp <- function(x, ...) cor(imp.box(x))
cor.pwc <- function(x, ...) cor(x, use='pair')
cor.std <- function(x, ...)
{
    x <- scale(x, TRUE, FALSE)
    x[is.na(x)] <- 0
    cor(x)
}


#' Correlation of association statistics with missing values
#'
#' this is method No. 2 under idea context that the true correlation is
#' known.
#'
#' @param x N * M data with missing values.
#' @param r M * M correlation that generated x
cor.md2 <- function(x, r, ...)
{
    ## ratio of missing
    m <- colMeans(is.na(x)) 

    ## missing both m_i and m_j.
    mij <- outer(m, m)
    diag(mij) <- m

    ## missing m_i but not m_j
    xij <- outer(m, 1 - m)
    diag(xij) <- 0

    ## missing m_j but not m_i
    xji <- outer(1 - m, m)
    diag(xji) <- 0

    ## not missing m_i and m_j
    aij <- outer(1 - m, 1 - m)
    diag(aij) <- 1 - m
    
    ## sanity check
    stopifnot(abs(mij + xij + xji + aij - 1) < 1e-6)
    
    ## proportion  of one-complete  one missing  among in  data points  entered
    ## association analysis
    p <- (xij + xji) / (1 - mij)
    
    ## shrink the correlation
    r * (1 - p)
}


try1 <- function()
{
    library(MASS)
    L <- 2
    rho <- -0.9
    oo <- 6e6
    S <- (1 - rho) * diag(1, L, L) + rho * matrix(1, L, L)
    mu <- rep(0, L)
    d <- mvrnorm(n=oo, mu=mu, Sigma=S)
    mi <- 0.25 # proportion of missing at SNP i
    mj <- 0.50 # proportion of missing at SNP j

    d1 <- d[, 1]; s1 <- sd(d1); m1 <- mean(d1)
    d2 <- d[, 2]; s2 <- sd(d2); m2 <- mean(d2)

    ## expected proportion of pairs where  one value is missing excluding pairs
    ## where both values are missing
    p <- (mi * (1 - mj) + mj * (1 - mi)) / (1 - mi * mj)

    ## expected fractions of missing values at SNP i and j
    mmi <- round(oo * mi)
    mmj <- round(oo * mj)
    idx <- seq(1, oo)

    ## indices of data where values are missing
    idx.i <- sort(sample(idx, mmi, replace=FALSE))
    idx.j <- sort(sample(idx, mmj, replace=FALSE))
    y <- d
    y[,1][idx.i] <- rnorm(mmi)
    y[,2][idx.j] <- rnorm(mmj)

    ## remove double missing values
    idx.b <- intersect(idx.i, idx.j) # length(idx.b)/oo --> mi*mj
    y <- y[-idx.b,]
    cat("empirical and expected reduced correlation", cor(y)[1,2], (1-p)*rho, "\n")
    rm(y); rm(d)
###
### output
### empirical and expected reduced correlation -0.3854851 -0.3857143
}

try2 <- function()
{
    library(MASS)
    L <- 2
    N <- 6e6
    rho <- -0.9
    S <- (1 - rho) * diag(1, L, L) + rho * matrix(1, L, L)
    d <- mvrnorm(N, rep(0, L), S)

    mi <- 0.25 # proportion of missing at SNP i
    mj <- 0.50 # proportion of missing at SNP j

    d1 <- d[, 1]; s1 <- sd(d1); m1 <- mean(d1)
    d2 <- d[, 2]; s2 <- sd(d2); m2 <- mean(d2)

    ## expected proportion of pairs where  one value is missing excluding pairs
    ## where both values are missing
    p <- (mi * (1 - mj) + mj * (1 - mi)) / (1 - mi * mj)

    ## expected fractions of missing values at SNP i and j
    mmi <- round(N * mi)
    mmj <- round(N * mj)
    idx <- seq(1, N)

    ## indices of data where values are missing
    idx.i <- sort(sample(idx, mmi, replace=FALSE))
    idx.j <- sort(sample(idx, mmj, replace=FALSE))
    y <- d
    y[,1][idx.i] <- rnorm(mmi)
    y[,2][idx.j] <- rnorm(mmj)

    ## remove double missing values
    idx.b <- intersect(idx.i, idx.j) # length(idx.b)/N --> mi*mj
    y <- y[-idx.b,]
    cat("empirical and expected reduced correlation", cor(y)[1,2], (1-p)*rho, "\n")
    rm(y); rm(d)
###
### output
### empirical and expected reduced correlation -0.3854851 -0.3857143
}
