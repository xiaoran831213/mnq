#' impute nothing
imp.nul <- function(g, ...) g

#' Imputation by Mode
#'
#' @param g genotype matrix, N row samples and M column variants
#' @return complete matrix
imp.mod <- function(g, ...)
{
    apply(g, 2L, function(x)
    {
        i <- which(is.na(x))
        u <- unique(x[-i])
        a <- tabulate(match(x[-i], u))
        x[i] <- rep(u[a == max(a)], l=length(i))
        x
    })
}

imp.avg <- function(g, ...)
{
    apply(g, 2L, function(x)
    {
        i <- which(is.na(x))
        x[i] <- mean(x, na.rm=TRUE)
        x
    })
}

imp.rnd <- function(g, ...)
{
    apply(g, 2L, function(x)
    {
        i <- which(is.na(x))
        x[i] <- sample(x[-i], length(i), TRUE)
        x
    })
}


#' Impute missing genotype calls
#'
#' Soft imputation of missing values considering the correlation among variants
#'
#' @details
#' The imputation can be seen as predicting  each of the \code{M} variant by the
#' rest \code{M  - 1} variants,  and taking  the weighted average  prediction of
#' missed calls.
#'
#' @param g genotype matrix in allele dosage format
#' @return the same matrix with imputed dosage values
#' @noRd
imp.box <- function(g)
{
    N <- nrow(g); M <- ncol(g); Z <- is.na(g); C <- !Z

    .zs <- function(g)
    {
        ## column center; setting NA to 0 after opt out incomplete pairs.
        x <- as.matrix(scale(g, TRUE, FALSE))
        x[Z] <- 0

        ## sum_i(x_ir * x_is * i_rs), sum of product between of x_r and x_s,
        ## complete pairs only
        xy <- crossprod(x)

        ## sum_i(x_ir * x_ir * i_rs), sum of square of x_r in the complete
        ## pairs of x_r and x_s
        xx <- crossprod(x^2, C)
        
        ## sum_i(x_ir * x_is * i_rs) / sum_i(x_ir  * x_ir * i_rs) = xy / xx_rs,
        ## the regression coeficient
        rc <- xy / xx

        ## sum_i(x_ir * i_rs) / sum_i(i_rs), mean  of x_r in the complete pairs
        ## of x_r and x_s
        rs <- crossprod(x, C) # the sum
        nn <- crossprod(C)    # the non-NA count
        mu <- rs / nn         # the mean

        ## sum_i(x_is * x_is * i_rs) - 2 * rc sum_i(x_ir * x_is * i_rs) + rc^2 *
        ## sum_i(x_ir * x_ir * i_rs), squared residual
        e2 <- t(xx) - xy * xy / xx

        ## denominator
        d2 <- xx - 2 * rs * mu + mu^2
        
        ## squared standard error
        s2 <- e2 / d2 / (nn - 2)
        diag(s2) <- (N - diag(nn)) / ((diag(nn) - 2) * (diag(nn) - 1))
        ## z-scores
        s2[s2 <= 0 ] <- min(s2[s2 > 0]) # avoid s2=0 caused by perfect fit
        rc / sqrt(s2)
    }
    
    p <- g
    x <- array(c(g == 0 & C, g == 1 & C, g == 2 & C), c(N, M, 3)) * 1
    
    ## g <- imp.mod(g)
    c_0 <- Inf # consistancy
    w <- 1

    ## complete pairs, for each genotype {0, 1, 2}
    n <- array(apply(x, 3, crossprod, C), c(M, M, 3))

    ## weights
    w <- .zs(p)^2
    w <- w / mean(diag(w))

    ## transition
    y <- array(0, c(N, M, 3))
    for(i in 1:3)
    {
        n_i <- crossprod(x[, , i], C) # pairwise complete count for x == 0, 1, or 2
        r_i <- 1/ n_i
        r_i[is.infinite(r_i)] <- 0
        for(j in 1:3)
        {
            y[, , j] <- y[, , j] + x[, , i] %*% (w * crossprod(x[, , i], x[, , j]) * r_i)
        }
    }

    ## balance the contribution of predictors.
    y <- y / rowSums(C)
    y <- y / array(rowSums(y, dims=2), c(N, M, 3))
    
    ## imputation
    p <- 0 * y[, , 1] + 1* y[, , 2] + 2 * y[, , 3]

    g[Z] <- p[Z]
    g
}

#' mean and covariance of incomplete multivariate normal data
#'
#' Assuming the  row samples in an  imcomplete \code{N * M}  matrix \code{x}
#' follows multivate  normal, using expectation conditional  maximization (ECM)
#' to estimate its \code{N * 1} mean vector and \code{M * M} covariance matrix.
#'
#' @param x \code{N * M} matrix with \code{N} samples of a \code{M}-dimensional
#' random vector;
#' @param itr Maximum number of iterations for ECM algorithm (def=50);
#' @param tol Convergence tolerance for ECM (def=1e-8)
#' @param msg print messages on screen?
#'
#' @return the miximum likelihood 
cor.ema <- function(x, itr=NULL, tol=NULL, msg=0)
{
    if(is.null(itr))
        itr <- 50
    if(is.null(tol))
        tol <- sqrt(.Machine$double.eps)

    ## step 2 - initialization
    N <- nrow(x)
    M <- ncol(x)
    Nan <- is.na(x)

    DOF <- N - M - sum(rowSums(Nan) == M)
    stopifnot(DOF > 0)

    NanCols <- colSums(Nan)
    stopifnot(all(NanCols <= N - 2))

    os <- list()
    x0 <- x
    m0 <- colMeans(x0, na.rm=TRUE)
    v0 <- cov(x0, use="pair")
    ## o0 <- obj(x0, m0, v0)
    ## step 3 - main loop
    while(itr > 0)
    {
        ## Step 4 - mean expectation and conditional maximization
        m1 <- rep(0, M)
        c1 <- 0
        for(n in seq(N))
        {
            d <- x0[n, ]
            a <- is.na(d)
            b <- !a
            if(any(a))
            {
                .s <- try(solve(v0[b, b], (d - m0)[b]), TRUE)
                if(inherits(.s, "try-error"))
                    next
                d[a] <- m0[a] + v0[a, b, drop=FALSE] %*% .s
            }
            m1 <- m1 + d
            c1 <- c1 + 1
        }
        m1 <- m1 / c1

        ## step 5 covariance expectation and conditional maximization
        v1 <- matrix(0, M, M)
        c1 <- 0
        for(n in seq(N))
        {
            d <- x0[n, ]
            a <- is.na(d)
            b <- !a
            e <- matrix(0, M, M)
            if(any(a))
            {
                .s <- try(solve(v0[b, b, drop=FALSE], (d - m1)[b]), TRUE)
                .t <- try(solve(v0[b, b, drop=FALSE], v0[b, a]), TRUE)
                if(inherits(.s, "try-error") || inherits(.t, "try-error"))
                    next
                if(inherits(try(v0[a, b] %*% .s), "try-error"))
                {
                    print("aaa")
                }
                d[a] <- m0[a] + v0[a, b] %*% .s
                e[a, a] <- v0[a, a] - v0[a, b] %*% .t
            }
            v1 <- v1 + tcrossprod(d - m1) + e
            c1 <- c1 + 1
        }
        v1 <- v1 / c1
        
        ## Step 6 - evaluate objective and test for convergence
        dif <- mce(v0, v1)
        if(msg) # print out?
            cat(sprintf("%04d: %-7.1e", itr, dif), "\n", sep="")
        if(abs(dif) < tol) # stop?
            break
        itr <- itr - 1 # the next iteration
        m0 <- m1
        v0 <- v1
    }
    cov2cor(v1)
}
