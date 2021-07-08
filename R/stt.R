## score test

tr <- function(x) sum(diag(x))

es <- function(x) sum(eigen(x, TRUE, TRUE)$values)

qst <- function(Y, K, X, fix, vcs, ...)
{
    m <- X %*% fix
    e <- Y - m
    V <- vsm(K, vcs)
    G <- vsm(K[-1], vcs[-1])
    A <- chol2inv(chol(V))
    h <- m +  G %*% A  %*% e # prediction
    
    r <- Y - h       # residuel, y - mu
    mu2 <- mean(r^2) # mu_2
    mu4 <- mean(r^4) # mu_4

    R <- ksc(K[[2]])
    Q <- t(r) %*% R %*% r / mu2

    EQ <- sum(diag(R))
    VQ <- 2 * sum(diag(R)^2) + (mu4/mu2^2 - 3) * sum(diag(R)^2)
    
    zsc <- drop((Q - EQ) / sqrt(VQ))
    pvl <- 2 * (1 - pnorm(abs(zsc)))

    pvl
    ## list(Z=zsc, P=pvl)
}

## qst <- function(Y, K, X, fix, vcs, ...)
## {
##     N <- NROW(Y)     # sample size
##     muh <- X %*% fix # mu hat
##     sg2 <- vcs[1]    # sigma^2
    
##     ## Q statistics
##     R <- K[[2]]
##     Q <- t(Y - muh) %*% R %*% (Y - muh) / sg2

##     R <- ksc(R)
##     EQ <- sum(diag(R))
##     VQ <- 2 / (N+1) * ((N - 1) * sum(diag(R)^2) - sum(diag(R))^2)
    
##     zsc <- drop((Q - EQ) / sqrt(VQ))
##     pvl <- 2 * (1 - pnorm(abs(zsc)))

##     list(Z=zsc, P=pvl)
## }
