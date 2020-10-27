## score test

tr <- function(x) sum(diag(x))

qst <- function(Y, K, X, fix, vcs, ...)
{
    mu <- X %*% fix         # mu - vector
    mu2 <- mean((Y - mu)^2) # mu_2
    mu4 <- mean((Y - mu)^4) # mu_4

    ## R <- Reduce("+", mapply("*", K, vcs, SIMPLIFY=FALSE))
    ## R <- Reduce("+", K)
    R <- cov2cor(K[[1]])
    Q <- t(Y - mu) %*% R %*% (Y - mu) / mu2

    EQ <- sum(diag(R))
    ## VQ <- 2 * sum(R * R) + (mu4/mu2^2 - 3) * sum(diag(R)^2)
    VQ <- 2 * sum(diag(R)^2) + (mu4/mu2^2 - 3) * sum(diag(R)^2)
    
    zsc <- drop((Q - EQ) / sqrt(VQ))
    pvl <- 1 - pnorm(zsc)

    list(Z=zsc, P=pvl)
}
