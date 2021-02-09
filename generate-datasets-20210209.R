#####################################################
## Sample generation functions
#####################################################

posdef.matrix <- function(n, ev = runif(n, 0, 10), seed=500){
  ## Function written by Ravi Varadhan: https://stat.ethz.ch/pipermail/r-help/2008-February/153708
  ## Generating a random positive-definite matrix with user-specified positive eigenvalues
  ## If eigenvalues are not specified, they are generated from a uniform distribution
  set.seed(seed)
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

dataset.normal.X2 <- function(s, N, M, beta){
  set.seed(s)
  sd <- 2.0
  r <- 0.8
  Sigma <- sd^2 * r^abs(outer(1:M, 1:M , "-"))
  
  X_1 <- MASS::mvrnorm(N, rep(0,M), posdef.matrix(M, seed=s))
  X_2 <- MASS::mvrnorm(N, rep(0,M), posdef.matrix(M, seed=s*2000))
  epsilon <- MASS::mvrnorm(N, rep(0,M), Sigma, tol=1e-6)
  Y <- matrix(beta[1], N, M) + beta[2]*X_1 + beta[3]*X_2 + epsilon
  data <- as.data.frame(cbind(response=c(t(Y)), X1=c(t(X_1)), X2=c(t(X_2))))
  return(data)
}

dataset.normal.X4 <- function(s, N, M, beta){
  set.seed(s)
  sd <- 2.0
  r <- 0.8
  Sigma <- sd^2 * r^abs(outer(1:M, 1:M , "-"))
  
  X_1 <- MASS::mvrnorm(N, rep(0,M), posdef.matrix(M, seed=s))
  X_2 <- MASS::mvrnorm(N, rep(0,M), posdef.matrix(M, seed=s*2000))
  X_3 <- MASS::mvrnorm(N, rep(0,M), posdef.matrix(M, seed=s*4000))
  X_4 <- MASS::mvrnorm(N, rep(0,M), posdef.matrix(M, seed=s*6000))
  epsilon <- MASS::mvrnorm(N, rep(0,M), Sigma, tol=1e-6)
  Y <- matrix(beta[1], N, M) + beta[2]*X_1 + beta[3]*X_2 + beta[4]*X_3 + beta[5]*X_4 + epsilon
  data <- as.data.frame(cbind(response=c(t(Y)), X1=c(t(X_1)), X2=c(t(X_2)), X3=c(t(X_3)), X4=c(t(X_4))))
  return(data)
}
