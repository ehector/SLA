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
  
  X_1 <- matrix(rnorm(N*M), N, M)
  X_2 <- matrix(rnorm(N*M), N, M)
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
  
  X_1 <- matrix(rnorm(N*M), N, M)
  X_2 <- matrix(rnorm(N*M), N, M)
  X_3 <- matrix(rnorm(N*M), N, M)
  X_4 <- matrix(rnorm(N*M), N, M)
  epsilon <- MASS::mvrnorm(N, rep(0,M), Sigma, tol=1e-6)
  Y <- matrix(beta[1], N, M) + beta[2]*X_1 + beta[3]*X_2 + beta[4]*X_3 + beta[5]*X_4 + epsilon
  data <- as.data.frame(cbind(response=c(t(Y)), X1=c(t(X_1)), X2=c(t(X_2)), X3=c(t(X_3)), X4=c(t(X_4))))
  return(data)
}


dataset.binomial.X2 <- function(s, N, M, beta){
  set.seed(s)
  r <- 0.8
  Sigma <- r^abs(outer(1:M, 1:M , "-"))
  
  X_1 <- matrix(rnorm(N*M), N, M)
  X_2 <- matrix(rnorm(N*M), N, M)
  Y <- matrix(0, N, M)
  for(n in 1:N){
    Y[n,] <- SimCorMultRes::rbin(clsize=M, intercepts=beta[1], betas=beta[-1], xformula= ~X_1[n,] + X_2[n,], link="logit", cor.matrix=Sigma)$simdata$y
  }
  data <- as.data.frame(cbind(response=c(t(Y)), X1=c(t(X_1)), X2=c(t(X_2))))
  return(data)
}


dataset.binomial.X4 <- function(s, N, M, beta){
  set.seed(s)
  r <- 0.8
  Sigma <- r^abs(outer(1:M, 1:M , "-"))
  
  X_1 <- matrix(rnorm(N*M), N, M)
  X_2 <- matrix(rnorm(N*M), N, M)
  X_3 <- matrix(rnorm(N*M), N, M)
  X_4 <- matrix(rnorm(N*M), N, M)
  Y <- matrix(0, N, M)
  for(n in 1:N){
    Y[n,] <- SimCorMultRes::rbin(clsize=M, intercepts=beta[1], betas=beta[-1], xformula= ~X_1[n,] + X_2[n,] + X_3[n,] + X_4[n,], link="logit", cor.matrix=Sigma)$simdata$y
  }
  data <- as.data.frame(cbind(response=c(t(Y)), X1=c(t(X_1)), X2=c(t(X_2)), X3=c(t(X_3)), X4=c(t(X_4))))
  return(data)
}
