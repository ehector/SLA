#####################################################
## Sample generation functions
#####################################################

dataset.normal.X2.b1sin <- function(s, N, M, B, beta){
  set.seed(s)
  sd <- 2.0
  r <- 0.8
  Sigma <- sd^2 * r^abs(outer(1:M, 1:M , "-"))
  
  X_1 <- matrix(rnorm(N*M), N, M)
  X_2 <- matrix(rnorm(N*M), N, M)
  epsilon <- MASS::mvrnorm(N, rep(0,M), Sigma, tol=1e-6)
  beta_continuous <- c(sapply(sin((1:B)*2*pi/B), function(x) rep(x, M/B)))
  Y <- matrix(beta[1], N, M) + t(apply(X_1, 1, function(x) x*beta_continuous)) + beta[2]*X_2 + epsilon
  data <- as.data.frame(cbind(response=c(t(Y)), X1=c(t(X_1)), X2=c(t(X_2))))
  return(data)
}
