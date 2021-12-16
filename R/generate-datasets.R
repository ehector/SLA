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

dataset.normal.X4.baseline <- function(s, N, M, beta){
  set.seed(s)
  sd <- 2.0
  r <- 0.8
  Sigma <- sd^2 * r^abs(outer(1:M, 1:M , "-"))
  
  X_1 <- matrix(rnorm(N), N, M)
  X_2 <- matrix(rnorm(N), N, M)
  X_3 <- matrix(rnorm(N), N, M)
  X_4 <- matrix(rnorm(N), N, M)
  epsilon <- MASS::mvrnorm(N, rep(0,M), Sigma, tol=1e-6)
  Y <- matrix(beta[1], N, M) + beta[2]*X_1 + beta[3]*X_2 + beta[4]*X_3 + beta[5]*X_4 + epsilon
  data <- as.data.frame(cbind(response=c(t(Y)), X1=c(t(X_1)), X2=c(t(X_2)), X3=c(t(X_3)), X4=c(t(X_4))))
  return(data)
}

dataset.normal.X4.baseline.b1cos.b2sqrt <- function(s, N, M, B, beta){
  set.seed(s)
  sd <- 2.0
  r <- 0.8
  Sigma <- sd^2 * r^abs(outer(1:M, 1:M , "-"))
  
  X_1 <- matrix(rnorm(N), N, M)
  X_2 <- matrix(rnorm(N), N, M)
  X_3 <- matrix(rnorm(N), N, M)
  X_4 <- matrix(rnorm(N), N, M)
  epsilon <- MASS::mvrnorm(N, rep(0,M), Sigma, tol=1e-6)
  beta1_continuous <- c(sapply(cos((1:B)*2*pi/B), function(x) rep(x, M/B)))
  beta2_continuous <- c(sapply(sqrt(1:B), function(x) rep(x, M/B)))
  Y <- matrix(beta[1], N, M) + t(apply(X_1, 1, function(x) x*beta1_continuous)) + t(apply(X_2, 1, function(x) x*beta2_continuous)) + beta[2]*X_3 + beta[3]*X_4 + epsilon
  data <- as.data.frame(cbind(response=c(t(Y)), X1=c(t(X_1)), X2=c(t(X_2)), X3=c(t(X_3)), X4=c(t(X_4))))
  return(data)
}

dataset.normal.X4.baseline.FFT <- function(s, N, M, beta){
  set.seed(s)
  sd <- 2.0
  r <- 0.8
  A.row <- sd^2 * r^abs(1 - c(1:M, (M-1):2))
  K_p <- length(A.row)
  lambda.sqrt <- sqrt(fft(A.row))
  
  X_1 <- matrix(rnorm(N), N, M)
  X_2 <- matrix(rnorm(N), N, M)
  X_3 <- matrix(rnorm(N), N, M)
  X_4 <- matrix(rnorm(N), N, M)
  epsilon <- t(sapply(1:N, function(x){
    Re(fft(lambda.sqrt * fft(complex(K_p, rnorm(K_p), rnorm(K_p)), inverse = TRUE) / sqrt(K_p))[1:M]) / sqrt(K_p)
  }))
  Y <- matrix(beta[1], N, M) + beta[2]*X_1 + beta[3]*X_2 + beta[4]*X_3 + beta[5]*X_4 + epsilon
  data <- as.data.frame(cbind(response=c(t(Y)), X1=c(t(X_1)), X2=c(t(X_2)), X3=c(t(X_3)), X4=c(t(X_4))))
  return(data)
}

dataset.normal.X4.baseline.b1cos.b2sqrt.FFT <- function(s, N, M, B, beta){
  set.seed(s)
  sd <- 2.0
  r <- 0.8
  A.row <- sd^2 * r^abs(1 - c(1:M, (M-1):2))
  K_p <- length(A.row)
  lambda.sqrt <- sqrt(fft(A.row))
  
  X_1 <- matrix(rnorm(N), N, M)
  X_2 <- matrix(rnorm(N), N, M)
  X_3 <- matrix(rnorm(N), N, M)
  X_4 <- matrix(rnorm(N), N, M)
  epsilon <- t(sapply(1:N, function(x){
    Re(fft(lambda.sqrt * fft(complex(K_p, rnorm(K_p), rnorm(K_p)), inverse = TRUE) / sqrt(K_p))[1:M]) / sqrt(K_p)
  }))
  beta1_continuous <- c(sapply(cos((1:B)*2*pi/B), function(x) rep(x, M/B)))
  beta2_continuous <- c(sapply(sqrt(1:B), function(x) rep(x, M/B)))
  Y <- matrix(beta[1], N, M) + t(apply(X_1, 1, function(x) x*beta1_continuous)) + 
    t(apply(X_2, 1, function(x) x*beta2_continuous)) + beta[2]*X_3 + beta[3]*X_4 + epsilon
  data <- as.data.frame(cbind(response=c(t(Y)), X1=c(t(X_1)), X2=c(t(X_2)), X3=c(t(X_3)), X4=c(t(X_4))))
  return(data)
}

dataset.normal.X4.baseline.b1cos.FFT <- function(s, N, M, B, beta){
  set.seed(s)
  sd <- 2.0
  r <- 0.8
  A.row <- sd^2 * r^abs(1 - c(1:M, (M-1):2))
  K_p <- length(A.row)
  lambda.sqrt <- sqrt(fft(A.row))
  
  X_1 <- matrix(rnorm(N), N, M)
  X_2 <- matrix(rnorm(N), N, M)
  X_3 <- matrix(rnorm(N), N, M)
  X_4 <- matrix(rnorm(N), N, M)
  epsilon <- t(sapply(1:N, function(x){
    Re(fft(lambda.sqrt * fft(complex(K_p, rnorm(K_p), rnorm(K_p)), inverse = TRUE) / sqrt(K_p))[1:M]) / sqrt(K_p)
  }))
  beta1_continuous <- c(sapply(cos((1:B)*2*pi/B), function(x) rep(x, M/B)))
  Y <- matrix(beta[1], N, M) + t(apply(X_1, 1, function(x) x*beta1_continuous)) + beta[2]*X_2 + beta[3]*X_3 + beta[4]*X_4 + epsilon
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
    X_1n <- X_1[n,]
    X_2n <- X_2[n,]
    Y[n,] <- SimCorMultRes::rbin(clsize=M, intercepts=beta[1], betas=beta[-1], xformula= ~X_1n + X_2n, link="logit", cor.matrix=Sigma)$simdata$y
  }
  data <- as.data.frame(cbind(response=c(t(Y)), X1=c(t(X_1)), X2=c(t(X_2))))
  return(data)
}

dataset.binomial.X2.b1sin <- function(s, N, M, B, beta){
  set.seed(s)
  r <- 0.8
  Sigma <- r^abs(outer(1:M, 1:M , "-"))
  
  X_1 <- matrix(rnorm(N*M), N, M)
  X_2 <- matrix(rnorm(N*M), N, M)
  Y <- matrix(0, N, M)
  
  beta_continuous <- c(sapply(4*(1:B)/B*(1-(1:B)/B), function(x) rep(x, M/B)))
  beta_1_M <- rep(beta[1],M)
  beta_2_M <- cbind(beta_continuous, rep(beta[2],M))
  
  for(n in 1:N){
    X_1n <- X_1[n,]
    X_2n <- X_2[n,]
    Y[n,] <- SimCorMultRes::rbin(clsize=M, intercepts=beta_1_M, betas=beta_2_M, xformula= ~X_1n + X_2n, link="logit", cor.matrix=Sigma)$simdata$y
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
    X_1n <- X_1[n,]
    X_2n <- X_2[n,]
    X_3n <- X_3[n,]
    X_4n <- X_4[n,]
    Y[n,] <- SimCorMultRes::rbin(clsize=M, intercepts=beta[1], betas=beta[-1], xformula= ~X_1n + X_2n + X_3n + X_4n, link="logit", cor.matrix=Sigma)$simdata$y
  }
  data <- as.data.frame(cbind(response=c(t(Y)), X1=c(t(X_1)), X2=c(t(X_2)), X3=c(t(X_3)), X4=c(t(X_4))))
  return(data)
}

dataset.binomial.X4.b1cos.b2sqrt <- function(s, N, M, B, beta){
  set.seed(s)
  r <- 0.8
  Sigma <- r^abs(outer(1:M, 1:M , "-"))
  
  X_1 <- matrix(rnorm(N*M), N, M)
  X_2 <- matrix(rnorm(N*M), N, M)
  X_3 <- matrix(rnorm(N*M), N, M)
  X_4 <- matrix(rnorm(N*M), N, M)
  Y <- matrix(0, N, M)
  
  beta1_continuous <- c(sapply(cos((1:B)*2*pi/B), function(x) rep(x, M/B)))
  beta2_continuous <- c(sapply(sqrt(1:B), function(x) rep(x, M/B)))
  beta_1_M <- rep(beta[1],M)
  beta_2_M <- cbind(beta1_continuous, beta2_continuous, rep(beta[2],M), rep(beta[3],M))
  
  for(n in 1:N){
    X_1n <- X_1[n,]
    X_2n <- X_2[n,]
    X_3n <- X_3[n,]
    X_4n <- X_4[n,]
    Y[n,] <- SimCorMultRes::rbin(clsize=M, intercepts=beta_1_M, betas=beta_2_M, xformula= ~X_1n + X_2n + X_3n + X_4n, link="logit", cor.matrix=Sigma)$simdata$y
  }
  data <- as.data.frame(cbind(response=c(t(Y)), X1=c(t(X_1)), X2=c(t(X_2)), X3=c(t(X_3)), X4=c(t(X_4))))
  return(data)
}

dataset.poisson.X2 <- function(s, N, M, beta){
  set.seed(s)
  sd <- 2.0
  r <- 0.8
  A.row <- sd^2 * r^abs(1 - c(1:M, (M-1):2))
  K_p <- length(A.row)
  lambda.sqrt <- sqrt(fft(A.row))
  Sigma <- sd^2 * r^abs(outer(1:M, 1:M , "-"))
  
  X_1 <- matrix(rnorm(N*M), N, M)
  X_2 <- matrix(rnorm(N*M), N, M)
  epsilon <- t(sapply(1:N, function(x){
    Re(fft(lambda.sqrt * fft(complex(K_p, rnorm(K_p), rnorm(K_p)), inverse = TRUE) / sqrt(K_p))[1:M]) / sqrt(K_p)
  }))
  
  X1 <- c(t(X_1))
  X2 <- c(t(X_2))
  Xbeta <- cbind(1,X1,X2)%*%beta
  Y <- qpois(pnorm(c(t(epsilon)), sd = sqrt(diag(Sigma))), exp(Xbeta))
  data <- as.data.frame(cbind(response=Y, X1=X1, X2=X2))
  return(data)
}

dataset.poisson.X4 <- function(s, N, M, beta){
  set.seed(s)
  sd <- 2.0
  r <- 0.8
  A.row <- sd^2 * r^abs(1 - c(1:M, (M-1):2))
  K_p <- length(A.row)
  lambda.sqrt <- sqrt(fft(A.row))
  Sigma <- sd^2 * r^abs(outer(1:M, 1:M , "-"))
  
  X_1 <- matrix(rnorm(N), N, M)
  X_2 <- matrix(rnorm(N), N, M)
  X_3 <- matrix(rnorm(N), N, M)
  X_4 <- matrix(rnorm(N), N, M)
  epsilon <- t(sapply(1:N, function(x){
    Re(fft(lambda.sqrt * fft(complex(K_p, rnorm(K_p), rnorm(K_p)), inverse = TRUE) / sqrt(K_p))[1:M]) / sqrt(K_p)
  }))
  
  X1 <- c(t(X_1))
  X2 <- c(t(X_2))
  X3 <- c(t(X_3))
  X4 <- c(t(X_4))
  Xbeta <- cbind(1,X1,X2,X3,X4)%*%beta
  Y <- qpois(pnorm(c(t(epsilon)), sd = sqrt(diag(Sigma))), exp(Xbeta))
  data <- as.data.frame(cbind(response=Y, X1=X1, X2=X2, X3=X3, X4=X4))
  return(data)
}
