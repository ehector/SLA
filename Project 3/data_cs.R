

sigmagenerator <- function(corst, p, rho = 0){
  if(corst == "ind") { Sigma <- diag(rep(1,p)) }
  if(corst == "cs")  { Sigma <- matrix(rho,p,p); diag(Sigma) <- 1 }
  if(corst == "AR-1") { Sigma <- rho^(abs(outer(1:p, 1:p, "-"))) ; diag(Sigma) <- 1 }
  return(Sigma)
}

datagenerator_cs <- function(n, m, p, B, tempdatadir, type, beta0, intercept=FALSE, corst_x="ind", rho_x=0, 
corst_y="cs", rho_y=0, seed=NULL){
  n          # sample size in each data batch
  m          # number of visits for each individual, a number or a vector of length n
  p          # number of covariates (not including intercept)
  type       # c("gaussian", "binomial", "poisson")
  beta0      # coefficients, length equals to p if intercept=FALSE, and p+1 if intercept=TRUE
  intercept  # if TRUE, the beta0[1] is the coefficient of intercept
  corst_x    # c("ind", "cs", "ar1")
  rho_x      # parameter of correlation
  corst_y    # c("ind", "cs", "ar1")
  rho_y      # parameter of correlation
  seed       # random seed
  
  if(length(n)!=1 & length(n)!=B){ stop("n must be a number or a vector of length B.") }
  if(!is.null(seed)){ set.seed(seed) }
  seed.list <- sample(1:1e8, 1, replace=FALSE)
  if(length(m)!=1 & length(m)!=n){ stop("m must be a positive integer or a vector of length n.") }
  
  dir.create(tempdatadir)
  M <- m * B; # total number of visits for each individual

  set.seed(seed.list)
   
  ix <- matrix(rnorm(n * M), M, n)
  
  id_com <- rep(1:n, each = M)
  visit_com <- rep(1:M, times = n)

  Sigma_y <- sigmagenerator(corst = corst_y, p = M, rho = rho_y)
  ste <- svd(Sigma_y)
  vari_inv <- ste$v %*% diag(sqrt(ste$d)) %*% t(ste$u)
  error <- t(vari_inv %*% ix)
  x <- array(rnorm(n * p * M), c(p, M, n))
  x[1, ,]<-1

  y_m <- matrix(NA, n, M)

    for (j in 1:M){
      y_m[,j]= t(x[, j, ]) %*% beta0 + error[,j]
    } 


  X_com <- matrix(NA, n * M, p)
  y_com <- matrix(NA, n * M, 1)
    for(i in 1:n){
        X_com[(1+(i-1)*M):(M*i),] <- t(x[,,i])
        y_com[(1+(i-1)*M):(M*i)] <- t(y_m[i,])
    }

  X_com <- as.matrix(X_com[,-1])
  
  save(y_com, X_com, id_com, visit_com, file = paste(tempdatadir, "/full", ".RData", sep=""))

   
  for (b in 1:B) {
      id <- rep(1:n, each = m) ## subjects are fixed, but there are more visits in later data batches
      visit <- rep( ((b-1) * m +1) : (b*m), times = n);
      y <- y_com[visit_com == visit]
      X <- X_com[which(visit_com == visit), ]
      save(y, X, id, visit, file = paste(tempdatadir, "/", b, ".RData", sep=""))
  }

}
