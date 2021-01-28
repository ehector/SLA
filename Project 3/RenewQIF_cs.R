
### Sequential updating in QIF
renewqif_cs <- function (B, tempdatadir, family, intercept=FALSE, corstr="exchangeable"){
  
  if (!(family %in% c("gaussian", "binomial", "poisson"))) { stop("'family' value not supported") }
  if (!corstr %in% c("exchangeable")) { stop("'corstr' value not supported") }

  tol <- 1e-4;
  maxit <- 200;
  time_load <- 0;
  Q_stat <- c();
  beta<-c();

  load(paste(tempdatadir,"/",1,".RData",sep=""));
  obs <- lapply(split(id, id), "length")
  nobs <- as.numeric(obs)

  ptm <- proc.time()

  beta_new <- Beta_save <- coef(glm(y ~ X  ,family=family));
  
  if(intercept==TRUE){X<-cbind(1,X)}

  p <- ncol(X)
  n <- length(nobs)
  # for the first data batch, there is no historical data, initialize by 0
  X_save <- matrix(rep(0, n*p), nrow=n);
  Y_save <- matrix(rep(0, n), nrow=n);
  ### initial inference statistics ### 
  N_accum <- 0;

  g_accum <- matrix(rep(0, 2*p),nrow=2*p);
  G_accum <- matrix(rep(0, 2*p*p), nrow=2*p);
  C_accum <- matrix(rep(0, 2*p*2*p), nrow=2*p);
  
  phi_new <- 0;
  
  for(b in 1:B){
    
   	load<-proc.time()
    load(paste(tempdatadir,"/",b,".RData",sep=""))
    time_load <- time_load + (proc.time()- load)[3]

    if(intercept==TRUE){X<-cbind(1,X)}
    obs <- lapply(split(id, id), "length")
    # number of repeated measurements for each subject
    nobs <- as.numeric(obs)

    beta_old <- beta_new;
    beta <- c(beta,beta_old[2]);
    # update beta with the current data batch 
    estimate <- increQIF_cs(X, y, b, X_save, Y_save, Beta_save, nobs, family, beta_old,
      g_accum, G_accum, C_accum, maxit, tol)
      
          beta_new <- estimate$beta;
          g_accum <- estimate$g_accum;
          G_accum <- estimate$G_accum;
          C_accum <- estimate$C_accum;
          
          phi_sub <- estimate$phi;

          N_accum <- N_accum + nrow(X);
          w <- (N_accum-nrow(X)) / N_accum;
          phi_new <- phi_new * w + phi_sub * (1-w);   
      
      x_save_new <- estimate$x_save_new;
      y_save_new <- estimate$y_save_new;
      
      if(b==1){
        X_save <- x_save_new
        Y_save <- y_save_new
        Beta_save <- beta_new
      }else{
        X_save <- rbind(X_save, x_save_new)
        Y_save <- c(Y_save, y_save_new)
        Beta_save <- c(Beta_save, beta_new)
      }
    
  }


  J_accum <- t(G_accum) %*% ginv(C_accum) %*% G_accum;
  varb <- sqrt(diag(ginv(J_accum)));

  out_beta <- cbind(Estimate = drop(beta_new), StdErr = varb, Zscore = beta_new / varb)
  
  time_total <- (proc.time()-ptm)[3]
  time_run <- time_total-time_load 
  out <- list()
  out$beta <- out_beta
  out$phi <- phi_new
  out$time <- time_total
  out$run <- time_run
  out$trace<-beta
  return(out)
}