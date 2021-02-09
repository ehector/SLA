SSLA <- function(data, N, M, B, K, corstr, family){
  if(K!=1) stop("'K' value not currently supported")
  m_b <- M/B
  m_kb <- m_b/K
  all_inds <- matrix(1:(M*N), M, N)
  
  beta_new <- coef(glm(response ~ X1+X2, data=data, family=family))
  
  p <- ncol(data)
  
  # for the first data batch, there is no historical data, initialize by 0
  x_save <- matrix(rep(0, N*M*p), nrow=N*M)
  y_save <- matrix(rep(0,N*M), nrow=N*M)
  
  ### initial inference statistics ### 
  g_accum <- matrix(rep(0, 2*p),nrow=2*p);
  S_accum <- matrix(rep(0, 2*p*p), nrow=2*p);
  C_accum <- matrix(rep(0, 2*p*2*p), nrow=2*p);
  
  for(b in 1:B){
    for(k in 1:K){
      ind_kb <- c(all_inds[((b-1)*m_b+1):(b*m_b),][((k-1)*m_kb+1):(k*m_kb),])
      
      block_y <- as.matrix(data$response[ind_kb])
      block_x <- as.matrix(cbind(intercept=rep(1,m_kb*N), data[ind_kb,-match("response",colnames(data))]))
      
      beta_old <- beta_new
      
      # update beta with the current data batch 
      estimate <- increQIF_ar1(block_x, block_y, x_save, y_save, nobs=rep(m_kb,N), family, beta_old,
                               g_accum, S_accum, C_accum, maxit=10000, tol=1e-6)
      
      beta_new <- estimate$beta
      g_accum <- estimate$g_accum
      S_accum <- estimate$S_accum
      C_accum <- estimate$C_accum
      
      ## save the record for the last visit before loading a new data batch
      x_save <- block_x[seq(m_kb,nrow(block_x),m_kb),]
      y_save <- block_y[seq(m_kb,nrow(block_x),m_kb)]
    }
  }
  J_accum <- t(S_accum) %*% solve(C_accum) %*% S_accum;
  varb <- sqrt(diag(solve(J_accum)));
  
  out_beta <- cbind(Estimate = drop(beta_new), StdErr = varb, z.score = beta_new / varb, p.value=2*pnorm(-abs(beta_new / varb)))
}