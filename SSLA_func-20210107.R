SSLA <- function(data, N, M, B, K, corstr, family){
  m_b <- M/B
  m_kb <- m_b/K
  all_inds <- matrix(1:(M*N), M, N)
  
  for(b in 1:B){
    for(k in 1:K){
      ind_kb <- c(all_inds[((b-1)*m_b+1):(b*m_b),][((k-1)*m_kb+1):(k*m_kb),])
      
      block_y <- as.matrix(data$response[ind_kb])
      block_x <- as.matrix(cbind(intercept=rep(1,m_kb*N), data[ind_kb,-match("response",colnames(data))]))

      init_betas <- coef(glm(block_y ~ 0 + block_x, family=family))
      qif_block_fit <- QIF_sub(block_x, block_y, nobs=rep(m_kb,N), family, corstr, init_betas, maxit=10000, tol=1e-6)
      
      u_kb <- do.call(cbind, qif_block_fit$gi_list)
      V_kb <- u_kb %*% t(u_kb)
      S_kb <- qif_block_fit$G_sub
    }
  }
}