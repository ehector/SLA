sla <- function(x, ...) UseMethod("sla")

print.sla <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nVariance:\n")
  print(x$vcov)
}

summary.sla <- function(object, ...){
  se <- sqrt(diag(object$vcov))
  zval <- coef(object) / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               z.value = zval,
               p.value = 2*pnorm(-abs(zval)))
  res <- list(call=object$call,
              coefficients=TAB)
  class(res) <- "summary.sla"
  return(res) 
}

print.summary.sla <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
}

sla <- function(formula, data, N, M, B, family, q=NULL){
  cl <- match.call()
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if(is.null(q)) q <- seq(0, 1, length.out=100)
  response <- model.response(mf, "numeric")
  covariates <- model.matrix(mt, mf)
  colnames(covariates)[match("X.Intercept.",colnames(covariates))] <- "(Intercept)"
  
  output <- sla_compute(response, covariates, N, M, B, family, q)
  output <- c(output, list(call=cl, formula=formula))
  names(output$coefficients) <- colnames(covariates)
  colnames(output$vcov) <- rownames(output$vcov) <- colnames(covariates)
  class(output) <- "sla"
  return(output)
}

sla_compute <- function(response, covariates, N, M, B, family, q){
  p <- ncol(covariates)
  m_b <- M/B
  all_inds <- matrix(1:(M*N), M, N)
  Q <- length(q)
  
  first_batch_ind <- c(all_inds[1:m_b,])
  beta_new <- suppressWarnings(as.vector(coef(glm(response[first_batch_ind] ~ 0 + covariates[first_batch_ind,], family=family))))
  beta_save <- matrix(beta_new, nrow=p, ncol=Q)
  
  # for the first data batch, there is no historical data, initialize by NA
  x_save <- matrix(NA, N, p)
  y_save <- rep(NA,N)
  # fitted.values <- rep(NA, length(response))
  
  ### initial inference statistics ### 
  g_accum <- matrix(0,nrow=2*p, ncol=1)
  g_all_accum <- matrix(0, nrow=2*p, ncol=N)
  S_accum <- matrix(0, nrow=2*p, ncol=p)
  S_i_accum <- array(0, dim=c(N, 2*p, p))
  
  # save the beta and optimal q's estimates from each batch
  hist_beta <- matrix(0,B,p)
  hist_vcov <- list()
  q_opt <- rep(NA,B)
  
  for(b in 1:B){
    ind_b <- c(all_inds[((b-1)*m_b+1):(b*m_b),])
    
    block_y <- as.matrix(response[ind_b])
    block_x <- as.matrix(covariates[ind_b,,drop=FALSE])
    
    beta_old <- beta_new
    beta_init <- suppressWarnings(c(coef(glm(block_y ~ 0 + block_x, family=family))))
    
    # update beta with the current data batch 
    estimates <- list()
    for(j in 1:Q){
      estimates[[j]] <- tryCatch(try_increQIF_ar1(block_x, block_y, x_save, y_save, nobs=rep(m_b,N), family, beta_old, beta_init, 
                                     g_accum, g_all_accum, S_i_accum, S_accum, q=q[j], maxit=10000, tol=1e-6),
                                 error = function(c) list(convergence=FALSE)
      )
      if(!estimates[[j]]$convergence & j!=Q) {
        df <- data.frame(convergence=FALSE)
        estimates <- c(estimates, lapply(1:(Q-j), function(x) df))
        break
      }
    }
    converged <- sapply(estimates, function(est) est$convergence)
    QIF_opt <- unlist(sapply(estimates, function(x) x$Objective))
    optimum_candidates <- which(QIF_opt == min(QIF_opt, na.rm=T))
    optimum <- optimum_candidates[1] 
    q_opt[b] <- q[optimum]

    beta_new <- estimates[[optimum]]$beta
    g_accum <- estimates[[optimum]]$g_accum
    g_all_accum <- estimates[[optimum]]$g_all_accum
    S_accum <- estimates[[optimum]]$S_accum
    S_i_accum <- estimates[[optimum]]$S_i_accum
    
    hist_beta[b,] <- beta_new
    hist_vcov[[b]] <- solve(t(estimates[[optimum]]$S_accum)%*%matrix_inv(estimates[[optimum]]$g_all_accum%*%t(estimates[[optimum]]$g_all_accum))%*%estimates[[optimum]]$S_accum)
    
    ## save the record for the last visit before loading a new data batch
    x_save <- block_x[seq(m_b,nrow(block_x),m_b),]
    y_save <- block_y[seq(m_b,nrow(block_x),m_b)]
    # fitted.values[ind_b] <- block_x%*%beta_new
  }
  J_accum <- matrix_inv(t(S_accum) %*% matrix_inv(g_all_accum%*%t(g_all_accum)) %*% S_accum)
  
  out_beta <- list(coefficients = drop(beta_new), vcov = J_accum, #fitted.values=fitted.values, 
                   hist_beta=hist_beta, hist_vcov=hist_vcov, q=q, q_opt=q_opt)
}

online_score <- function(response, covariates, N, M, B, family, q){
  p <- ncol(covariates)
  m_b <- M/B
  all_inds <- matrix(1:(M*N), M, N)
  Q <- length(q)
  
  first_batch_ind <- c(all_inds[1:m_b,])
  beta_new <- suppressWarnings(as.vector(coef(glm(response[first_batch_ind] ~ 0 + covariates[first_batch_ind,], family=family))))
  beta_save <- matrix(beta_new, nrow=p, ncol=Q)
  
  # for the first data batch, there is no historical data, initialize by NA
  x_save <- matrix(NA, N, p)
  y_save <- rep(NA,N)
  # fitted.values <- rep(NA, length(response))
  
  ### initial inference statistics ### 
  g_accum <- matrix(0,nrow=p, ncol=1)
  g_all_accum <- matrix(0, nrow=p, ncol=N)
  S_accum <- matrix(0, nrow=p, ncol=p)
  S_i_accum <- array(0, dim=c(N, p, p))
  
  # save the beta and optimal q's estimates from each batch
  hist_beta <- matrix(0,B,p)
  hist_vcov <- list()
  q_opt <- rep(NA,B)
  
  for(b in 1:B){
    ind_b <- c(all_inds[((b-1)*m_b+1):(b*m_b),])
    
    block_y <- as.matrix(response[ind_b])
    block_x <- as.matrix(covariates[ind_b,,drop=FALSE])
    
    beta_old <- beta_new
    beta_init <- suppressWarnings(c(coef(glm(block_y ~ 0 + block_x, family=family))))
    
    # update beta with the current data batch 
    estimates <- list()
    for(j in 1:Q){
      estimates[[j]] <- tryCatch(online_weighted_score(block_x, block_y, nobs=rep(m_b,N), family, beta_old, beta_init, 
                                                  g_accum, g_all_accum, S_i_accum, S_accum, q=q[j], maxit=10000, tol=1e-6),
                                 error = function(c) list(convergence=FALSE)
      )
      if(!estimates[[j]]$convergence & j!=Q) {
        df <- data.frame(convergence=FALSE)
        estimates <- c(estimates, lapply(1:(Q-j), function(x) df))
        break
      }
    }
    converged <- sapply(estimates, function(est) est$convergence)
    QIF_opt <- unlist(sapply(estimates, function(x) x$Objective))
    optimum_candidates <- which(QIF_opt == min(QIF_opt, na.rm=T))
    optimum <- optimum_candidates[1] 
    q_opt[b] <- q[optimum]
    
    beta_new <- estimates[[optimum]]$beta
    g_accum <- estimates[[optimum]]$g_accum
    g_all_accum <- estimates[[optimum]]$g_all_accum
    S_accum <- estimates[[optimum]]$S_accum
    S_i_accum <- estimates[[optimum]]$S_i_accum
    
    hist_beta[b,] <- beta_new
    hist_vcov[[b]] <- solve(t(estimates[[optimum]]$S_accum)%*%matrix_inv(estimates[[optimum]]$g_all_accum%*%t(estimates[[optimum]]$g_all_accum))%*%estimates[[optimum]]$S_accum)
  }
  J_accum <- matrix_inv(t(S_accum) %*% matrix_inv(g_all_accum%*%t(g_all_accum)) %*% S_accum)
  
  out_beta <- list(coefficients = drop(beta_new), vcov = J_accum, #fitted.values=fitted.values, 
                   hist_beta=hist_beta, hist_vcov=hist_vcov, q=q, q_opt=q_opt)
}
