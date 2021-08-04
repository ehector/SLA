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

sla <- function(formula, data, N, M, B, family){
  cl <- match.call()
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  response <- model.response(mf, "numeric")
  covariates <- model.matrix(mt, mf)
  colnames(covariates)[match("X.Intercept.",colnames(covariates))] <- "(Intercept)"
  
  output <- sla_compute(response, covariates, N, M, B, family)
  output <- c(output, list(call=cl, formula=formula))
  names(output$coefficients) <- colnames(covariates)
  colnames(output$vcov) <- rownames(output$vcov) <- colnames(covariates)
  class(output) <- "sla"
  return(output)
}

sla_compute <- function(response, covariates, N, M, B, family){
  m_b <- M/B
  all_inds <- matrix(1:(M*N), M, N)
  
  first_batch_ind <- c(all_inds[1:m_b,])
  beta_new <- as.vector(coef(glm(response[first_batch_ind] ~ 0 + covariates[first_batch_ind,], family=family)))
  
  p <- ncol(covariates)
  
  # for the first data batch, there is no historical data, initialize by NA
  x_save <- matrix(NA, N, p)
  y_save <- rep(NA,N)
  
  ### initial inference statistics ### 
  g_accum <- matrix(0,nrow=2*p, ncol=1)
  g_all_accum <- matrix(0, nrow=2*p, ncol=N)
  S_accum <- matrix(0, nrow=2*p, ncol=p)
  S_i_accum <- array(0, dim=c(N, 2*p, p))
  
  for(b in 1:B){
    ind_b <- c(all_inds[((b-1)*m_b+1):(b*m_b),])
    
    block_y <- as.matrix(response[ind_b])
    block_x <- as.matrix(covariates[ind_b,,drop=FALSE])
    
    beta_old <- beta_new
    
    # update beta with the current data batch 
    estimate <- increQIF_ar1(block_x, block_y, x_save, y_save, nobs=rep(m_b,N), family, beta_old,
                             g_accum, g_all_accum, S_i_accum, S_accum, maxit=10000, tol=1e-6)
    
    beta_new <- estimate$beta
    g_accum <- estimate$g_accum
    g_all_accum <- estimate$g_all_accum
    S_i_accum <- estimate$S_i_accum
    S_accum <- estimate$S_accum
    
    ## save the record for the last visit before loading a new data batch
    x_save <- block_x[seq(m_b,nrow(block_x),m_b),]
    y_save <- block_y[seq(m_b,nrow(block_x),m_b)]
  }
  J_accum <- matrix_inv(t(S_accum) %*% matrix_inv(g_all_accum%*%t(g_all_accum)) %*% S_accum)
  
  out_beta <- list(coefficients = drop(beta_new), vcov = J_accum)
}

offline_qif <- function(response, covariates, N, M, family, corstr){
  beta_old <- as.vector(coef(glm(response ~ 0 + covariates, family=family)))
  
  estimate <- offlineQIF(covariates, response, nobs=rep(M,N), family, corstr=corstr, beta_old, maxit=10000, tol=1e-6)
  
  beta <- estimate$beta
  vcov <- estimate$vcov
  
  out_beta <- list(coefficients = drop(beta), vcov = vcov)
}