work_dir <- "/Users/ehector/Dropbox/Projects/Stream-and-divide/simulations/"
#work_dir <- "/share/research1/ehector/"
#sim <- as.numeric(Sys.getenv('LSB_JOBINDEX'))
simulations <- 1:100

library(Rcpp)
source(paste0(work_dir, "generate-datasets-20210209.R"))
source(paste0(work_dir, "SSLA_func-20210316.R"))
sourceCpp(paste0(work_dir, "increQIF_AR1-20210316.cpp"))

N <- 500
M <- 400
#beta <- c(0.2,-0.2,0.2)
beta <- c(0.2,-0.2,0.2)

B <- 2
K <- 2
corstr <- "AR-1"
family <- "gaussian"

results <- list()
time <- matrix(0, length(simulations), 5)
for(sim in simulations){
  print(paste0("sim=",sim))
  data <- dataset.normal.X2(sim, N, M, beta)
  time[sim,] <- system.time( results[[sim]] <- SSLA(data, N, M, B, K, corstr, family) )
}

z_half_alpha <- qnorm(1-0.05/2) 
P <- length(beta)
results$MSE <- rowMeans(sapply(results[simulations], function(x) (x[,1]-beta)^2))
results$mean_bias <- rowMeans(sapply(results[simulations], function(x) (x[,1]-beta)))
results$beta_mean_var <- rep(0, P)
results$beta_var <- rep(0, P)
results$beta_cov_prob <- rep(0, P)
results$beta_cov_length <- rep(0, P)
results$under_beta <- rep(0, P)
results$type_1_error <- rep(0,P)
for(p in 1:P){
  results$beta_mean_var[p] <- mean(sapply(results[simulations], function(x) x[p,2]^2))
  results$beta_var[p] <- var(sapply(results[simulations], function(x) x[p,1]))
  results$beta_cov_prob[p] <- mean(sapply(results[simulations], function(x) x[p,1]-1.96*x[p,2] < beta[p] & beta[p] < x[p,1]+1.96*x[p,2]))
  results$beta_cov_length[p] <- mean(sapply(results[simulations], function(x) x[p,1]+1.96*x[p,2]-(x[p]-1.96*x[p,2])))
  results$under_beta[p] <- mean(sapply(results[simulations], function(x) x[p,1] - beta[p] < 0))
  results$type_1_error[p] <- mean(sapply(results[simulations], function(x) abs((x[p,1]-beta[p])/x[p,2]) > z_half_alpha))
}
results$mean_time <- colMeans(time)

save.image(paste0(work_dir, "N", N, "M", M, "B", B, "K", K, "X2.RData"))
