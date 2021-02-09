work_dir <- "/Users/ehector/Dropbox/Projects/Stream-and-divide/simulations/"
#work_dir <- "/share/research1/ehector/"
#sim <- as.numeric(Sys.getenv('LSB_JOBINDEX'))
sim <- 1

library(Rcpp)
source(paste0(work_dir, "generate-datasets-20210107.R"))
source(paste0(work_dir, "SSLA_func-20210209.R"))
sourceCpp(paste0(work_dir, "increQIF_AR1-20210209.cpp"))
set.seed(sim)

N <- 100
M <- 50
beta <- c(0.5, 0.8, -1.2)
data <- dataset.normal.X2(sim, N, M, beta)

B <- 5
K <- 1
corstr <- "AR-1"
family <- "gaussian"

estimates <- list()
estimates[[sim]] <- SSLA(data, N, M, B, K, corstr, family)
