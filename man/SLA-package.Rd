\name{SLA-package}
\alias{SLA-package}
\alias{SLA}
\docType{package}
\title{\packageTitle{SLA}}
\description{\packageDescription{SLA}}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{SLA}
  \packageIndices{SLA}

  This section provides an example of how to use the package.
}
\author{
  \packageAuthor{SLA}
  
  Maintainer: \packageMaintainer{SLA}
}
\references{
  L. Luo and P. X.-K. Song (2020). Renewable estimation and incremental inference in generalized linear models withstreaming datasets. Journal of the Royal Statistical Society, Series B, 82:69–97.

  A. Qu, B. G. Lindsay and B. Li (2000). Improving generalised estimating equations using quadratic inference functions. Biometrika, 87(4):823–836.

  The posdef.matrix function was written by Ravi Varadhan: https://stat.ethz.ch/pipermail/r-help/2008-February/153708.
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{package}
\seealso{
}
\examples{
########### Generate Data ########### 
N <- 100
B <- 100
M <- B*30
beta <- c(0.2,0.5)

family <- "gaussian"
data <- dataset.normal.X2.b1sin(s=1, N, M, B, beta)

########### Analysis with SLA ###########

q_seq <- sort( exp(-seq(0.1, 1, length.out = 20)*B^0.3) )
results <- sla(response ~ X1 + X2, data=data, N=N, M=M, B=B, family=family, q=q_seq)

View(results$hist_beta) ## This is the matrix of beta estimates for each batch

View(do.call(rbind, lapply(results$hist_vcov, diag))) ## This is the matrix of beta variance estimates for each batch
}
