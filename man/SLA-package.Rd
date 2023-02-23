\name{SLA-package}
\alias{SLA-package}
\alias{SLA}
\docType{package}
\title{\packageTitle{SLA}}
\description{\packageDescription{SLA}}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{SLA}
  \packageIndices{SLA}

  This is a repository for the R package to estimate mean regression parameters for longitudinal data using an online streaming procedure. The R package's main files are:
- src/increQIF_AR1.cpp: this file defines the Rcpp functions that computes the streaming updates for the mean regression parameters.
- R/SLA_func.R: this file defines the R function for the SLA estimation of model parameters.
}
\author{
  \packageAuthor{SLA}
  
  Maintainer: \packageMaintainer{SLA}
}
\references{
# Installation

The SLA R package can be installed in one of two ways:
- from the downloaded gzipped tarball as R CMD INSTALL SLA_1.0-1.tar.gz
- from the downloaded and renamed SLA folder as R CMD build SLA and R CMD INSTALL SLA_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed. If you encounter a library not found error for lgfortran, please try installing gfortran from here: https://cran.r-project.org/bin/macosx/tools/.

# Citation

If you use the SLA R package, please consider citing the relevant manuscript: L. Luo, J. Wang and E.C. Hector (2023). Statistical inference for streamed longitudinal data. Biometrika. doi: 10.1093/biomet/asad010.

# References

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
B <- 50
M <- B*20
beta <- c(0.2,0.5)

data <- dataset.normal.X2.b1sin(s=1, N, M, B, beta)

########### Analysis with SLA ###########

q_seq <- sort( exp(-seq(0.1, 1, length.out = 20)*B^0.3) )
results <- sla(response ~ X1 + X2, data=data, N=N, M=M, B=B, family=family, q=q_seq)

View(results$hist_beta) ## This is the matrix of beta estimates for each batch

View(do.call(rbind, lapply(results$hist_vcov, diag))) ## This is the matrix of beta variance estimates for each batch
}
