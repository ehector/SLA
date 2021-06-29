# About
Super Speedy Longitudinal Analysis (SSLA)

This is a repository for the R package to estimate mean regression parameters for longitudinal data using an online streaming procedure. The R package's main files are:
- src/increQIF_AR1.cpp: this file defines the Rcpp functions that computes the streaming updates for the mean regression parameters.
- R/SSLA_func.R: this file defines the R function for the SSLA estimation of model parameters.

The SSLA man file contains an example for running the regression models from the paper.

Please email ehector@ncsu.edu with any questions or bug-reports.

# Installation

The SSLA R package can be installed in one of two ways:
- from the downloaded gzipped tarball as R CMD INSTALL SSLA_1.0-1.tar.gz
- from the downloaded and renamed SSLA folder as R CMD build SSLA and R CMD INSTALL SSLA_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed.

# Citation

If you use the SSLA R package, please consider citing the relevant manuscript: Hector & Luo.

# References

E. C. Hector and P. X.-K. Song (2021). A distributed and integrated method of moments for high-dimensional correlated data analysis. Journal of the American Statistical Association, 116(534):805–818.

L. Luo and P. X.-K. Song (2020). Renewable estimation and incremental inference in generalized linear models withstreaming datasets. Journal of the Royal Statistical Society, Series B, 82:69–97.

A. Qu, B. G. Lindsay and B. Li (2000). Improving generalised estimating equations using quadratic inference functions. Biometrika, 87(4):823–836.

The posdef.matrix function was written by Ravi Varadhan: https://stat.ethz.ch/pipermail/r-help/2008-February/153708.
