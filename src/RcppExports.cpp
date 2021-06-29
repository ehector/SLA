// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// increQIF_ar1
List increQIF_ar1(arma::mat X, arma::vec y, arma::mat x_save, arma::vec y_save, arma::vec nobs, String family, arma::vec beta_old, arma::vec g_accum, arma::mat S_accum, arma::mat C_accum, int maxit, double tol);
RcppExport SEXP _SSLA_increQIF_ar1(SEXP XSEXP, SEXP ySEXP, SEXP x_saveSEXP, SEXP y_saveSEXP, SEXP nobsSEXP, SEXP familySEXP, SEXP beta_oldSEXP, SEXP g_accumSEXP, SEXP S_accumSEXP, SEXP C_accumSEXP, SEXP maxitSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_save(x_saveSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y_save(y_saveSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type nobs(nobsSEXP);
    Rcpp::traits::input_parameter< String >::type family(familySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_old(beta_oldSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type g_accum(g_accumSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S_accum(S_accumSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type C_accum(C_accumSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(increQIF_ar1(X, y, x_save, y_save, nobs, family, beta_old, g_accum, S_accum, C_accum, maxit, tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SSLA_increQIF_ar1", (DL_FUNC) &_SSLA_increQIF_ar1, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_SSLA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
