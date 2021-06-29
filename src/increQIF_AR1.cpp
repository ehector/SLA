// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;

//[[Rcpp::export]]
List increQIF_ar1(arma::mat X, arma::vec y, arma::mat x_save, arma::vec y_save, arma::vec nobs, String family, arma::vec beta_old, 
                  arma::vec g_accum, arma::mat S_accum, arma::mat C_accum, int maxit, double tol){
    
    int niter = 0;
    bool stop_flag = FALSE;
    bool converged = FALSE;
    
    int N = X.n_rows;
    int p = X.n_cols;
    int n = nobs.n_elem;

   
    arma::vec g_sum; 
    arma::mat S_temp;
    arma::mat C_temp;

    arma::vec gb_new;
    arma::mat Sb_new;
    arma::mat Cb_new;
   
    //initialization by the beta estimated from previous data
    arma::vec beta_new = beta_old;

    arma::vec mu;
    arma::vec vu;
    arma::vec mu_old;
    arma::vec vu_old;
    arma::mat qif1dev;
    arma::mat qif2dev;

    while (!stop_flag){
    	niter += 1;
        //reset gb to 0 after each iteration

        gb_new = zeros<vec>(p * 2);
        Sb_new = zeros<mat>(p * 2, p);
        Cb_new = zeros<mat>(p * 2, p * 2);
        
        // update gb_new with beta_new over iterations
        arma::vec eta = X * beta_new;
        arma::vec eta_old = x_save * beta_new;

        if(family == "gaussian") {
            mu = eta ;
            vu = ones<vec>(N) ;
            mu_old = eta_old ;
            vu_old = ones<vec>(n) ;
        } else if(family == "binomial") {
            mu = exp(eta)/( 1 + exp(eta) ) ;
            vu = exp(eta)/(pow( (1 + exp(eta)), 2 )) ;
            mu_old = exp(eta_old)/( 1 + exp(eta_old) ) ;
            vu_old = exp(eta_old)/(pow( (1 + exp(eta_old)), 2 )) ;
        } else if(family == "poisson") {
            mu = exp(eta) ;
            vu = exp(eta) ;
            mu_old = exp(eta_old) ;
            vu_old = exp(eta_old) ;
        } else{Rcpp::stop("Unknown distribution family\n");}
    

        int loc1 = 0;
        int loc2 = -1;
        for (int i = 0; i < n; i++){
            loc1 = loc2 + 1 ;
            loc2 = loc1 + nobs[i] -1;
            
            arma::mat Xi = X.rows(loc1, loc2);
            arma::vec yi = y.subvec(loc1, loc2);
            arma::vec mui = mu.subvec(loc1, loc2);
            arma::vec vui = vu.subvec(loc1, loc2);
            
            mat xi_save = x_save.row(i) ; // p * 1 covariate for last visit, this is a row vector!
            double yi_save = y_save(i); // 1 element
            double mui_save = mu_old(i);
            double vui_save = vu_old(i);
            
            int mi = nobs[i];
            
            arma::mat Ai_half = diagmat(sqrt(vui)) ;
            arma::mat Ai_inv_half = diagmat(1/sqrt(vui)) ;
            
            mat M1 = zeros<mat>(mi, mi);
            for (int j = 0; j < mi; j++ ){
                for (int k = 0; k < mi; k++){
                    if(abs(j-k)==1){ M1(j, k) = 1; }
                }
            }        
            
            vec ui_new12 = xi_save.t() * (yi(0) - mui(0)) * sqrt(vui_save / vui(0));
            vec ui_new21 = Xi.row(0).t() * (yi_save - mui_save) * sqrt(vui(0) / vui_save);
            
            mat Si_new = xi_save.t() * Xi.row(0) * sqrt (vui(0) * vui_save);
            
            vec gi_new = join_cols(Xi.t() * (yi - mui), 
                                   Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui) + ui_new12 + ui_new21 );
            
            gb_new += gi_new;
            Sb_new += join_cols(Xi.t() * Ai_half * Ai_half * Xi, 
                                Xi.t() * Ai_half * M1 * Ai_half * Xi + Si_new + Si_new.t()) ;
            Cb_new += gi_new * gi_new.t();
            
        } 

        g_sum = g_accum + S_accum * (beta_old - beta_new) + gb_new;

        S_temp = S_accum + Sb_new;
        C_temp = C_accum + Cb_new;

        qif1dev = S_temp.t() * pinv(C_temp) * g_sum ;

        qif2dev = S_temp.t() * pinv(C_temp) * S_temp ;

        vec d_beta = solve(qif2dev, qif1dev);

        double df_beta = as_scalar(qif1dev.t() * d_beta);

        beta_new += d_beta;

        if(fabs(df_beta) < tol) {converged = TRUE; stop_flag = TRUE;}
        if (niter > maxit) {stop_flag = TRUE;}
    }
    if (converged==FALSE) {Rcpp::stop("algorithm reached 'maxit' but did not converge\n"); }
    
    vec res = (y - mu)/sqrt(vu);
    double phi = as_scalar(res.t() * res / ( N - p ) );
    
    return List::create(Named("beta") = beta_new,
                    Named("g_accum") = g_sum, 
                    Named("S_accum") = S_temp, 
                    Named("C_accum") = C_temp,
                    Named("phi_sub") = phi          
                    );

}
