// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;

//[[Rcpp::export]]
List increQIF_ind(arma::mat X, arma::vec y, arma::vec id, String family,
	arma::vec beta_old, arma::vec g_accum, arma::mat G_accum, arma::mat C_accum, 
    int maxit, double tol){
 
    int niter = 0;
    bool stop_flag = FALSE;
    bool converged = FALSE;
    
    int N = X.n_rows;
    int p = X.n_cols;
    arma::vec unique_id = unique(id);
    int n = unique_id.n_elem;

   
    arma::vec g_sum; 
    arma::mat G_temp;
    arma::mat C_temp;

    arma::vec gb_new;
    arma::mat Gb_new;
    arma::mat Cb_new;
   
    //initialization by the beta estimated from previous data
    arma::vec beta_new = beta_old;

    arma::vec mu;
    arma::vec vu;

    arma::mat qif1dev;
    arma::mat qif2dev;

    while (!stop_flag){
    	niter += 1;
        //reset gb to 0 after each iteration

        gb_new = zeros<vec>(p * 2);
        Gb_new = zeros<mat>(p * 2, p);
        Cb_new = zeros<mat>(p * 2, p * 2);
        
        // update gb_new with beta_new over iterations
        arma::vec eta = X * beta_new;

        if(family == "gaussian") {
            mu = eta ;
            vu = ones<vec>(N) ;
        } else if(family == "binomial") {
            mu = exp(eta)/( 1 + exp(eta) ) ;
            vu = exp(eta)/(pow( (1 + exp(eta)), 2 )) ;
        } else if(family == "poisson") {
            mu = exp(eta) ;
            vu = exp(eta) ;
        } else{Rcpp::stop("Unknown distribution family\n");}
    
    
    for (int i = 0; i < n; i++){
        arma::mat Xi = X.rows( find (id == unique_id(i)) ) ;
    	arma::vec yi = y.elem( find (id == unique_id(i)) );
        arma::vec mui = mu.elem( find (id == unique_id(i)) );
        arma::vec vui = vu.elem( find(id == unique_id(i)) ) ;

        int mi = mui.n_elem ;

        arma::mat Ai_half = diagmat(sqrt(vui)) ;
        arma::mat Ai_inv_half = diagmat(1/sqrt(vui)) ;

        mat M1 = zeros<mat>(mi, mi);
        for (int j = 0; j < mi; j++ ){
            for (int k = 0; k < mi; k++){
                if(abs(j-k)==1){ M1(j, k) = 1; }
            }
        }        
       
        vec gi_new = join_cols(Xi.t() * (yi - mui), 
            Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui));
        
        gb_new += gi_new;
        Gb_new += join_cols(Xi.t() * Ai_half * Ai_half * Xi, 
            Xi.t() * Ai_half * M1 * Ai_half * Xi ) ;
        Cb_new += gi_new * gi_new.t();

    }       

        g_sum = g_accum + G_accum * (beta_old - beta_new) + gb_new;

        G_temp = G_accum + Gb_new;
        C_temp = C_accum + Cb_new;

        qif1dev = G_temp.t() * pinv(C_temp) * g_sum ;

        qif2dev = G_temp.t() * pinv(C_temp) * G_temp ;

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
                    Named("G_accum") = G_temp, 
                    Named("C_accum") = C_temp,
                    Named("phi_sub") = phi          
                    );

}