// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;

//[[Rcpp::export]]
List increQIF_cs1(arma::mat X, arma::vec y, int b, arma::mat X_save, arma::vec Y_save, 
    arma::vec nobs, String family,
	arma::vec beta_old, arma::vec g_accum, arma::mat G_accum, arma::mat C_accum, 
    int maxit, double tol){
 
    int niter = 0;
    bool stop_flag = FALSE;
    bool converged = FALSE;
    
    int N = X.n_rows;
    int p = X.n_cols;
    int n = nobs.n_elem;
    
    arma::mat x_save_new = zeros<mat>(n, p);
    arma::vec y_save_new = zeros<vec>(n);
   
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
    arma::vec mu_old;
    arma::vec vu_old;

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
        arma::vec eta_old = X_save * beta_new;
  
        if(family == "gaussian") {
            mu = eta ;
            vu = ones<vec>(N) ;
            mu_old = eta_old ;
            vu_old = ones<vec>(N) ;
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
        
        arma::mat Xi_old = X_save.rows(loc1, loc2);
        arma::vec yi_old = Y_save.subvec(loc1, loc2);
        arma::vec mui_old = mu_old.subvec(loc1, loc2);
        arma::vec vui_old = vu_old.subvec(loc1, loc2);
        
        arma::mat Ai_half_old = diagmat(sqrt(vui_old)) ;
        arma::mat Ai_inv_half_old = diagmat(1/sqrt(vui_old)) ;

        int mi = nobs[i] ;

        arma::mat Ai_half = diagmat(sqrt(vui)) ;
        arma::mat Ai_inv_half = diagmat(1/sqrt(vui)) ;
        
        mat X_temp = Xi.t() * Ai_half * ones(mi, mi);  // p * m, m columns are the same
        vec y_temp = ones(mi, mi) * Ai_inv_half * (yi - mui);  // m * 1, m rows are the same
        
        x_save_new.row(i) = X_temp.col(1).t();  // a row vector
        y_save_new(i) = y_temp(1);

        arma::mat M1 = ones(mi, mi) - eye<mat>(mi, mi);    
    
    int loc3 = 0;
    int loc4 = -1;
    vec ui_new = zeros<vec>(p);   // to aggregate cross batch term for score function
    mat Si_new = zeros<mat>(p, p); // to aggregate cross batch term for sensitivity matrix
 //   for(int k = 0; k < b-1; k++){
    	loc3 = loc4 + 1;
    	loc4 = loc3 + n - 1;
    	arma::mat x_save = X_save.rows(loc3, loc4);
    //	arma::vec y_save = Y_save.subvec(loc3, loc4);  // outcome in data batch k

    	mat xi_save = x_save.row(i) ; // 1 * p covariate for the sum of all visits, this is a row vector !!
    //	double mu_adjust = as_scalar(xi_save * (beta_new - beta_save));
    //    double yi_save = y_save(i); // 1 element
        
        vec ui_newkb = Xi_old.t() * Ai_half_old * ones(mi,mi) * Ai_inv_half * (yi - mui);
        vec ui_newbk = Xi.t() * Ai_half * ones(mi, mi) * Ai_inv_half_old * (yi_old - mui_old);
        ui_new += ui_newkb + ui_newbk;  //  sum over all previous data batches

        Si_new += (x_save_new.row(i)).t() * xi_save + xi_save.t() * x_save_new.row(i);
//    }

    
        vec gi_new = join_cols(Xi.t() * (yi - mui), 
            Xi.t() * Ai_half * M1 * Ai_inv_half * (yi - mui) + ui_new);
        
        gb_new += gi_new;
        Gb_new += join_cols(Xi.t() * Ai_half * Ai_half * Xi, 
            Xi.t() * Ai_half * M1 * Ai_half * Xi + Si_new) ;
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
                    Named("x_save_new") = x_save_new,
                    Named("y_save_new") = y_save_new,
                    Named("phi_sub") = phi          
                    );

}