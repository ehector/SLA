// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>
#include <sstream>

using namespace Rcpp; 
using namespace arma;

// [[Rcpp::export]]
arma::mat matrix_inv(const arma::mat& X){
    arma::mat X_inv = inv(X);
    return(X_inv);
}

//[[Rcpp::export]]
double obj_eval_ar1(const arma::mat& X, const arma::vec& y, const arma::vec& nobs, const String& family, const arma::vec& beta){
    
    int N = X.n_rows;
    int p = X.n_cols;
    int n = nobs.n_elem;
    
    arma::vec gb_new = zeros<vec>(p * 2);
    arma::mat g_all_new = zeros<mat>(2 * p, n);
    
    arma::vec mu;
    arma::vec vu;
    
    arma::vec eta = X * beta;
    
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
    
    int loc1 = 0;
    int loc2 = -1;
    for (int i = 0; i < n; i++){
        loc1 = loc2 + 1 ;
        loc2 = loc1 + nobs[i] -1;
        
        arma::mat Xi = X.rows(loc1, loc2);
        arma::vec yi = y.subvec(loc1, loc2);
        arma::vec mui = mu.subvec(loc1, loc2);
        arma::vec vui = vu.subvec(loc1, loc2);
        
        int mi = nobs[i];
        
        arma::mat Ai_half = diagmat(sqrt(vui)) ;
        arma::mat Ai_inv_half = diagmat(1/sqrt(vui)) ;
        
        mat M1 = zeros<mat>(mi, mi);
        for (int j = 0; j < mi; j++ ){
            for (int k = 0; k < mi; k++){
                if(abs(j-k)==1){ M1(j, k) = 1; }
            }
        }        
        
        vec gi_new = join_cols(Xi.t() * (yi - mui), 
                               Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui) );
        gb_new += gi_new;
        g_all_new.col(i) = gi_new;
    } 
    
    double obj_new = as_scalar(gb_new.t() * inv(g_all_new * g_all_new.t()) * gb_new);
    
    return(obj_new);
}

//[[Rcpp::export]]
List increQIF_ar1(const arma::mat& X, const arma::vec& y, const arma::mat& x_save, const arma::vec& y_save, const arma::vec& nobs, 
                  const String& family, const arma::vec& beta_old, const arma::vec& g_accum, const arma::mat& g_all_accum, 
                  const arma::cube& S_i_accum, const arma::mat& S_accum, const double& q,
                  const int& maxit, const double& tol){
    std::ostringstream __my_cerr;
    arma::set_cerr_stream(__my_cerr);
    
    int niter = 0;
    bool stop_flag = FALSE;
    bool converged = FALSE;
    
    int N = X.n_rows;
    int p = X.n_cols;
    int n = nobs.n_elem;

   
    arma::vec g_sum; 
    arma::mat g_all_sum;
    arma::mat S_temp;
    arma::mat C_temp;

    arma::vec gb_new;
    arma::mat g_all_new;
    arma::cube Sb_i_new;
    arma::mat Sb_new;
   
    //initialization by the beta estimated from previous data
    arma::vec beta_new = beta_old;

    arma::vec mu;
    arma::vec vu;
    arma::vec mu_old;
    arma::vec vu_old;
    arma::mat qif1dev;
    arma::mat qif2dev;
    arma::vec residual = zeros<vec>(maxit+1);

    while (!stop_flag){
    	niter += 1;
        //reset gb to 0 after each iteration

        gb_new = zeros<vec>(p * 2);
        g_all_new = zeros<mat>(2 * p, n);
        Sb_i_new = zeros<cube>(n, 2*p, p);
        Sb_new = zeros<mat>(p * 2, p);
        
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
            vec ui_new12 = zeros<vec>(p);
            vec ui_new21 = zeros<vec>(p); 
            mat Si_new = zeros<mat>(p,p);
            bool include = NumericVector::is_na(yi_save);
            if( !include ) {
                ui_new12 = xi_save.t() * (yi(0) - mui(0)) * sqrt(vui_save / vui(0));
                ui_new21 = Xi.row(0).t() * (yi_save - mui_save) * sqrt(vui(0) / vui_save);
                Si_new = xi_save.t() * Xi.row(0) * sqrt (vui(0) * vui_save);
            }
            
            vec gi_new = join_cols(Xi.t() * (yi - mui), 
                                   Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui) + ui_new12 + q*ui_new21 );
            mat S_i_new = join_cols(Xi.t() * Ai_half * Ai_half * Xi, 
                                    Xi.t() * Ai_half * M1 * Ai_half * Xi + Si_new + q*Si_new.t());
            gb_new += gi_new;
            g_all_new.col(i) = gi_new;
            Sb_new += S_i_new ;
            Sb_i_new.row(i) = S_i_new;
            
        } 

        mat S_beta_prod = zeros<mat>(2*p, n);
        for(int r=0; r<2*p; r++){
            mat sub_mat = S_i_accum.col(r);
            S_beta_prod.row(r) = (sub_mat * (beta_old - beta_new)).t();
        }
        g_sum = q*g_accum + q*S_accum * (beta_old - beta_new) + gb_new;
        g_all_sum = q*g_all_accum + q*S_beta_prod + g_all_new;

        S_temp = q*S_accum + Sb_new;
        C_temp = g_all_sum * g_all_sum.t();

        qif1dev = S_temp.t() * inv(C_temp) * g_sum ;

        qif2dev = S_temp.t() * inv(C_temp) * S_temp ;

        vec d_beta = solve(qif2dev, qif1dev);

        double df_beta = as_scalar(qif1dev.t() * d_beta);
        residual(niter-1) = df_beta;

        beta_new += d_beta;

        if(fabs(df_beta) < tol) {converged = TRUE; stop_flag = TRUE;}
        if (niter > maxit) {stop_flag = TRUE;}
    }
    // if (converged==FALSE) {Rcpp::stop("algorithm reached 'maxit' but did not converge\n"); }
    
    double obj_new = as_scalar(gb_new.t() * inv(g_all_new * g_all_new.t()) * gb_new);
    
    std::string text(__my_cerr.str());
    
    return List::create(Named("convergence") = converged,
                        Named("beta") = beta_new,
                        Named("g_accum") = g_sum, 
                        Named("g_all_accum") = g_all_sum,
                        Named("S_i_accum") = S_i_accum + Sb_i_new,
                        Named("S_accum") = S_temp,
                        Named("residual") = residual,
                        Named("Objective") = obj_new,
                        Named("convergence") = text
                    );

}

//[[Rcpp::export]]
List try_increQIF_ar1(const arma::mat& X, const arma::vec& y, const arma::mat& x_save, const arma::vec& y_save, const arma::vec& nobs, 
                  const String& family, const arma::vec& beta_old, const arma::vec& g_accum, const arma::mat& g_all_accum, 
                  const arma::cube& S_i_accum, const arma::mat& S_accum, const double& q,
                  const int& maxit, const double& tol){
    List increQIF;
    bool f = FALSE;
    try{
        increQIF = increQIF_ar1(X, y, x_save, y_save, nobs, family, beta_old, g_accum, g_all_accum, S_i_accum, S_accum, q, maxit, tol);
    } 
    catch (const std::exception &__ex) {
        return(List::create(Named("convergence") = f));
    }
    catch (...) {
        return(List::create(Named("convergence") = f));
    }
    return(increQIF);
}

//[[Rcpp::export]]
List offlineQIF(arma::mat X, arma::vec y, arma::vec nobs, String family, String corstr,
                  arma::vec beta_old, int maxit, double tol){
    
    int niter2= 0;
    bool stop_flag2= FALSE;
    bool converged2= FALSE;
    
    int N = X.n_rows;
    int p = X.n_cols;
    int n = nobs.n_elem;
    
    arma::vec gb_sub;
    arma::mat Gb_sub;
    arma::mat Cb_sub;
    
    //initialization by the beta estimated from previous data
    
    arma::vec beta_sub = beta_old;
    
    arma::vec mu;
    arma::vec vu;
    arma::mat M0;
    arma::mat M1;
    arma::mat M2;
    
    arma::mat qif1dev_sub;
    arma::mat qif2dev_sub;
    
    
    while (!stop_flag2){
        niter2 += 1;
        //reset gb to 0 after each iteration
        if(corstr=="independence"){
            gb_sub = zeros<vec>(p );
            Gb_sub = zeros<mat>(p , p);
            Cb_sub = zeros<mat>(p , p);
        }else if(corstr=="CS+AR1" | corstr=="AR-1(full)"){
            gb_sub = zeros<vec>(p * 3);
            Gb_sub = zeros<mat>(p * 3, p);
            Cb_sub = zeros<mat>(p * 3, p * 3);
        }else{
            gb_sub = zeros<vec>(p * 2);
            Gb_sub = zeros<mat>(p * 2, p);
            Cb_sub = zeros<mat>(p * 2, p * 2);
        }
        
        // update gb_new with beta_new over iterations
        arma::vec eta2 = X * beta_sub;
        
        if(family == "gaussian") {
            mu = eta2 ;
            vu = ones<vec>(N) ;
        } else if(family == "binomial") {
            mu = exp(eta2)/( 1 + exp(eta2) ) ;
            vu = exp(eta2)/(pow( (1 + exp(eta2)), 2 )) ;
        } else if(family == "poisson") {
            mu = exp(eta2) ;
            vu = exp(eta2) ;
        } else{Rcpp::stop("Unknown distribution family\n");}
        
        int loc1 = 0;
        int loc2 = -1;
        for (int i = 0; i < n; i++){
            loc1 = loc2 + 1 ;
            loc2 = loc1 + nobs[i] -1;
            
            arma::mat Xi = X.rows(loc1, loc2) ;
            arma::vec yi = y.subvec(loc1, loc2);
            arma::vec mui = mu.subvec(loc1, loc2);
            arma::vec vui = vu.subvec(loc1, loc2) ;
            
            int mi = nobs[i] ;
            
            arma::mat Ai_half = diagmat(sqrt(vui)) ;
            arma::mat Ai_inv_half = diagmat(1/sqrt(vui)) ;
            
            if (corstr == "independence") {
                M0 = eye<mat>(mi, mi);
                
            } else if (corstr == "exchangeable") {
                M0 = eye<mat>(mi, mi);
                // M1 is a matrix with 1 on off-diagonal elements
                M1 = ones(mi, mi) - M0;  
            } else if (corstr == "AR-1") {
                M0 = eye<mat>(mi, mi);
                M1 = zeros<mat>(mi, mi);
                for (int j = 0; j < mi; j++ ){
                    for (int k = 0; k < mi; k++){
                        if(abs(j-k)==1){ M1(j, k) = 1; }
                    }
                }           
            } else if (corstr == "AR-1(full)"){
                M0 = eye<mat>(mi, mi);
                M1 = zeros<mat>(mi, mi);
                for (int j = 0; j < mi; j++ ){
                    for (int k = 0; k < mi; k++){
                        if(abs(j-k)==1){ M1(j, k) = 1; }
                    }
                }   
                M2 = zeros<mat>(mi, mi);
                M2(0,0) = 1;
                M2(mi-1,mi-1) = 1;
            } else if (corstr == "unstructured"){
                int m = N/n;
                M0 = eye<mat>(m, m);
                M1 = zeros<mat>(m, m);
                
                int loc3 = 0;
                int loc4 = -1;
                for (int i = 0; i < n; i++){
                    loc3 = loc4 +1;
                    loc4 = loc3 + nobs[i] -1;
                    arma::mat Xi = X.rows(loc3, loc4);
                    arma::vec yi = y.subvec(loc3, loc4);
                    arma::vec mui = mu.subvec(loc3, loc4);
                    M1 += (yi - mui) * (yi - mui).t(); 
                }
                M1 = M1 / n;
            } else if(corstr=="CS+AR1"){
                M0 = eye<mat>(mi, mi);
                M1 = ones(mi, mi) - M0;  
                M2 = zeros<mat>(mi, mi);
                for (int j = 0; j < mi; j++ ){
                    for (int k = 0; k < mi; k++){
                        if(abs(j-k)==1){ M2(j, k) = 1; }
                    }
                }     
            }else {Rcpp::stop("Unknown correlation structure\n");}   
            
            if(corstr=="independence"){
                vec gi_sub = Xi.t() * (yi - mui);
                gb_sub += gi_sub;
                Gb_sub += Xi.t() * Ai_half * M0 * Ai_half * Xi ;
                Cb_sub += gi_sub * gi_sub.t();
                
            }else if(corstr == "CS+AR1" | corstr=="AR-1(full)"){
                vec gi_sub = join_cols( join_cols( Xi.t() * (yi - mui), 
                                                   Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui)),
                                                   Xi.t() * Ai_half * M2 * Ai_inv_half * (yi-mui) );
                gb_sub += gi_sub;
                Gb_sub += join_cols( join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
                                               Xi.t() * Ai_half * M1 * Ai_half * Xi),
                                               Xi.t() * Ai_half * M2 * Ai_half * Xi) ;
                Cb_sub += gi_sub * gi_sub.t();
                
            }else{
                arma::vec gi_sub = join_cols(Xi.t() * (yi - mui), Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui) );
                gb_sub += gi_sub;
                Gb_sub += join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
                                    Xi.t() * Ai_half * M1 * Ai_half * Xi) ;
                Cb_sub += gi_sub * gi_sub.t();
            }
            
        }       
        
        qif1dev_sub = Gb_sub.t() * inv(Cb_sub) * gb_sub ;
        
        qif2dev_sub = Gb_sub.t() * inv(Cb_sub) * Gb_sub ;
        
        vec d_beta_sub = solve(qif2dev_sub, qif1dev_sub);
        
        double df_beta_sub = as_scalar(qif1dev_sub.t() * d_beta_sub);
        
        beta_sub += d_beta_sub;
        
        if(fabs(df_beta_sub) < tol) {converged2 = TRUE; stop_flag2 = TRUE;}
        if (niter2 > maxit) {stop_flag2 = TRUE;}
    }
    
    if (converged2==FALSE) {Rcpp::stop("algorithm for single data batch reached 'maxit' but did not converge\n"); }
    
    mat beta_cov = inv(qif2dev_sub);
    
    return List::create(Named("beta") = beta_sub,
                        Named("vcov") = beta_cov,
                        Named("g_sub") = gb_sub,
                        Named("G_sub") = Gb_sub,
                        Named("C_sub") = Cb_sub
    );
}
