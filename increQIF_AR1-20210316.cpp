// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp; 
using namespace arma;

//[[Rcpp::export]]
List QIF_sub(arma::mat X, arma::vec y, int K, arma::vec nobs, String family, String corstr, arma::vec beta){
    
    int N = X.n_rows;
    int p = X.n_cols;
    int n = nobs.n_elem;
    
    arma::vec gb_sub;
    arma::mat Gb_sub;
    arma::mat gi_list;
    
    if(corstr=="independence"){
        gi_list = zeros<mat>(n, p*K);
    } else if(corstr=="CS+AR1"){
        gi_list = zeros<mat>(n, p*3*K);
    } else{
        gi_list = zeros<mat>(n, p*2*K);
    }
   
    arma::vec mu;
    arma::vec vu;
    arma::mat M0;
    arma::mat M1;
    arma::mat M2;
    
    if(corstr=="independence"){
        gb_sub = zeros<vec>(p*K);
        Gb_sub = zeros<mat>(p*K, p);
    }else if(corstr=="CS+AR1"){
        gb_sub = zeros<vec>(p*3*K);
        Gb_sub = zeros<mat>(p*3*K, p);
    }else{
        gb_sub = zeros<vec>(p*2*K);
        Gb_sub = zeros<mat>(p*2*K, p);
    }

    arma::vec eta2 = X * beta;
    
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
        int mi = nobs[i] ;
        for(int k = 0; k < K; k++){
            loc1 = loc2 + 1 ;
            loc2 = loc1 + nobs[i] -1;
            
            arma::mat Xi = X.rows(loc1, loc2) ;
            arma::vec yi = y.subvec(loc1, loc2);
            arma::vec mui = mu.subvec(loc1, loc2);
            arma::vec vui = vu.subvec(loc1, loc2) ;
            
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
            }else if(corstr=="CS+AR1"){
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
                arma::vec gi_sub = zeros<vec>(p*K);
                gi_sub.rows(k*p, (k+1)*p-1) = Xi.t() * (yi - mui);
                Gb_sub.rows(k*p, (k+1)*p-1) += Xi.t() * Ai_half * M0 * Ai_half * Xi ;
                gb_sub += gi_sub;
                gi_list.row(i) += gi_sub.t();
                
            }else if(corstr == "CS+AR1"){
                arma::vec gi_sub = zeros<vec>(p*3*K);
                gi_sub.rows(k*p*3, (k+1)*p*3-1) = join_cols( join_cols( Xi.t() * (yi - mui), 
                                                   Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui)),
                                                   Xi.t() * Ai_half * M2 * Ai_inv_half * (yi-mui) );
                Gb_sub.rows(k*p*3, (k+1)*p*3-1) += join_cols( join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
                                               Xi.t() * Ai_half * M1 * Ai_half * Xi),
                                               Xi.t() * Ai_half * M2 * Ai_half * Xi) ;
                gb_sub += gi_sub;
                gi_list.row(i) += gi_sub.t();
                
            }else{
                arma::vec gi_sub = zeros<vec>(p*2*K);
                gi_sub.rows(k*p*2, (k+1)*p*2-1) = join_cols(Xi.t() * (yi - mui), Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui) );
                Gb_sub.rows(k*p*2, (k+1)*p*2-1) += join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
                                    Xi.t() * Ai_half * M1 * Ai_half * Xi) ;
                gb_sub += gi_sub;
                gi_list.row(i) += gi_sub.t();
            }   
        }
    }  
    
    arma::mat Cb_sub = gi_list.t() * gi_list;
    
    arma::mat A = Gb_sub.t() * pinv(Cb_sub) * Gb_sub;
    //arma::mat S = Gb_sub.t() * pinv(Cb_sub) * Gb_sub;
    //arma::vec ones = zeros<vec>(n);
    //return List::create(Named("psi") = psi,
    //                    Named("S") = S,
    //                    Named("V") = Cb_sub,
    //                    Named("A") = A
    //);
    return List::create(Named("S") = Gb_sub,
                        Named("V") = Cb_sub,
                        Named("A") = A
    );
}

//[[Rcpp::export]]
List increQIF_ar1(arma::mat X, arma::vec y, int K, List psi_save, List D_save, List A_save, 
                  arma::vec nobs, String family, String corstr, arma::vec beta_old, 
                  arma::vec G_accum, arma::mat T_accum, arma::mat W_accum, int maxit, double tol){
    
    int niter = 0;
    bool stop_flag = FALSE;
    bool converged = FALSE;
    
    int p = X.n_cols;
    int N = X.n_rows;
    int n = nobs.n_elem;

    arma::mat onep = ones<mat>(p,p);
    arma::vec onep_vec = ones<vec>(p);
    arma::vec onen = ones<vec>(nobs.n_elem);
    
    arma::vec u_ib;
    arma::mat Sb_sub;
    //arma::mat Vb_sub;
    arma::vec G_b;
    arma::mat T_b;
    arma::mat W_b;
    List psi_save_new(n);
    List D_save_new(n);
    List A_save_new(n);
    arma::mat T_temp;
    arma::mat W_temp;
    arma::vec G_temp; 
   
    //initialization by the beta estimated from previous data
    arma::vec beta_new = beta_old;

    arma::vec mu;
    arma::vec vu;
    arma::mat M0;
    arma::mat M1;
    arma::mat M2;
    
    arma::mat qif1dev;
    arma::mat qif2dev;

    while (!stop_flag){
    	niter += 1;
        
        //reset gb to 0 after each iteration
        if(corstr=="independence"){
            u_ib = zeros<vec>(p*K);
            Sb_sub = zeros<mat>(p*K, K);
            //Vb_sub = zeros<mat>(p*K, p*K);
        }else if(corstr=="CS+AR1"){
            u_ib = zeros<vec>(p*3*K);
            Sb_sub = zeros<mat>(p*3*K, p);
            //Vb_sub = zeros<mat>(p*3*K, p * 3*K);
        }else{
            u_ib = zeros<vec>(p*2*K);
            Sb_sub = zeros<mat>(p*2*K, p);
            //Vb_sub = zeros<mat>(p*2*K, p*2*K);
        }
        G_b = zeros<vec>(p*2);
        T_b = zeros<mat>(p*2,p);
        W_b = zeros<mat>(p*2,p*2);
        
        // update gb_new with beta_new over iterations
        arma::vec eta2 = X * beta_new;
        
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
        
        List new_update = QIF_sub(X, y, K, nobs, family, "AR-1", beta_new);
        
        arma::mat Sb_new = new_update[0];
        arma::mat Vb_new = new_update[1];
        arma::mat A_ib = new_update[2];
        
        int loc1 = 0;
        int loc2 = -1;
        for (int i = 0; i < n; i++){
            int mi = nobs[i];
            arma::vec psi_save_i = psi_save[i];
            arma::mat A_save_i = A_save[i];
            arma::mat D_save_i = D_save[i];
            
            arma::mat inv_chol_A_save;
            
            if(as_scalar(onep_vec.t()*A_save_i*onep_vec)==0){
                inv_chol_A_save = zeros<mat>(p,p);
            } else{
                inv_chol_A_save = chol(inv(A_save_i));
            }
            for(int k = 0; k < K; k++){
                loc1 = loc2 + 1 ;
                loc2 = loc1 + mi -1;
                
                arma::mat Xi = X.rows(loc1, loc2) ;
                arma::vec yi = y.subvec(loc1, loc2);
                arma::vec mui = mu.subvec(loc1, loc2);
                arma::vec vui = vu.subvec(loc1, loc2) ;
                
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
                }else if(corstr=="CS+AR1"){
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
                    u_ib.rows(k*p, (k+1)*p-1) = Xi.t() * (yi - mui);
                    Sb_sub.rows(k*p, (k+1)*p-1) = Xi.t() * Ai_half * M0 * Ai_half * Xi ;
                }else if(corstr == "CS+AR1"){
                    u_ib.rows(k*p*3, (k+1)*p*3-1) = join_cols( join_cols( Xi.t() * (yi - mui), 
                                                 Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui)),
                                                 Xi.t() * Ai_half * M2 * Ai_inv_half * (yi-mui) );
                    Sb_sub.rows(k*p*3, (k+1)*p*3-1) = join_cols( join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
                                                  Xi.t() * Ai_half * M1 * Ai_half * Xi),
                                                  Xi.t() * Ai_half * M2 * Ai_half * Xi) ;
                }else{
                    u_ib.rows(k*p*2, (k+1)*p*2-1) = join_cols(Xi.t() * (yi - mui), Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui) );
                    Sb_sub.rows(k*p*2, (k+1)*p*2-1) = join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
                                       Xi.t() * Ai_half * M1 * Ai_half * Xi) ;
                }   
            }
            //Vb_sub = u_ib * u_ib.t();

            arma::vec psib_sub = Sb_new.t() * inv(Vb_new) * u_ib;
            
            //arma::mat A_ib = Sb_new.t() * inv(Vb_new) * Vb_sub * inv(Vb_new)*Sb_new;
            arma::mat D_ib = Sb_new.t()*inv(Vb_new)*Sb_sub;
            arma::vec g_ib_1 = D_ib.t() * inv(A_ib) * psib_sub;
            arma::vec g_ib_2 = D_ib.t()*chol(inv(A_ib))*onep*inv_chol_A_save*psi_save_i + 
                D_save_i.t()*inv_chol_A_save*onep*chol(inv(A_ib))*psib_sub;
            arma::vec g_ib = join_cols(g_ib_1, g_ib_2);
            G_b += g_ib;
            arma::mat T_b_1 = D_ib.t() * inv(A_ib) * D_ib;
            arma::mat T_b_2 = D_ib.t() * chol(inv(A_ib))*onep*inv_chol_A_save*D_save_i + 
                D_save_i.t()*inv_chol_A_save*onep*chol(inv(A_ib))*D_ib;
            T_b += join_cols(T_b_1, T_b_2);
            W_b += g_ib*g_ib.t();
            
            psi_save_new[i] = psib_sub;
            D_save_new[i] = D_ib;
            A_save_new[i] = A_ib;
        } 
        G_temp = G_accum + T_accum * (beta_old - beta_new) + G_b;
        T_temp = T_accum + T_b;
        W_temp = W_accum + W_b;

        qif1dev = T_temp.t() * pinv(W_temp) * G_temp ;

        qif2dev = T_temp.t() * pinv(W_temp) * T_temp ;

        vec d_beta = solve(qif2dev, qif1dev);

        double df_beta = as_scalar(qif1dev.t() * d_beta);

        beta_new += d_beta;

        if(fabs(df_beta) < tol) {converged = TRUE; stop_flag = TRUE;}
        if (niter > maxit) {stop_flag = TRUE;}
    }
    if (converged==FALSE) {Rcpp::stop("algorithm reached 'maxit' but did not converge\n"); }
    
    return List::create(Named("beta") = beta_new,
                    Named("G_accum") = G_temp, 
                    Named("T_accum") = T_temp, 
                    Named("W_accum") = W_temp,
                    Named("psi_save") = psi_save_new,
                    Named("A_save") = A_save_new,
                    Named("D_save") = D_save_new
                    );

}