## use transformation log(exp(eta)+1) in Rcpp for all positive parameters
## varComp4 also use bivariate normal covariance matrix structure instead of cholesky decomposition
## profiled likelihood and different parameter transformation



# outputs:
# profiled likelihood value: Lp
library(Rcpp)
library(inline)


rcpp_inc = '
#include<iostream>
using namespace Rcpp;
using namespace arma;


double lp(vec& eta, vec& tiid, vec& tiidc, vec& Y1, vec& Y2, vec& w1, vec& w2, List& G1, List& G2, List& Z, List& E, List& T);


double lp(vec& eta, vec& tiid, vec& tiidc, vec& Y1, vec& Y2, vec& w1, vec& w2, List& G1, List& G2, List& Z, List& E, List& T){

  double eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9, Lp, t1A, t2A, t1, t2, idsL, ni, w1_i, w2_i, det_V1_log, det_V2_log;
  vec Y1i, Y2i;
  mat V1_i, V2_i;

  t1A = 0;                // initialize the 1st summation result in (29), (30)
  t2A = 0;                // initialize the 2nd summation result in (29), (30)
  eta1 = exp(eta(0));  // sigma2_intercept_1 
  eta2 = 2/(datum::pi) * atan(eta(1));              // rho1
  eta3 = exp(eta(2));  // sigma2_slope_1 
  eta4 = exp(eta(3));  // sigma2_nonparametric_1 
  
  eta5 = exp(eta(4));  // sigma2_intercept_2 
  eta6 = 2/(datum::pi) * atan(eta(5));              // rho2 
  eta7 = exp(eta(6));  // sigma2_slope_2 
  eta8 = exp(eta(7));  // sigma2_nonparametric_2 
  
  eta9 = 1/(1+exp(-eta(8)));  // phi in CAR1, should be between 0 and 1

  idsL = tiidc.n_elem;
  
  

  for(int i = 0; i < idsL; i++){
  
    //std::cout << i << std::endl;
  

    // residual for subject i
    if(i==0){
      int from = 0;
      int to = tiidc(0)-1;
      Y1i = Y1.subvec(from,to);
      Y2i = Y2.subvec(from,to);
    }else{
      int from = tiidc(i-1);
      int to = tiidc(i)-1;
      Y1i = Y1.subvec(from,to);
      Y2i = Y2.subvec(from,to);
    }


    // construct V1_i and V2_i
    ni = tiid(i);
    mat G1_i = G1[i];
    mat G2_i = G2[i];
    mat Z_i = Z[i];
    mat T_i = T[i];

    

    for(int row = 0; row < T_i.n_rows; row++){
      for(int col = 0; col < T_i.n_cols; col++){
        double temppow = pow(eta9,T_i(row,col));
        //if(temppow < 1e-4){
        //  temppow = 0;
        //}
        T_i(row,col) = temppow;
      }
    }    

    mat I(ni,ni,fill::eye);

    // X.submat( first_row, first_col, last_row, last_col )
    G1_i.submat( 0, 0, 0, 0 ) = eta1; // (1,1) the element
    G1_i.submat( 0, 1, 0, 1 ) = eta2 * sqrt(eta1) * sqrt(eta3); // (1,2) the element
    G1_i.submat( 1, 0, 1, 0 ) = eta2 * sqrt(eta1) * sqrt(eta3); // (2,1) the element
    G1_i.submat( 1, 1, 1, 1 ) = eta3; // (2,2) the element
    
    G1_i.cols(2,G1_i.n_cols - 1) = eta4 * G1_i.cols(2,G1_i.n_cols - 1); // lower right matrix need to multiply sigma2_nonparametric_1 
    
    //if(eta9 < 1e-4){
    //  V1_i = Z_i * G1_i * Z_i.t() + I;
    //}else{
    //  V1_i = Z_i * G1_i * Z_i.t() + T_i;
    //}
    
    V1_i = Z_i * G1_i * Z_i.t() + T_i;

    
    
    G2_i.submat( 0, 0, 0, 0 ) = eta5; // (1,1) the element
    G2_i.submat( 0, 1, 0, 1 ) = eta6 * sqrt(eta5) * sqrt(eta7); // (1,2) the element
    G2_i.submat( 1, 0, 1, 0 ) = eta6 * sqrt(eta5) * sqrt(eta7); // (2,1) the element
    G2_i.submat( 1, 1, 1, 1 ) = eta7; // (2,2) the element
    
    G2_i.cols(2,G2_i.n_cols - 1) = eta8 * G2_i.cols(2,G2_i.n_cols - 1); // lower right matrix need to multiply sigma2_nonparametric_2 


    //if(eta9 < 1e-4){
    //  V2_i = Z_i * G2_i * Z_i.t() + I;
    //}else{
    //  V2_i = Z_i * G2_i * Z_i.t() + T_i;
    //}
    
    V2_i = Z_i * G2_i * Z_i.t() + T_i;

    
    if(V1_i.is_sympd()==0){
      vec eigval;
      mat eigvec;
      eig_sym(eigval,eigvec,V1_i);
      if( (abs(eigval(0))<1e-4) & eigval(1)>1e-4 ){
        mat I(ni,ni,fill::eye);
        V1_i = V1_i + 1e-4 * I;
      }
    }
    

    if(V2_i.is_sympd()==0){
      vec eigval;
      mat eigvec;
      eig_sym(eigval,eigvec,V2_i);
      if( (abs(eigval(0))<1e-4) & eigval(1)>1e-4 ){
        mat I(ni,ni,fill::eye);
        V2_i = V2_i + 1e-4 * I;
      }
    } 
    


    // calculate summations
    w1_i = w1(i);
    w2_i = w2(i);
    

    det_V1_log = log_det_sympd( V1_i );
    
    det_V2_log = log_det_sympd( V2_i );

    
    t1 = as_scalar(w1_i * det_V1_log + w2_i * det_V2_log);
    t2 = as_scalar(w1_i * Y1i.t() * inv_sympd( V1_i ) * Y1i + w2_i * Y2i.t() * inv_sympd( V2_i ) * Y2i);

    
    t1A = t1A + t1;
    t2A = t2A + t2;

  }

  
  Lp = tiidc(idsL-1) * log(t2A) + t1A;

  
  return Lp;
}
'


src = '
vec eta, tiid, tiidc, Y1, Y2, w1, w2;
List G1, G2, Z, E, T;

eta = as<vec>(eta_in);
tiid = as<vec>(tiid_in);
tiidc = as<vec>(tiidc_in);
Y1 = as<vec>(Y1_in);
Y2 = as<vec>(Y2_in);
w1 = as<vec>(w1_in);
w2 = as<vec>(w2_in);
G1 = as<List>(G1_in);
G2 = as<List>(G2_in);
Z = as<List>(Z_in);
E = as<List>(E_in);
T = as<List>(T_in);

return wrap(lp(eta, tiid, tiidc, Y1, Y2, w1, w2, G1, G2, Z, E,T));
'


lp4V2_CAR1 = cxxfunction(signature( eta_in="numeric",tiid_in='numeric',tiidc_in='numeric',Y1_in='numeric',Y2_in='numeric',
                               w1_in='numeric',w2_in='numeric', G1_in='numeric',G2_in='numeric', Z_in='numeric', 
                               E_in='numeric',T_in='numeric'), 
                    includes = rcpp_inc, 
                    body = src, 
                    plugin = "RcppArmadillo")






