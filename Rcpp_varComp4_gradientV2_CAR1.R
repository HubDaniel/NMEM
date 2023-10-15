# This provides the gradient function for the optim
# profiled likelihood and different parameter transformation



# outputs:
# profiled likelihood value: Lp
library(Rcpp)
library(inline)


rcpp_inc = '
#include<iostream>
using namespace Rcpp;
using namespace arma;


vec lp(vec& eta, vec& tiid, vec& tiidc, vec& Y1, vec& Y2, vec& w1, vec& w2, List& G1, List& G2, List& Z, List& E, List& T);


vec lp(vec& eta, vec& tiid, vec& tiidc, vec& Y1, vec& Y2, vec& w1, vec& w2, List& G1, List& G2, List& Z, List& E, List& T){

  double K1, K, eta1, eta2, eta3, eta4, eta5, eta6, eta7, eta8, eta9, idsL, ni, w1_i, w2_i, d1, d2, d3, d4, d5, d6, d7, d8, d9;
  vec Y1i, Y2i, gr;
  mat V1_i, V2_i;
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

  // calculate K first
  
  K1 = 0;
  K = 0; 
  for(int i = 0; i < idsL; i++){
    
    //std::cout << eta << std::endl;
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
    mat invV1_i = inv_sympd(V1_i);
    mat invV2_i = inv_sympd(V2_i);
    K1 = K1 + as_scalar(w1_i * Y1i.t() * invV1_i * Y1i + w2_i * Y2i.t() * invV2_i * Y2i);
    

  }
  
  K = tiidc(idsL-1) / K1;


  d1 = 0;
  d2 = 0;
  d3 = 0;
  d4 = 0;
  d5 = 0;
  d6 = 0;
  d7 = 0;
  d8 = 0;
  d9 = 0;
  
  for(int i = 0; i < idsL; i++){

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
    

    //std::cout << "gr- b" << std::endl;

    w1_i = w1(i);
    w2_i = w2(i);
    mat invV1_i = inv_sympd(V1_i);
    mat invV2_i = inv_sympd(V2_i);
    //X.submat( first_row, first_col, last_row, last_col )
    mat E1 = E[i]; // E matrix for eta(0)
    E1.submat( 0, 0, 0, 0 ) = exp(eta(0));
    E1.submat( 0, 1, 0, 1 ) = 2/(datum::pi) * atan(eta(1)) * sqrt(exp(eta(2))) * exp(eta(0)) / (2*sqrt(exp(eta(0))));
    E1.submat( 1, 0, 1, 0 ) = 2/(datum::pi) * atan(eta(1)) * sqrt(exp(eta(2))) * exp(eta(0)) / (2*sqrt(exp(eta(0))));
    E1.submat( 2, 2, E1.n_rows-1, E1.n_cols-1 ) = E1.submat( 2, 2, E1.n_rows-1, E1.n_cols-1 ) * 0;
    mat ZE1Z = Z_i * E1 * Z_i.t();

    
    mat E2 = E[i]; // E matrix for eta(1)
    E2.submat( 0, 1, 0, 1 ) = 2/(datum::pi) * sqrt(exp(eta(0))*exp(eta(2))) / (1+pow(eta(1),2));
    E2.submat( 1, 0, 1, 0 ) = 2/(datum::pi) * sqrt(exp(eta(0))*exp(eta(2))) / (1+pow(eta(1),2));
    E2.submat( 2, 2, E2.n_rows-1, E2.n_cols-1 ) = E2.submat( 2, 2, E2.n_rows-1, E2.n_cols-1 ) * 0;
    mat ZE2Z = Z_i * E2 * Z_i.t();
    
    mat E3 = E[i]; // E matrix for eta(2)
    E3.submat( 0, 1, 0, 1 ) = 2/(datum::pi) * atan(eta(1)) * sqrt(exp(eta(0))) * exp(eta(2)) / (2*sqrt(exp(eta(2))));
    E3.submat( 1, 0, 1, 0 ) = 2/(datum::pi) * atan(eta(1)) * sqrt(exp(eta(0))) * exp(eta(2)) / (2*sqrt(exp(eta(2))));
    E3.submat( 1, 1, 1, 1 ) = exp(eta(2));
    E3.submat( 2, 2, E3.n_rows-1, E3.n_cols-1 ) = E3.submat( 2, 2, E3.n_rows-1, E3.n_cols-1 ) * 0;
    mat ZE3Z = Z_i * E3 * Z_i.t();
    
    mat E4 = E[i]; // E matrix for expeta4
    E4.submat( 2, 2, E4.n_rows-1, E4.n_cols-1 ) = E4.submat( 2, 2, E4.n_rows-1, E4.n_cols-1 ) * exp(eta(3));
    mat ZE4Z = Z_i * E4 * Z_i.t();
    

    
    mat E5 = E[i]; // E matrix for eta(4)
  
    E5.submat( 0, 0, 0, 0 ) = exp(eta(4));
    E5.submat( 0, 1, 0, 1 ) = 2/(datum::pi) * atan(eta(5)) * sqrt(exp(eta(6))) * exp(eta(4)) / (2*sqrt(exp(eta(4))));
    E5.submat( 1, 0, 1, 0 ) = 2/(datum::pi) * atan(eta(5)) * sqrt(exp(eta(6))) * exp(eta(4)) / (2*sqrt(exp(eta(4))));
    E5.submat( 2, 2, E5.n_rows-1, E5.n_cols-1 ) = E5.submat( 2, 2, E5.n_rows-1, E5.n_cols-1 ) * 0;
    mat ZE5Z = Z_i * E5 * Z_i.t();
    
    mat E6 = E[i]; // E matrix for eta(5)
    E6.submat( 0, 1, 0, 1 ) = 2/(datum::pi) * sqrt(exp(eta(4))*exp(eta(6))) / (1+pow(eta(5),2));
    E6.submat( 1, 0, 1, 0 ) = 2/(datum::pi) * sqrt(exp(eta(4))*exp(eta(6))) / (1+pow(eta(5),2));
    E6.submat( 2, 2, E6.n_rows-1, E6.n_cols-1 ) = E6.submat( 2, 2, E6.n_rows-1, E6.n_cols-1 ) * 0;
    mat ZE6Z = Z_i * E6 * Z_i.t();
    
    mat E7 = E[i]; // E matrix for eta(6)
    E7.submat( 0, 1, 0, 1 ) = 2/(datum::pi) * atan(eta(5)) * sqrt(exp(eta(4))) * exp(eta(6)) / (2*sqrt(exp(eta(6))));
    E7.submat( 1, 0, 1, 0 ) = 2/(datum::pi) * atan(eta(5)) * sqrt(exp(eta(4))) * exp(eta(6)) / (2*sqrt(exp(eta(6))));
    E7.submat( 1, 1, 1, 1 ) = exp(eta(6));
    E7.submat( 2, 2, E7.n_rows-1, E7.n_cols-1 ) = E7.submat( 2, 2, E7.n_rows-1, E7.n_cols-1 ) * 0;
    mat ZE7Z = Z_i * E7 * Z_i.t();
    
    mat E8 = E[i]; // E matrix for expeta8
    E8.submat( 2, 2, E8.n_rows-1, E8.n_cols-1 ) = E8.submat( 2, 2, E8.n_rows-1, E8.n_cols-1 ) * exp(eta(7));
    mat ZE8Z = Z_i * E8 * Z_i.t();
    
    
    mat E9 = T[i];

    
    for(int row = 0; row < E9.n_rows; row++){
      for(int col = 0; col < E9.n_cols; col++){
        if(row!=col){
          E9(row,col) = E9(row,col) * pow(eta9,E9(row,col)-1) * eta9 / (1 + exp(eta(8)));
        }
      }
    }

    // derivatives
    d1 = d1 + as_scalar(-K * (w1_i * Y1i.t() * invV1_i * ZE1Z * invV1_i * Y1i) + w1_i * trace(invV1_i * ZE1Z));
    d2 = d2 + as_scalar(-K * (w1_i * Y1i.t() * invV1_i * ZE2Z * invV1_i * Y1i) + w1_i * trace(invV1_i * ZE2Z));
    d3 = d3 + as_scalar(-K * (w1_i * Y1i.t() * invV1_i * ZE3Z * invV1_i * Y1i) + w1_i * trace(invV1_i * ZE3Z));
    d4 = d4 + as_scalar(-K * (w1_i * Y1i.t() * invV1_i * ZE4Z * invV1_i * Y1i) + w1_i * trace(invV1_i * ZE4Z));
    

    d5 = d5 + as_scalar(-K * (w2_i * Y2i.t() * invV2_i * ZE5Z * invV2_i * Y2i) + w2_i * trace(invV2_i * ZE5Z));
    d6 = d6 + as_scalar(-K * (w2_i * Y2i.t() * invV2_i * ZE6Z * invV2_i * Y2i) + w2_i * trace(invV2_i * ZE6Z));
    d7 = d7 + as_scalar(-K * (w2_i * Y2i.t() * invV2_i * ZE7Z * invV2_i * Y2i) + w2_i * trace(invV2_i * ZE7Z));
    d8 = d8 + as_scalar(-K * (w2_i * Y2i.t() * invV2_i * ZE8Z * invV2_i * Y2i) + w2_i * trace(invV2_i * ZE8Z));
    
    d9 = d9 + as_scalar(-K * (w1_i * Y1i.t() * invV1_i * E9 * invV1_i * Y1i + w2_i * Y2i.t() * invV2_i * E9 * invV2_i * Y2i) + w1_i * trace(invV1_i * E9) + w2_i * trace(invV2_i * E9));
    
  }
  
  

  mat DA1 = zeros( 1,1 ); 
  DA1.submat( 0, 0, 0, 0 ) = d1;
  mat DA2 = zeros( 1,1 ); 
  DA2.submat( 0, 0, 0, 0 ) = d2;
  mat DA3 = zeros( 1,1 ); 
  DA3.submat( 0, 0, 0, 0 ) = d3;
  mat DA4 = zeros( 1,1 ); 
  DA4.submat( 0, 0, 0, 0 ) = d4;
  mat DA5 = zeros( 1,1 ); 
  DA5.submat( 0, 0, 0, 0 ) = d5;
  mat DA6 = zeros( 1,1 ); 
  DA6.submat( 0, 0, 0, 0 ) = d6;
  mat DA7 = zeros( 1,1 ); 
  DA7.submat( 0, 0, 0, 0 ) = d7;
  mat DA8 = zeros( 1,1 ); 
  DA8.submat( 0, 0, 0, 0 ) = d8;  
  mat DA9 = zeros( 1,1 ); 
  DA9.submat( 0, 0, 0, 0 ) = d9;
  

  DA1.insert_cols(DA1.n_cols, DA2);
  DA1.insert_cols(DA1.n_cols, DA3);
  DA1.insert_cols(DA1.n_cols, DA4);
  DA1.insert_cols(DA1.n_cols, DA5);
  DA1.insert_cols(DA1.n_cols, DA6);
  DA1.insert_cols(DA1.n_cols, DA7);
  DA1.insert_cols(DA1.n_cols, DA8);
  DA1.insert_cols(DA1.n_cols, DA9);

  gr = vectorise(DA1);

  return gr;
  
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

return wrap(lp(eta, tiid, tiidc, Y1, Y2, w1, w2, G1, G2, Z, E, T));
'


lp4_grV2_CAR1 = cxxfunction(signature( eta_in="numeric",tiid_in='numeric',tiidc_in='numeric',Y1_in='numeric',Y2_in='numeric',
                                  w1_in='numeric',w2_in='numeric', G1_in='numeric',G2_in='numeric', Z_in='numeric', 
                                  E_in='numeric',T_in='numeric'), 
                       includes = rcpp_inc, 
                       body = src, 
                       plugin = "RcppArmadillo")






