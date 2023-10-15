# This is Rcpp version of the function for estimating var. comp.
# ONLY random intercept

# inputs: 
# two parameters in a vector: eta (eta1, eta2)


# outputs:
# profiled likelihood value: Lp
library(Rcpp)
library(inline)


rcpp_inc = '
#include<iostream>
using namespace Rcpp;
using namespace arma;


double lp(vec& tiid, vec& tiidc, vec& Y1, vec& Y2, vec& w_1, vec& w_2, List& G1, List& G2, List& Z, List& T, double& Eta9);


double lp(vec& tiid, vec& tiidc, vec& Y1, vec& Y2, vec& w_1, vec& w_2, List& G1, List& G2, List& Z, List& T, double& Eta9){

  double Lp, t1A, t1, idsL, ni, w1_i, w2_i;
  vec Y1i, Y2i;
  mat V1_i, V2_i;


  t1A = 0;                // initialize the 2nd summation result in (29), (30)
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

    // construct V_i
    ni = tiid(i);
    mat G1_i = G1[i];
    mat G2_i = G2[i];
    mat Z_i = Z[i];
    mat T_i = T[i];

    for(int row = 0; row < T_i.n_rows; row++){
      for(int col = 0; col < T_i.n_cols; col++){
        double temppow = pow(Eta9,T_i(row,col));
        //if(temppow < 1e-4){
        //  temppow = 0;
        //}
        T_i(row,col) = temppow;
      }
    }
    

    mat I(ni,ni,fill::eye);
    //if(eta9 < 1e-4){
    //  V1_i = Z_i * G1_i * Z_i.t() + I;
    //}else{
    //  V1_i = Z_i * G1_i * Z_i.t() + T_i;
    //}
    
    V1_i = Z_i * G1_i * Z_i.t() + T_i;
    
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
    w1_i = w_1(i);
    w2_i = w_2(i);

    t1 = as_scalar(w1_i * Y1i.t() * inv_sympd( V1_i ) * Y1i + w2_i * Y2i.t() * inv_sympd( V2_i ) * Y2i) / tiidc(idsL-1);
    
    t1A = t1A + t1;
  }

  Lp = t1A;
  
  return Lp;
}
'


src = '
vec tiid, tiidc, Y1, Y2, w1, w2;
List G1, G2, Z, T;
double Eta9;

tiid = as<vec>(tiid_in);
tiidc = as<vec>(tiidc_in);
Y1 = as<vec>(Y1_in);
Y2 = as<vec>(Y2_in);
w1 = as<vec>(w1_in);
w2 = as<vec>(w2_in);
G1 = as<List>(G1_in);
G2 = as<List>(G2_in);
Z = as<List>(Z_in);
T = as<List>(T_in);
Eta9 = as<double>(Eta9_in);

return wrap(lp(tiid, tiidc, Y1, Y2, w1, w2, G1, G2, Z, T, Eta9));
'


lpev3_CAR1 = cxxfunction(signature( tiid_in='numeric',tiidc_in='numeric',Y1_in='numeric',Y2_in='numeric',
                               w1_in='numeric',w2_in='numeric', G1_in='numeric',G2_in='numeric', 
                               Z_in='numeric',T_in='numeric',Eta9_in='numeric'), 
                    includes = rcpp_inc, 
                    body = src, 
                    plugin = "RcppArmadillo")




# original function
# lp <- function(eta){
#   t1A <- 0
#   t2A <- 0
#   eta1 <- exp(eta[1])
#   eta2 <- exp(eta[2])
#   
#   for(i in 1:length(ids)){
#     if(i==1){
#       Yi1 <- Y1[1:counts[i]]
#       Yi2 <- Y2[1:counts[i]]
#     }else{
#       Yi1 <- Y1[(counts[i-1]+1):counts[i]]
#       Yi2 <- Y2[(counts[i-1]+1):counts[i]]
#     }
#     
#     ni <- count[i]
#     Zi <- rep(1,ni)
#     detVi1 <- eta1 * ni + 1
#     detVi2 <- eta2 * ni + 1
#     invVi1 <- diag(ni) - eta1/(1+eta1*ni) * Zi %*% t(Zi)
#     invVi2 <- diag(ni) - eta2/(1+eta2*ni) * Zi %*% t(Zi)
#     
#     
#     
#     t1 <- w_1[i] * t(Yi1) %*% invVi1 %*% Yi1 + w_2[i] * t(Yi2) %*% invVi2 %*% Yi2
#     t2 <- w_1[i] * log(detVi1) + w_2[i] * log(detVi2)
#     
#     t1A <- t1A + t1
#     t2A <- t2A + t2
#     
#   }
#   Lp <- nrow(yalld) * log(t1A) + t2A
#   return(as.numeric(Lp))
# }






