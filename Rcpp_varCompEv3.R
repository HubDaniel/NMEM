library(Rcpp)
library(inline)


rcpp_inc = '
#include<iostream>
using namespace Rcpp;
using namespace arma;


double lp(vec& tiid, vec& tiidc, vec& Y1, vec& Y2, vec& w_1, vec& w_2, List& G1, List& G2, List& Z);


double lp(vec& tiid, vec& tiidc, vec& Y1, vec& Y2, vec& w_1, vec& w_2, List& G1, List& G2, List& Z){

  double Lp, t1A, t1, idsL, ni, w1_i, w2_i, det_V1_log, det_V2_log;
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
    mat I(ni, ni, fill::eye);
    V1_i = Z_i * G1_i * Z_i.t() + I;
    V2_i = Z_i * G2_i * Z_i.t() + I;
    
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
List G1, G2, Z;

tiid = as<vec>(tiid_in);
tiidc = as<vec>(tiidc_in);
Y1 = as<vec>(Y1_in);
Y2 = as<vec>(Y2_in);
w1 = as<vec>(w1_in);
w2 = as<vec>(w2_in);
G1 = as<List>(G1_in);
G2 = as<List>(G2_in);
Z = as<List>(Z_in);

return wrap(lp(tiid, tiidc, Y1, Y2, w1, w2, G1, G2, Z));
'


lpev3 = cxxfunction(signature( tiid_in='numeric',tiidc_in='numeric',Y1_in='numeric',Y2_in='numeric',
                               w1_in='numeric',w2_in='numeric', G1_in='numeric',G2_in='numeric', Z_in='numeric'), 
                   includes = rcpp_inc, 
                   body = src, 
                   plugin = "RcppArmadillo")


