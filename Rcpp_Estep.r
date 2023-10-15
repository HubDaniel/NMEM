library(Rcpp)
library(inline)
rcpp_inc = '
#include<iostream>
using namespace Rcpp;
using namespace arma;

List rcpp_loop_w(vec& tiidc, vec& ids, vec& pred_mu1, vec& pred_mu2,
                 vec& p1, vec& p2, vec& yalld_y, vec& yalld_t, List& V1, List& V2);


List rcpp_loop_w(vec& tiidc, vec& ids, vec& pred_mu1, vec& pred_mu2,
                 vec& p1, vec& p2, vec& yalld_y, vec& yalld_t, List& V1, List& V2){
  
  
  vec indxt, predmu1, predmu2, w_1, w_2, exp1, exp2, exp_w1, exp_w2, w_temp1, w_temp2, times, ys, Sigma1_det_log, Sigma2_det_log;
  int n_sub, n_t_obs;
  mat Sigma_inv1, Sigma_inv2;
  
  n_sub         = ids.n_elem;
  w_1           = zeros( n_sub );
  w_2           = zeros( n_sub );

  //std::cout << "a" << std::endl;

  for(int i = 0; i < n_sub; i++){

    if(i==0){
      indxt     = linspace(0,tiidc(i)-1,tiidc(i));
    }else{
      indxt     = linspace(tiidc(i-1),tiidc(i)-1,tiidc(i)-tiidc(i-1));
    }
    
    
    
    // num of observations ,predicted values, observations, time points for the i-th patient
    n_t_obs     = indxt.n_elem;
    int from    = as_scalar(indxt.head(1));
    int to      = as_scalar(indxt.tail(1));
    predmu1     = pred_mu1.subvec(from,to);
    predmu2     = pred_mu2.subvec(from,to);
    times       = yalld_t.subvec(from,to);
    ys          = yalld_y.subvec(from,to);
    
    
    // Sigma_k 
    mat Sigma1      = V1[i];  // list elements have to be declared inline
    mat Sigma2      = V2[i];
    Sigma_inv1      = inv_sympd( Sigma1 );
    Sigma_inv2      = inv_sympd( Sigma2 );
    Sigma1_det_log  = log_det_sympd( Sigma1 );
    Sigma2_det_log  = log_det_sympd( Sigma2 );
    
    
    // ws - ifelse is used to prevent over/underflow
    exp1    = -0.5 * (ys - predmu1).t() * Sigma_inv1 * (ys - predmu1);
    exp2    = -0.5 * (ys - predmu2).t() * Sigma_inv2 * (ys - predmu2);
    
    if(as_scalar(exp2 - exp1) > 650){
      exp_w1  = exp(650 + 0.5 * Sigma1_det_log - 0.5 * Sigma2_det_log);
    }else if(as_scalar(exp2 - exp1) < -650){
      exp_w1  = exp(-650 + 0.5 * Sigma1_det_log - 0.5 * Sigma2_det_log);
    }else{
      exp_w1  = exp(exp2 - exp1 + 0.5 * Sigma1_det_log - 0.5 * Sigma2_det_log);
    }
    
    if(as_scalar(exp1 - exp2) > 650){
      exp_w2  = exp(650 + 0.5 * Sigma2_det_log - 0.5 * Sigma1_det_log);
    }else if(as_scalar(exp1 - exp2) < -650){
      exp_w2  = exp(-650 + 0.5 * Sigma2_det_log - 0.5 * Sigma1_det_log);
    }else{
      exp_w2  = exp(exp1 - exp2 + 0.5 * Sigma2_det_log - 0.5 * Sigma1_det_log);
    }
    
    w_temp1 = 1 / (1 + as_scalar(p2(i))/as_scalar(p1(i)) * exp_w1);
    w_temp2 = 1 / (1 + as_scalar(p1(i))/as_scalar(p2(i)) * exp_w2);
    
    //if(i==10){
    //  std::cout << i << std::endl;
    //  std::cout << w_temp1 << std::endl;
    //  std::cout << exp_w1 << std::endl;
    //  std::cout << as_scalar(p2(i))/as_scalar(p1(i)) << std::endl;
    //}
    
    w_1(i)  = as_scalar(w_temp1);
    w_2(i)  = as_scalar(w_temp2);

  }
  
  return List::create(Named("w_1") = w_1 , Named("w_2") = w_2);
}
'


src = '
vec tiidc, ids, pred_mu1, pred_mu2, p1, p2, yalld_y, yalld_t;
List V1, V2;

tiidc      = as<vec>(tiidc_in);
ids        = as<vec>(ids_in); 
pred_mu1   = as<vec>(pred_mu1_in);
pred_mu2   = as<vec>(pred_mu2_in);
p1         = as<vec>(p1_in);
p2         = as<vec>(p2_in);
yalld_y    = as<vec>(yalld_y_in);
yalld_t    = as<vec>(yalld_t_in);
V1         = as<List>(V1_in);
V2         = as<List>(V2_in);

return wrap(rcpp_loop_w(tiidc, ids, pred_mu1, pred_mu2, p1, p2, yalld_y, yalld_t, V1, V2));
'


rcpp_update_w = cxxfunction(signature(tiidc_in = "numeric", ids_in = "numeric", 
                                      pred_mu1_in = "numeric", pred_mu2_in = "numeric", p1_in = "numeric", 
                                      p2_in = "numeric",yalld_y_in="numeric",yalld_t_in="numeric",
                                      V1_in = "numeric",V2_in = "numeric"), 
                            includes = rcpp_inc, 
                            body = src, 
                            plugin = "RcppArmadillo")





