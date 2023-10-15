library(Rcpp)
library(inline)
rcpp_inc = '
#include<iostream>
using namespace Rcpp;
using namespace arma;


List Rcpp_GData(double& n_k, vec& group_k, vec& tiV, vec& tiV2);


List Rcpp_GData(double& n_k, vec& group_k, vec& tiV, vec& tiV2){

  vec allt, alln, t, n_days_temp;
  double tiL, n_days;
  
  tiL = tiV.size();
  
  
  for(int i = 0; i < n_k; i++){
  
    n_days_temp = shuffle(tiV2);  // equivalent to randomly sample one number of days observed
    n_days = n_days_temp(0);
    t = shuffle(tiV);             // this and next line is equivalent to randomly choose n_days from -60:0
    t = t.subvec(0,n_days-1);
    t = sort(t);
    allt = join_cols(allt,t);
    int sz = alln.size();
    alln.resize(sz+1);
    alln(sz) = n_days;

  }
  
  List L = List::create(Named("allt") = allt , Named("alln") = alln);
  
  return L;
}
'


src = '
vec group_k,tiV,tiV2;
double n_k;

n_k = as<double>(n_k_in);
group_k = as<vec>(group_k_in);
tiV = as<vec>(tiV_in);
tiV2 = as<vec>(tiV2_in);


return wrap(Rcpp_GData(n_k, group_k,  tiV,  tiV2));
'


Rcpp_GD_loop = cxxfunction(signature( n_k_in="numeric",  group_k_in="numeric",  tiV_in="numeric", tiV2_in="numeric"), 
                            includes = rcpp_inc, 
                            body = src, 
                            plugin = "RcppArmadillo")




