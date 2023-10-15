library(Rcpp)
library(inline)



rcpp_inc = '
#include<iostream>
using namespace Rcpp;
using namespace arma;


List chfs(vec& countsubs, mat& y, mat& s, mat& r, List& ww);


List chfs(vec& countsubs, mat& y, mat& s, mat& r, List& ww){

  double sumLogD = 0;
  mat YiTemp,SiTemp,RiTemp,ywk,swk,rwk;

  for(int j = 0; j < countsubs.n_elem ; j++){
    
    
    if(j==0){
      //std::cout << j << std::endl;
      int from = 0;
      int to = countsubs(j) - 1;
      YiTemp = y.rows(from, to);
      SiTemp = s.rows(from, to);
      RiTemp = r.rows(from, to);
      mat wwiChol = ww[j];
      wwiChol = chol(wwiChol);
      sumLogD = sumLogD + sum(log(wwiChol.diag()));
  
      
      ywk = solve(trimatl(wwiChol.t()),YiTemp);
      swk = solve(trimatl(wwiChol.t()),SiTemp);
      rwk = solve(trimatl(wwiChol.t()),RiTemp);
      
    }else{
      //std::cout << j << std::endl;
      int from = countsubs(j-1);
      int to = countsubs(j) - 1;
      YiTemp = y.rows(from, to);
      SiTemp = s.rows(from, to);
      RiTemp = r.rows(from, to);
      mat wwiChol = ww[j];
      wwiChol = chol(wwiChol);
      sumLogD = sumLogD + sum(log(wwiChol.diag()));
  
      ywk = join_cols(ywk,solve(trimatl(wwiChol.t()),YiTemp));
      swk = join_cols(swk,solve(trimatl(wwiChol.t()),SiTemp));
      rwk = join_cols(rwk,solve(trimatl(wwiChol.t()),RiTemp));
      
  
    }
  }

  
  return List::create(Named("y.wk") = ywk , Named("s.wk") = swk, Named("r.wk") = rwk, Named("sumLogD") = sumLogD);
}
'


src = '
vec countsubs;
mat y, s, r;
List ww;

countsubs = as<vec>(countsubs_in);
y = as<mat>(y_in);
s = as<mat>(s_in);
r = as<mat>(r_in);
ww = as<List>(ww_in);





return wrap(chfs(countsubs, y, s, r, ww));
'


cholfs.arma = cxxfunction(signature( countsubs_in="numeric",y_in="numeric",s_in="numeric",r_in="numeric",
                                   ww_in="List"), 
                        includes = rcpp_inc, 
                        body = src, 
                        plugin = "RcppArmadillo")


