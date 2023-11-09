# NMEM
The code used for the paper "A Nonparametric Mixed-Effects Mixture Model for Patterns of Clinical Measurements Associated with COVID-19"

https://arxiv.org/abs/2305.04140

"A_main_simulation.R" is the main code for the independent random error case. The results of NMEM in table 3, figure 3 and figure 4 are from here. Please contact Lu and Song to get the code used for Lu and Song's method in table 3, figure 3 and 4. This code contains the estimate of random effects.

"A_main_simulation_CAR1.R" is the main code for the CAR1 random error case. This code is used for checking independence of random errors in section 3.

"util.R" and "util_CAR1.R" are files defining major functions used in the main file.

"gss_newFunctions2.R" is a revision on the original gss package to fit smoothing splines.


All other files are Rcpp files. 

"Rcpp_Estep.r" calculates the Estep of the algorithm.

"Rcpp_cholForwardInSSANOVA.r" calculates the block diagonal matrix cholesky decomposition and is called in "gss_newFunctions2.R".

"Rcpp_generateData.r" is used for generating simulation data.

"Rcpp_varComp4V2.R" and "Rcpp_varComp4V2_CAR1.R" are Cpp version of profiled likelihood function to be optimized.

"Rcpp_varComp4_gradientV2.R" and "Rcpp_varComp4_gradientV2_CAR1.R" are Cpp version of the gradient of profiled likelihood function to be optimized.

"Rcpp_varCompEv3.R" and "Rcpp_varCompEv3_CAR1.R" are used to get the estimate of sigma^2 since we profiled it out.
