source("util.R")     
source("Rcpp_varComp4V2.R")
source("Rcpp_varComp4_gradientV2.R")
source("Rcpp_varCompEv3.R")
source("Rcpp_generateData.r")
source("Rcpp_Estep.r")
source("Rcpp_cholForwardInSSANOVA.r")
source("gss_newFunctions2.R")

library(xtable)
library(nlme)
library(tidyr)
library(assist)
library(Rcpp)
library(inline)
library(dplyr)
library(gss)
library(Matrix)
library(MASS)
library(glmnet)
library(scales)
library(matrixcalc)

tm.all <- proc.time()

hypers = hyper()

DATA = data_generation(hypers)

EMresults = ALGORITHM(hypers,DATA)

results = RESULTS(hypers,DATA,EMresults,tm.all)

