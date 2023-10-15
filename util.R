hyper <- function(){
  
  sigma2_1_trueo    = sigma2_2_trueo = 0.4136
  sigma2_a_1_trueo  =  0.0955 
  sigma2_a_2_trueo  =  0.0620
  sigma2_as_1_trueo = -0.0809
  sigma2_as_2_trueo = -0.0135
  sigma2_s_1_trueo  =  0.5130
  sigma2_s_2_trueo  =  0.0956
  sigma2_np_1_trueo =  14.8065
  sigma2_np_2_trueo =  14#0.2195
  
  beta_vp0o =  -1.0028   #intercept
  beta_vp1o =   0        #genderM
  beta_vp2o =  -0.3734   #raceW
  beta_vp3o =  0         #ethniHL
  beta_vp4o =  0         #vintage
  beta_vp5o =  0         #diabeteY
  beta_vp6o =  0         #hyperY
  beta_vp7o =   0        #BMI
  beta_vp8o =  -0.0107   #age
  beta_vp9o =   0        #typeCVCATH2
  beta_vp10o =  0        #typeCVCATH4
  
  
  use_estimatedo <- T  # generate data using estimated mean functions
  
  no <- 1000#3293            # total number of subjects
  Tio <- -30:0         # observation window
  
  max_itero <- 100                 # max iteration for EM algorithm
  max_rel_changeo <- 1e-5           # convergence criteria for EM iteration
  
  
  max_iter_innero <- 5              # max iteration for inner algorithm
  max_rel_change_innero <- 1e-5     # convergence criteria for inner iteration
  
  
  
  return(list(sigma2_1_trueo=sigma2_1_trueo,sigma2_2_trueo=sigma2_2_trueo,
              sigma2_a_1_trueo=sigma2_a_1_trueo,sigma2_a_2_trueo=sigma2_a_2_trueo,
              sigma2_as_1_trueo=sigma2_as_1_trueo,sigma2_as_2_trueo=sigma2_as_2_trueo,
              sigma2_s_1_trueo=sigma2_s_1_trueo,sigma2_s_2_trueo=sigma2_s_2_trueo,
              sigma2_np_1_trueo=sigma2_np_1_trueo,sigma2_np_2_trueo=sigma2_np_2_trueo,
              beta_vp0o=beta_vp0o,beta_vp1o=beta_vp1o,beta_vp2o=beta_vp2o,beta_vp3o=beta_vp3o,
              beta_vp4o=beta_vp4o,beta_vp5o=beta_vp5o,beta_vp6o=beta_vp6o,beta_vp7o=beta_vp7o,
              beta_vp8o=beta_vp8o,beta_vp9o=beta_vp9o,beta_vp10o=beta_vp10o,
              use_estimatedo=use_estimatedo,no=no,Tio=Tio,max_itero=max_itero,
              max_rel_changeo=max_rel_changeo,max_iter_innero=max_iter_innero,
              max_rel_change_innero=max_rel_change_innero))
}


# data generation
data_generation <- function( hypers,
                             mu1_true=
                               c( -0.1054,-0.0966,-0.0777,-0.0559,-0.0425,-0.0729,-0.1023,
                                  -0.1126,-0.1192,-0.0915,-0.0603,-0.0341,-0.0091, 0.0241,
                                  0.0058,-0.0303,-0.0823,-0.1289,-0.1260,-0.1153,-0.0939,
                                  -0.0512,-0.0218, 0.0258, 0.0629, 0.1034, 0.1406, 0.2325,
                                  0.4808, 1.0039, 1.7526
                               ),
                             mu2_true=
                               c(  0.0107, 0.0113, 0.0114, 0.0105, 0.0088, 0.0065, 0.0038,
                                   0.0006,-0.0032,-0.0071,-0.0105,-0.0127,-0.0131,-0.0115,
                                   -0.0083,-0.0035, 0.0026, 0.0096, 0.0174, 0.0258, 0.0351,
                                   0.0451, 0.0559, 0.0671, 0.0787, 0.0908, 0.1029, 0.1147,
                                   0.1260, 0.1370, 0.1479
                               )){
  
  
  sigma2_1_true = hypers$sigma2_1_trueo
  sigma2_2_true = hypers$sigma2_2_trueo
  sigma2_a_1_true = hypers$sigma2_a_1_trueo
  sigma2_a_2_true = hypers$sigma2_a_2_trueo
  sigma2_as_1_true = hypers$sigma2_as_1_trueo
  sigma2_as_2_true = hypers$sigma2_as_2_trueo
  sigma2_s_1_true = hypers$sigma2_s_1_trueo
  sigma2_s_2_true = hypers$sigma2_s_2_trueo
  sigma2_np_1_true = hypers$sigma2_np_1_trueo
  sigma2_np_2_true = hypers$sigma2_np_2_trueo
  
  
  beta_vp0 <- hypers$beta_vp0o
  beta_vp1 <- hypers$beta_vp1o
  beta_vp2 <- hypers$beta_vp2o
  beta_vp3 <- hypers$beta_vp3o
  beta_vp4 <- hypers$beta_vp4o
  beta_vp5 <- hypers$beta_vp5o
  beta_vp6 <- hypers$beta_vp6o       
  beta_vp7 <- hypers$beta_vp7o        
  beta_vp8 <- hypers$beta_vp8o   
  beta_vp9 <- hypers$beta_vp9o
  beta_vp10 <- hypers$beta_vp10o
  
  
  use_estimated = hypers$use_estimatedo
  
  
  ### ti
  ti <- hypers$Tio
  tiscaled <- (ti-min(ti))/(max(ti)-min(ti))
  
  ### smooth function
  
  if(hypers$use_estimated){
    # copy and paste estimated mean functions here
    mu1 <- function(timein) {  # timein is a vector of time points, from -60:0, return corresponding values
      extractINDX <- match(timein,ti)
      MEANS <- mu1_true
      
      return(MEANS[extractINDX])
    }
    
    mu2 <- function(timein) {
      extractINDX <- match(timein,ti)
      MEANS <- mu2_true
      
      return(MEANS[extractINDX])
    }
    
  }else{
    mu1 <- function(x){
      97+1*exp(0.25*(x+4))
    }
    mu2 <- function(x) rep(97,length(x))
 
    
  }
  
  
  ### random error (2 groups could have different random error)
  e1 <- function(x) rnorm(n = x,mean = 0,sd = sqrt(sigma2_1_true))
  e2 <- function(x) rnorm(n = x,mean = 0,sd = sqrt(sigma2_2_true))
  
  ### num of patients
  n <- hypers$no
  
  ############################################################################################################################
  ###################################################### Data Generation #####################################################
  ############################################################################################################################
  
  
  ### generate exogenous variables - all these are from the real data distributions, fitdist function is used for exponential distribution
  vp1 <- sample(c(0,1),size = n,replace = T,prob = c(0.44,0.56))  # genderM = 1
  vp2 <- sample(c(0,1),size = n,replace = T,prob = c(0.36,0.64))  # raceW = 1
  vp3 <- sample(c(0,1),size = n,replace = T,prob = c(0.79,0.21))  # ethniHL = 1
  vp4 <- rexp(n,rate = 0.22)                                      # vintage
  vp5 <- sample(c(0,1),size = n,replace = T,prob = c(0.48,0.52))  # diabeteY = 1
  vp6 <- sample(c(0,1),size = n,replace = T,prob = c(0.21,0.79))  # hyperY = 1
  vp7 <- rgamma(n,shape = 15.58,rate = 0.51)                        # BMI
  vp8 <- rnorm(n,mean = 61.69,sd = 14.20)                               # age
  vp910 <- sample(c(1,2,4),size = n,replace = T,prob = c(0.67,0.13,0.20)) # typeCVCATH2 or typeCVCATH4
  vp9 <- ifelse(vp910==2,1,0) # typeCVCATH2 = 1
  vp10 <- ifelse(vp910==4,1,0) # typeCVCATH4 = 1
  
  
  
  
  prob <- exp(beta_vp0+beta_vp1*vp1+beta_vp2*vp2+beta_vp3*vp3+beta_vp4*vp4+beta_vp5*vp5+
                beta_vp6*vp6+beta_vp7*vp7+beta_vp8*vp8+beta_vp9*vp9+beta_vp10*vp10)/
    (1+exp(beta_vp0+beta_vp1*vp1+beta_vp2*vp2+beta_vp3*vp3+beta_vp4*vp4+beta_vp5*vp5+
             beta_vp6*vp6+beta_vp7*vp7+beta_vp8*vp8+beta_vp9*vp9+beta_vp10*vp10))
  

  # soft cluster (randomly assign each patient according to above probability)
  groups <- 1*(runif(n)<prob)
  group1 <- which(groups==1)
  group2 <- which(groups==0)
  n1 <- length(group1)
  n2 <- length(group2)
  
  # all data
  D <- data.frame(id=1:n, prob=prob, group=1, 
                  vp1=vp1,vp2=vp2,vp3=vp3,vp4=vp4,vp5=vp5,vp6=vp6,vp7=vp7,vp8=vp8,vp9=vp9,vp10=vp10)
  D$group[group2] <- 2
  
  
  ###############
  ### group 1 ###
  ###############
  RGDtemp <- Rcpp_GD_loop(n1, group1, ti, (6:23))
  allt <- as.numeric(RGDtemp$allt)
  alln <- as.numeric(RGDtemp$alln)
  
  ### add nonparametric random effect - stochastic process with 0 mean and RK1 evaluated at time points as covariance matrix
  nonparametric_re_all <- c()
  gen_time_re_all <- c()
  
  allncumsum <- cumsum(alln)
  for(i in 1:length(allncumsum)){
    if(i==1){
      timepoints <- allt[1:allncumsum[i]]
    }else{
      timepoints <- allt[(allncumsum[i-1]+1) : allncumsum[i]]
    }
    
    # scale timepoints to be [0,1]
    scaledtime <- (timepoints - ti[1])/(ti[length(ti)]-ti[1])
    # RK1 evaluated at timepoints
    covMat <- sigma2_np_1_true * cubic(scaledtime)
    # generate nonparametric random effects
    gen_np_re <- mvrnorm(n=1,mu=rep(0,length(scaledtime)),Sigma=covMat)
    nonparametric_re_all <- c(nonparametric_re_all,gen_np_re)
    # generate random intercept and random slope (w.r.t time)
    covMat <- matrix(c(sigma2_a_1_true,sigma2_as_1_true,sigma2_as_1_true,sigma2_s_1_true),ncol=2)
    bb <- mvrnorm(n=1,mu=rep(0,2),Sigma=covMat)
    dMat <- cbind(1,scaledtime)
    gen_time_re <- dMat %*% bb
    gen_time_re_all <- c(gen_time_re_all,gen_time_re)
  }
  
  alltscaled <- (allt-ti[1])/(ti[length(ti)]-ti[1])
  mu <- mu1(allt)
  e  <- e1(sum(alln))
  y  <- mu + e + gen_time_re_all + nonparametric_re_all
  y_k1 <- data.frame(y=y,t=allt,id=rep(group1,alln))
  
  
  ###############
  ### group 2 ###
  ###############
  
  RGDtemp <- Rcpp_GD_loop(n2, group2, ti, (6:23))
  allt <- as.numeric(RGDtemp$allt)
  alln <- as.numeric(RGDtemp$alln)
  
  ### add nonparametric random effect - stochastic process with 0 mean and RK1 evaluated at time points as covariance matrix
  nonparametric_re_all <- c()
  gen_time_re_all <- c()
  
  allncumsum <- cumsum(alln)
  for(i in 1:length(allncumsum)){
    if(i==1){
      timepoints <- allt[1:allncumsum[i]]
    }else{
      timepoints <- allt[(allncumsum[i-1]+1) : allncumsum[i]]
    }
    # scale timepoints to be [0,1] -- question: for each patient, using their max and min to scale OR use common max=0 and min=-60 to scale?
    scaledtime <- (timepoints - ti[1])/(ti[length(ti)]-ti[1])
    # RK1 evaluated at timepoints
    covMat <- sigma2_np_2_true * cubic(scaledtime)
    # generate random effect realizations
    gen_np_re <- mvrnorm(n=1,mu=rep(0,length(scaledtime)),Sigma=covMat)
    nonparametric_re_all <- c(nonparametric_re_all,gen_np_re)
    # generate random intercept and random slope (w.r.t time)
    covMat <- matrix(c(sigma2_a_2_true,sigma2_as_2_true,sigma2_as_2_true,sigma2_s_2_true),ncol=2)
    bb <- mvrnorm(n=1,mu=rep(0,2),Sigma=covMat)
    dMat <- cbind(1,scaledtime)
    gen_time_re <- dMat %*% bb
    gen_time_re_all <- c(gen_time_re_all,gen_time_re)
  }
  
  alltscaled <- (allt-ti[1])/(ti[length(ti)]-ti[1])
  mu <- mu2(allt)
  e  <- e2(sum(alln))
  y  <- mu + e + gen_time_re_all + nonparametric_re_all
  y_k2 <- data.frame(y=y,t=allt,id=rep(group2,alln))
  
  
  
  y_k1 <- merge(y_k1,D,by = "id")[,c("y","t","id","group","prob")]
  y_k2 <- merge(y_k2,D,by = "id")[,c("y","t","id","group","prob")]
  
  y_all <- rbind(y_k1,y_k2)
  
 
  
  return(list(D=D,y_all=y_all,n1=n1,n2=n2,ti=ti,tiscaled=tiscaled,
              mu1=mu1,mu2=mu2,
              mu1_true=mu1_true,mu2_true=mu2_true))
}






# algorithm
ALGORITHM <- function(hypers,DATA){
  library(dplyr)
  
  
  y_all = DATA$y_all
  D     = DATA$D
  ti    = DATA$ti
  tiscaled = DATA$tiscaled
  
  ############################################################################################################################
  #################################################### Assign to 2 groups ####################################################
  ############################################################################################################################
  
  ### randomly split into 2 groups
  random_indx <- sort(sample(1:nrow(D),round(nrow(D)/2),replace = F) )
  row_indx <- y_all$id %in% random_indx
  y1d <- y_all[row_indx,]
  y2d <- y_all[!row_indx,]
  
  y1d <- y1d[order(y1d$id),]
  y2d <- y2d[order(y2d$id),]
  
  options(dplyr.summarise.inform = F)
  GB1 <- y1d %>% group_by(id,group) %>% dplyr::summarise(n=n()) %>% ungroup()
  GB2 <- y2d %>% group_by(id,group) %>% dplyr::summarise(n=n()) %>% ungroup()
  
  true_group <- c(GB1$group,GB2$group)
  
  yalld <- rbind(y1d,y2d) # all patients
  ids <- unique(yalld$id)
  n = length(ids)
  y1d$id <- as.factor(y1d$id)
  y2d$id <- as.factor(y2d$id)
  
  
  if( sum(ids != c(GB1$id,GB2$id)) ){
    stop("ids is different from true label's id")
  }
  
  
  colnames(y1d)[2] <- "ti"
  colnames(y2d)[2] <- "ti"
  colnames(yalld)[2] <- "ti"
  
  
  
  ############################################################################################################################
  ###################################################### EM - algorithm ######################################################
  ############################################################################################################################
  ### vectors for storing estimated parameter values
  p1_r <- c()
  p2_r <- c()
  sigma2_r <- c()
  sigma2_a_1_r <- c()
  sigma2_a_2_r <- c()
  sigma2_as_1_r <- c()
  sigma2_as_2_r <- c()
  sigma2_s_1_r <- c()
  sigma2_s_2_r <- c()
  sigma2_np_1_r <- c()
  sigma2_np_2_r <- c()
  w_1_r <- c()
  w_2_r <- c()
  
  beta_vp0_r  <- c(0) # they do not have initial EM measured value, manually add one
  beta_vp1_r  <- c(0) # they do not have initial EM measured value, manually add one
  beta_vp2_r  <- c(0) # they do not have initial EM measured value, manually add one
  beta_vp3_r  <- c(0) # they do not have initial EM measured value, manually add one
  beta_vp4_r  <- c(0) # they do not have initial EM measured value, manually add one
  beta_vp5_r  <- c(0) # they do not have initial EM measured value, manually add one
  beta_vp6_r  <- c(0) # they do not have initial EM measured value, manually add one
  beta_vp7_r  <- c(0) # they do not have initial EM measured value, manually add one
  beta_vp8_r  <- c(0) # they do not have initial EM measured value, manually add one
  beta_vp9_r  <- c(0) # they do not have initial EM measured value, manually add one
  beta_vp10_r  <- c(0) # they do not have initial EM measured value, manually add one

  
  fitgss1_c_r <- c(0) # they do not have initial measured value, manually add one
  fitgss2_c_r <- c(0) # they do not have initial measured value, manually add one
  
  max_iter <- hypers$max_itero      # randomly assign groups require more steps to converge
  max_rel_change <- hypers$max_rel_changeo
  
  
  ##########################
  ### E-step: first time ###
  ##########################
  
  ### initial p_ki - simple proportion in each group
  p1 <- length(unique(y1d$id))/(length(unique(y1d$id))+length(unique(y2d$id))) # proportion of people never reached 99.9 OR group 1
  p2 <- length(unique(y2d$id))/(length(unique(y1d$id))+length(unique(y2d$id))) # proportion of people ever reached 99.9  OR group 2
  p1 <- rep(p1,n)
  p2 <- rep(p2,n)
  
  ### E-step requires an initialization of variance component and mean functions in each group: fit for each group
  ### if use_fit_init=T: we use approximation for smoothing spline & approximation for nonparametric random effect to est. var. comp. and mean func.
  ### if use_fit_init=F: we use approximation for smoothing spline & some guess for part of var. comp. (faster but less accurate)
  y1d$all <- 1
  y2d$all <- 1
  y1d$ti <- (y1d$ti-ti[1])/(ti[length(ti)]-ti[1])
  y2d$ti <- (y2d$ti-ti[1])/(ti[length(ti)]-ti[1])
  yalld$ti <- (yalld$ti-ti[1])/(ti[length(ti)]-ti[1])
  ### approx for smoothing spline
  # pre-selected points
  xg <- seq(0,1,length.out = 1000)
  # RK1 at selected points
  SigmaAt<-cubic(xg)
  ee <- eigen(SigmaAt)
  # num of eigen vectors to use
  K <- 50
  # approx Z for smoothing spline
  Z1 <- t(sqrt(1/ee$values[1:K])*t(cubic(y1d$ti,xg)%*%ee$vectors[,1:K]))
  Z2 <- t(sqrt(1/ee$values[1:K])*t(cubic(y2d$ti,xg)%*%ee$vectors[,1:K]))
  
  
  ### fit for group 1
  yo    <- y1d$y
  tio   <- y1d$ti
  allo  <- y1d$all
  ido   <- y1d$id
  
  if(length(unique(y1d$ti)) == length(ti)){
    
    fit1 <- try(lme(yo~tio, random=list(allo=pdIdent(~Z1-1),ido=~1),method = "REML", control = lmeControl(opt = "optim")))
    if(class(fit1)!="lme"){
      fit1 <- lme(yo~tio, random=list(allo=pdIdent(~Z1-1),ido=~1),method = "REML")
    }
    
    pred_mu1_t <- predict(fit1,level=0:2)[,4]
    ti_indx     <- c()
    for(tii in tiscaled){
      ti_indx   <- c(ti_indx,which(y1d$ti==tii)[1])
    }
    pred_mu1_t  <- pred_mu1_t[ti_indx]
    pred_mu1_df <- data.frame(ti=tiscaled,predValue=pred_mu1_t)
    pred_mu1_f  <- yalld
    pred_mu1_ff <- merge(x = pred_mu1_f,y = pred_mu1_df,by = "ti",all.x = T)
    
    #pred_mu1    <- pred_mu1_ff[order(pred_mu1_ff$id),"predValue"]
    
    pred_mu1 <- data.frame()
    for(idss in ids){
      tempdatause <- pred_mu1_ff[pred_mu1_ff$id==idss,c("id","ti","y","predValue")]
      tempdatause <- tempdatause[order(tempdatause$ti),]
      pred_mu1 <- rbind(pred_mu1,tempdatause)
    }
    pred_mu1 <- pred_mu1$predValue
    
    
    sigma2_1 <- fit1$sigma^2
    sigma2_a_1 <- try(as.numeric(intervals(fit1)$reStruct$id["est."])^2,silent = T)
    if(class(sigma2_a_1)!="numeric"){
      sigma2_a_1 <- as.numeric(VarCorr(fit1)["(Intercept)","StdDev"])^2
    }
    
  }else{
    print("warning: slow algorithm used, stop and check...")
    fit1 <- slm(y~ti,rk=list(cubic(ti)),random=list(id=~1),data=y1d)
    sigma2_1 <- fit1$lme.obj$sigma^2
    sigma2_a_1 <- as.numeric(substr(capture.output(fit1)[13],8,19))^2
    pred_mu1 <- predict(fit1,newdata = yalld)[,1]
  }
  
  # guess for the var. of random slope and nonparametric random effect
  sigma2_s_1 <- sigma2_a_1
  sigma2_as_1 <- 0
  sigma2_np_1 <- 1
  
  
  ### fit for group 2
  yo    <- y2d$y
  tio   <- y2d$ti
  allo  <- y2d$all
  ido   <- y2d$id
  
  if(length(unique(y2d$ti)) == length(ti)){
    
    fit2 <- try(lme(yo~tio, random=list(allo=pdIdent(~Z2-1),ido=~1),method = "REML", control = lmeControl(opt = "optim")))
    if(class(fit2)!="lme"){
      fit2 <- lme(yo~tio, random=list(allo=pdIdent(~Z2-1),ido=~1),method = "REML")
    }
    
    pred_mu2_t <- predict(fit2,level=0:2)[,4]
    ti_indx     <- c()
    for(tii in tiscaled){
      ti_indx   <- c(ti_indx,which(y2d$ti==tii)[1])
    }
    pred_mu2_t  <- pred_mu2_t[ti_indx]
    pred_mu2_df <- data.frame(ti=tiscaled,predValue=pred_mu2_t)
    pred_mu2_f  <- yalld
    pred_mu2_ff <- merge(x = pred_mu2_f,y = pred_mu2_df,by = "ti",all.x = T)
    
    #pred_mu2    <- pred_mu2_ff[order(pred_mu2_ff$id),"predValue"]
    
    pred_mu2 <- data.frame()
    for(idss in ids){
      tempdatause <- pred_mu2_ff[pred_mu2_ff$id==idss,c("id","ti","y","predValue")]
      tempdatause <- tempdatause[order(tempdatause$ti),]
      pred_mu2 <- rbind(pred_mu2,tempdatause)
    }
    pred_mu2 <- pred_mu2$predValue
    
    
    sigma2_2 <- fit2$sigma^2
    sigma2_a_2 <- try(as.numeric(intervals(fit2)$reStruct$id["est."])^2,silent = T)
    if(class(sigma2_a_2)!="numeric"){
      sigma2_a_2 <- as.numeric(VarCorr(fit2)["(Intercept)","StdDev"])^2
    }
  }else{
    print("warning: slow algorithm used, stop and check...")
    fit2 <- slm(y~ti,rk=list(cubic(ti)),random=list(id=~1),data=y2d)
    sigma2_2 <- fit2$lme.obj$sigma^2
    sigma2_a_2 <- as.numeric(substr(capture.output(fit2)[13],8,19))^2
    pred_mu2 <- predict(fit2,newdata = yalld)[,1]
  }
  
  sigma2_s_2 <- sigma2_a_2
  sigma2_as_2 <- 0
  sigma2_np_2 <- 1
  
  sigma2 = mean(c(sigma2_1,sigma2_2))
  
  ### store the estimations
  p1_r <- cbind(p1_r,p1)
  p2_r <- cbind(p2_r,p2)
  sigma2_r <- c(sigma2_r,sigma2)
  sigma2_a_1_r <- c(sigma2_a_1_r,sigma2_a_1)
  sigma2_a_2_r <- c(sigma2_a_2_r,sigma2_a_2)
  sigma2_as_1_r <- c(sigma2_as_1_r,sigma2_as_1)
  sigma2_as_2_r <- c(sigma2_as_2_r,sigma2_as_2)
  sigma2_s_1_r <- c(sigma2_s_1_r,sigma2_s_1)
  sigma2_s_2_r <- c(sigma2_s_2_r,sigma2_s_2)
  sigma2_np_1_r <- c(sigma2_np_1_r,sigma2_np_1)
  sigma2_np_2_r <- c(sigma2_np_2_r,sigma2_np_2)
  
  
  
  
  fittedgss1 = pred_mu1 # in the EM iteration, need this name
  fittedgss2 = pred_mu2
  
  
  ##################################################
  # we already have an initial estimate of var. comp.
  # random error: sigma2_1, sigma2_2
  # random intercept: sigma2_a_1, sigma2_a_2
  # random slope: sigma2_s_1, sigma2_s_2
  # nonpara rand eff: sigma2_np_1, sigma2_np_2
  ##################################################
  
  tiid <- table(yalld$id)
  tiid <- tiid[match(ids,names(tiid))] # table function automatically sorts the result, need to re-order
  tiidc <- cumsum(tiid)
  
  # ZmatL matrix used through out the later optimization, won't change through each iteration
  # GL matrix is the G matrix but with sigma2's as 1. This is used in opt for var. comp. and won't change for each subject
  ZmatL <- vector("list",length(ids)) # used later in optimization for var. comp.
  GL <- vector("list",length(ids))    # used later in optimization for var. comp.

  
  for(i in 1:length(ids)){
    alltempdata <- yalld[yalld$id == ids[i] ,]
    nelem <- nrow(alltempdata)
    # Z
    Zmat <- cbind( rep(1,nelem), alltempdata$ti, diag(nelem) )
    ZmatL[[i]] <- Zmat
    # G
    Gappend <- as.matrix(bdiag(list(matrix(0,nrow=2,ncol=2),cubic(alltempdata$ti))))
    GL[[i]] <- Gappend
  }
  
  
  
  # construct V1 and V2 - used in E_step to estimate w_1 and w_2
  # construct G1 and G2 - used in M_step to estimate var. comp.
  
  V1 <- vector("list",length(ids))
  V2 <- vector("list",length(ids))
  
  
  for(i in 1:length(ids)){
    alltempdata <- yalld[yalld$id == ids[i] ,]
    nelem <- nrow(alltempdata)
    # V1
    LR1 <- sigma2_np_1 * cubic(alltempdata$ti)
    UL1 <- matrix(c(sigma2_a_1,sigma2_as_1,sigma2_as_1,sigma2_s_1),nrow=2)
    G1 <- as.matrix(bdiag(list(UL1,LR1)))
    V1[[i]] <- ZmatL[[i]] %*% G1 %*% t(ZmatL[[i]]) + sigma2*diag(nelem)
    # V2
    LR2 <- sigma2_np_2 * cubic(alltempdata$ti)
    UL2 <- matrix(c(sigma2_a_2,sigma2_as_2,sigma2_as_2,sigma2_s_2),nrow=2)
    G2 <- as.matrix(bdiag(list(UL2,LR2)))
    V2[[i]] <- ZmatL[[i]] %*% G2 %*% t(ZmatL[[i]]) + sigma2*diag(nelem)
  }
  
  
  
  
  
  n_iter <- 1
  cur_change <- max_rel_change + 1 # make sure at least run once
  
  ## logistic regression - use glmnet
  # random samples for several times, each time fit a glmnet, take the average of beta's
  nFolds <- 10
  foldids <- sample(rep(seq(nFolds),length.out = length(ids)))

  while( (n_iter <= max_iter) & (cur_change > max_rel_change) ){
    
    if(n_iter == 1){
      print(paste0("max_iter: ", max_iter, " || max_rel_diff: ",max_rel_change))
    }
    
    #########################
    ### E-step: find w_ik ###
    #########################
    # only for numerical reasons, p1 and p2 cannot be strictly 0/1
    p1 = ifelse(p1 < 1e-80, 1e-80, p1)
    p2 = 1 - p1
    p2 = ifelse(p2 < 1e-80, 1e-80, p2) # in case some p1 strictly 1
    
    ws <- rcpp_update_w(as.numeric(tiidc), ids, fittedgss1, fittedgss2, p1, p2, yalld$y, yalld$ti, V1, V2)
    
    w_1 <- as.numeric(ws$w_1) # w_i1
    w_2 <- as.numeric(ws$w_2) # w_i2
    
    w_1_r <- cbind(w_1_r,w_1)
    w_2_r <- cbind(w_2_r,w_2)
    
    
    X <- data.frame(w_1=w_1,vp1=D$vp1[ids],vp2=D$vp2[ids],vp3=D$vp3[ids],vp4=D$vp4[ids],vp5=D$vp5[ids],
                    vp6=D$vp6[ids],vp7=D$vp7[ids],vp8=D$vp8[ids],vp9=D$vp9[ids],vp10=D$vp10[ids])
    tempfit <- glm(w_1~.,data=X, family="binomial")
    coefs <- coef(tempfit)
    
    Xm <- model.matrix(tempfit)[,-1]
    
    # responseY=as.integer(10000*w_1)
    # responseY0 = 10000-responseY
    # cvfit=cv.glmnet(Xm, cbind(responseY0,responseY), family="binomial",alpha=1,nfolds = 10,maxit=1.5e+05,foldid = foldids)
    # coefs <- coef(cvfit, s = "lambda.min")
    
    vp0hat <- coefs[1]
    vp1hat <- coefs[2]
    vp2hat <- coefs[3]
    vp3hat <- coefs[4]
    vp4hat <- coefs[5]
    vp5hat <- coefs[6]
    vp6hat <- coefs[7]
    vp7hat <- coefs[8]
    vp8hat <- coefs[9]
    vp9hat <- coefs[10]
    vp10hat <- coefs[11]

    
    
    betashat <- c(vp0hat,vp1hat,vp2hat,vp3hat,vp4hat,vp5hat,vp6hat,vp7hat,vp8hat,vp9hat,vp10hat)

    checkbetas = cbind(rep(1,nrow(Xm)),Xm) %*% betashat
    usedbetas = ifelse(abs(checkbetas)>=20,20*sign(checkbetas),checkbetas)
    
    p1 <- exp( usedbetas ) /  ( 1+exp( usedbetas ) )
    p2 <- 1 - p1
    
    
    
    
    #################################################################
    ### M-step: iterative estimation of var. comp. and mean func. ###
    #################################################################
    
    # ssanova use approximation by randomly choosing id.basis
    # in our case, the H1 is spanned by those 61 (or other chosen numbers) time points we have
    # so we do not need to randomly select id.basis, we can find id.basis corresponding to 61 time points
    
    id.basis_indx     <- c()
    for(tii in tiscaled){
      id.basis_indx   <- c(id.basis_indx,which(yalld$ti==tii)[1])
    }
    
    
    # inner loop convergence criteria
    max_iter_inner <- hypers$max_iter_innero      # the estimate of c and d in smoothing spline does not converge
    max_rel_change_inner <- hypers$max_rel_change_innero
    n_iter_inner <- 1
    
    
    sigma2_inner_r <- c(sigma2)
    sigma2_a_1_inner_r <- c(sigma2_a_1)
    sigma2_a_2_inner_r <- c(sigma2_a_2)
    sigma2_as_1_inner_r <- c(sigma2_as_1)
    sigma2_as_2_inner_r <- c(sigma2_as_2)
    sigma2_s_1_inner_r <- c(sigma2_s_1)
    sigma2_s_2_inner_r <- c(sigma2_s_2)
    sigma2_np_1_inner_r <- c(sigma2_np_1)
    sigma2_np_2_inner_r <- c(sigma2_np_2)
    if(n_iter==1){
      fitgss1_c_inner_r <- c(0) # they do not have initial measured value, manually add one
      fitgss2_c_inner_r <- c(0) # they do not have initial measured value, manually add one
    }else{
      fitgss1_c_inner_r <- c(fitgss1_c)
      fitgss2_c_inner_r <- c(fitgss2_c)
    }
    
    
    cur_change_inner <- max_rel_change_inner + 1 # make sure at least run once
    
    
    while( (n_iter_inner <= max_iter_inner) & (cur_change_inner > max_rel_change_inner) ){
      
      if(n_iter_inner == 1){
        print(paste0("max_iter_inner: ", max_iter_inner, " || max_rel_diff_inner: ",max_rel_change_inner))
      }
      
      ###------------------------------###
      ### fixed var. comp. and c and d ###
      ###------------------------------###
      
      # GROUP1: random intercept, random slope and nonparametric random effect
      WVW_all_g1l <- vector("list",n)
      for(i in 1:n){
        WVW_all_g1l[[i]] <- V1[[i]] / w_1[i]
      }
      fitgss1 <- ssanova9( y~ti,data=yalld,cov=list("known",WVW_all_g1l),id.basis = id.basis_indx, method="m", tiidin = tiid, tiidcin = tiidc)
      fittedgss1 <- fitted(fitgss1)
      
      
      
      # GROUP2: random intercept, random slope and nonparametric random effect
      WVW_all_g2l <- vector("list",n)
      for(i in 1:n){
        WVW_all_g2l[[i]] <- V2[[i]] / w_2[i]
      }
      fitgss2 <- ssanova9( y~ti,data=yalld,cov=list("known",WVW_all_g2l),id.basis = id.basis_indx, method="m", tiidin = tiid, tiidcin = tiidc)
      fittedgss2 <- fitted(fitgss2)
      
      
      fitgss1_c <- predict.ssanova(fitgss1,newdata=data.frame(ti=tiscaled),estSigma2in = sigma2)
      fitgss2_c <- predict.ssanova(fitgss2,newdata=data.frame(ti=tiscaled),estSigma2in = sigma2)
      
      ###---------------------------------------------###
      ### fixed c and d, estimate variance components ###
      ###---------------------------------------------###
      
      Y1 <- yalld$y - fittedgss1
      Y2 <- yalld$y - fittedgss2
      

      optr <- try(optim(par=rep(0,8),
                        fn=lp4V2,
                        gr=lp4_grV2,
                        tiid_in=tiid, tiidc_in=tiidc, Y1_in=Y1, Y2_in=Y2,
                        w1_in=w_1, w2_in=w_2, G1_in=GL, G2_in=GL, Z_in=ZmatL, E_in=GL,
                        method="L-BFGS-B",
                        upper=c(20,100,20,20,20,100,20,20),
                        control=list(maxit=10000,factr=1e2)))
      
      
      
      if(class(optr)=="try-error"){
        optr <- optim(par=rep(0,8),
                      fn=lp4V2,
                      gr=lp4_grV2,
                      tiid_in=tiid, tiidc_in=tiidc, Y1_in=Y1, Y2_in=Y2,
                      w1_in=w_1, w2_in=w_2, G1_in=GL, G2_in=GL, Z_in=ZmatL, E_in=GL,
                      method="L-BFGS-B",
                      upper=c( 20, 100, 20, 20, 20, 100, 20, 20),
                      lower=c(-30,-100,-30,-30,-30,-100,-30,-30),
                      control=list(maxit=10000,factr=1e2))
      }
      
     
      
      check_gr = lp4_grV2(optr$par,tiid_in=tiid, tiidc_in=tiidc, Y1_in=Y1, Y2_in=Y2,
                          w1_in=w_1, w2_in=w_2, G1_in=GL, G2_in=GL, Z_in=ZmatL,E_in = GL)
      

      Eta1 <- exp(optr$par[1])
      Eta2 <- 2/pi * atan(optr$par[2])
      Eta3 <- exp(optr$par[3])
      Eta4 <- exp(optr$par[4])
      
      Eta5 <- exp(optr$par[5])
      Eta6 <- 2/pi * atan(optr$par[6])
      Eta7 <- exp(optr$par[7])
      Eta8 <- exp(optr$par[8])
      
      
      
      # use evaluate function lpev to get estiamte of sigma
      # now G1 and G2 are the estimated G1 and G2
      G1EL <- vector("list",n) # used for evaluating sigma2
      UL1 <- matrix(c(Eta1,Eta2*sqrt(Eta1)*sqrt(Eta3),Eta2*sqrt(Eta1)*sqrt(Eta3),Eta3),ncol=2) 
      for(i in 1:n){
        alltempdata <- yalld[yalld$id == ids[i] ,]
        LR1 <- Eta4 * cubic(alltempdata$ti)
        G1 <- as.matrix(bdiag(list(UL1,LR1)))
        G1EL[[i]] <- G1
      }
      
      G2EL <- vector("list",n) # used for evaluating sigma2
      UL2 <- matrix(c(Eta5,Eta6*sqrt(Eta5)*sqrt(Eta7),Eta6*sqrt(Eta5)*sqrt(Eta7),Eta7),ncol=2)
      for(i in 1:n){
        alltempdata <- yalld[yalld$id == ids[i] ,]
        LR2 <- Eta8 * cubic(alltempdata$ti)
        G2 <- as.matrix(bdiag(list(UL2,LR2)))
        G2EL[[i]] <- G2
      }
      
      
      sigma2 <- lpev3(tiid, tiidc, Y1_in=Y1, Y2_in=Y2, w1_in=w_1, w2_in=w_2, G1_in=G1EL, G2_in=G2EL, Z_in=ZmatL)
      sigma2_a_1 <- as.numeric(UL1[1,1] * sigma2)
      sigma2_as_1 <- as.numeric(UL1[2,1] * sigma2)
      sigma2_s_1 <- as.numeric(UL1[2,2] * sigma2)
      sigma2_np_1 <- as.numeric(Eta4 * sigma2)
      sigma2_a_2 <- as.numeric(UL2[1,1] * sigma2)
      sigma2_as_2 <- as.numeric(UL2[2,1] * sigma2)
      sigma2_s_2 <- as.numeric(UL2[2,2] * sigma2)
      sigma2_np_2 <- as.numeric(Eta8 * sigma2)
      
      
      
      
      ###----------------------------------------------------------###
      ### now we need to update V1 and V2 for next inner iteration ###
      ###----------------------------------------------------------###
      
      V1 <- vector("list",n)
      V2 <- vector("list",n)
      
      for(i in 1:length(ids)){
        alltempdata <- yalld[yalld$id == ids[i] ,]
        nelem <- nrow(alltempdata)
        # V1
        LR1 <- sigma2_np_1 * cubic(alltempdata$ti)
        UL1 <- matrix(c(sigma2_a_1,sigma2_as_1,sigma2_as_1,sigma2_s_1),nrow=2)
        G1 <- as.matrix(bdiag(list(UL1,LR1)))
        V1[[i]] <- ZmatL[[i]] %*% G1 %*% t(ZmatL[[i]]) + sigma2*diag(nelem)
        # V2
        LR2 <- sigma2_np_2 * cubic(alltempdata$ti)
        UL2 <- matrix(c(sigma2_a_2,sigma2_as_2,sigma2_as_2,sigma2_s_2),nrow=2)
        G2 <- as.matrix(bdiag(list(UL2,LR2)))
        V2[[i]] <- ZmatL[[i]] %*% G2 %*% t(ZmatL[[i]]) + sigma2*diag(nelem)
      }
      
      
      ###-----------------------###
      ### store the estimations ###
      ###-----------------------###
      
      sigma2_inner_r <- c(sigma2_inner_r,sigma2)
      sigma2_a_1_inner_r <- c(sigma2_a_1_inner_r,sigma2_a_1)
      sigma2_a_2_inner_r <- c(sigma2_a_2_inner_r,sigma2_a_2)
      sigma2_as_1_inner_r <- c(sigma2_as_1_inner_r,sigma2_as_1)
      sigma2_as_2_inner_r <- c(sigma2_as_2_inner_r,sigma2_as_2)
      sigma2_s_1_inner_r <- c(sigma2_s_1_inner_r,sigma2_s_1)
      sigma2_s_2_inner_r <- c(sigma2_s_2_inner_r,sigma2_s_2)
      sigma2_np_1_inner_r <- c(sigma2_np_1_inner_r,sigma2_np_1)
      sigma2_np_2_inner_r <- c(sigma2_np_2_inner_r,sigma2_np_2)
      fitgss1_c_inner_r <- cbind(fitgss1_c_inner_r,fitgss1_c)
      fitgss2_c_inner_r <- cbind(fitgss2_c_inner_r,fitgss2_c)
      
      len_inner <- length(sigma2_a_1_inner_r)
      
      ###-----------------###
      ### relative change ###
      ###-----------------###
      
      ThetaNew <- c(sigma2_inner_r[len_inner],
                    sigma2_a_1_inner_r[len_inner],sigma2_a_2_inner_r[len_inner],
                    sigma2_as_1_inner_r[len_inner],sigma2_as_2_inner_r[len_inner],
                    sigma2_s_1_inner_r[len_inner],sigma2_s_2_inner_r[len_inner],
                    sigma2_np_1_inner_r[len_inner],sigma2_np_2_inner_r[len_inner],
                    fitgss1_c_inner_r[,len_inner],
                    fitgss2_c_inner_r[,len_inner])
      
      ThetaOld <- c(sigma2_inner_r[len_inner-1],
                    sigma2_a_1_inner_r[len_inner-1],sigma2_a_2_inner_r[len_inner-1],
                    sigma2_as_1_inner_r[len_inner-1],sigma2_as_2_inner_r[len_inner-1],
                    sigma2_s_1_inner_r[len_inner-1],sigma2_s_2_inner_r[len_inner-1],
                    sigma2_np_1_inner_r[len_inner-1],sigma2_np_2_inner_r[len_inner-1],
                    fitgss1_c_inner_r[,len_inner-1],
                    fitgss2_c_inner_r[,len_inner-1])
      
      
      
      cur_change_inner <- sum((ThetaNew - ThetaOld)^2)/(sum(ThetaOld^2)+0.00001)
      
      print(paste0("      (iter finished) cur_iter_inner: ",n_iter_inner, " || cur_diff_inner: ",round(cur_change_inner,9)))
      
      n_iter_inner <- n_iter_inner + 1
      
    }
    
    
    ### store the estimations
    p1_r <- cbind(p1_r,p1)
    p2_r <- cbind(p2_r,p2)
    beta_vp0_r <- c(beta_vp0_r,vp0hat)
    beta_vp1_r <- c(beta_vp1_r,vp1hat)
    beta_vp2_r <- c(beta_vp2_r,vp2hat)
    beta_vp3_r <- c(beta_vp3_r,vp3hat)
    beta_vp4_r <- c(beta_vp4_r,vp4hat)
    beta_vp5_r <- c(beta_vp5_r,vp5hat)
    beta_vp6_r <- c(beta_vp6_r,vp6hat)
    beta_vp7_r <- c(beta_vp7_r,vp7hat)
    beta_vp8_r <- c(beta_vp8_r,vp8hat)
    beta_vp9_r <- c(beta_vp9_r,vp9hat)
    beta_vp10_r <- c(beta_vp10_r,vp10hat)

    
    sigma2_r <- c(sigma2_r,sigma2)
    sigma2_a_1_r <- c(sigma2_a_1_r,sigma2_a_1)
    sigma2_a_2_r <- c(sigma2_a_2_r,sigma2_a_2)
    sigma2_as_1_r <- c(sigma2_as_1_r,sigma2_as_1)
    sigma2_as_2_r <- c(sigma2_as_2_r,sigma2_as_2)
    sigma2_s_1_r <- c(sigma2_s_1_r,sigma2_s_1)
    sigma2_s_2_r <- c(sigma2_s_2_r,sigma2_s_2)
    sigma2_np_1_r <- c(sigma2_np_1_r,sigma2_np_1)
    sigma2_np_2_r <- c(sigma2_np_2_r,sigma2_np_2)
    fitgss1_c_r <- cbind(fitgss1_c_r,fitgss1_c)
    fitgss2_c_r <- cbind(fitgss2_c_r,fitgss2_c)
    
    
    len <- length(sigma2_a_1_r)
    
    # relative change
    
    ThetaNew <- c(sigma2_r[len],
                  sigma2_a_1_r[len],sigma2_a_2_r[len],
                  sigma2_s_1_r[len],sigma2_s_2_r[len],
                  sigma2_as_1_r[len],sigma2_as_2_r[len],
                  sigma2_np_1_r[len],sigma2_np_2_r[len],
                  fitgss1_c_r[,len],
                  fitgss2_c_r[,len],
                  beta_vp0_r[len],beta_vp1_r[len],beta_vp2_r[len],beta_vp3_r[len],beta_vp4_r[len],
                  beta_vp5_r[len],beta_vp6_r[len],beta_vp7_r[len],beta_vp8_r[len],beta_vp9_r[len],
                  beta_vp10_r[len])
    
    ThetaOld <- c(sigma2_r[len-1],
                  sigma2_a_1_r[len-1],sigma2_a_2_r[len-1],
                  sigma2_s_1_r[len-1],sigma2_s_2_r[len-1],
                  sigma2_as_1_r[len-1],sigma2_as_2_r[len-1],
                  sigma2_np_1_r[len-1],sigma2_np_2_r[len-1],
                  fitgss1_c_r[,len-1],
                  fitgss2_c_r[,len-1],
                  beta_vp0_r[len-1],beta_vp1_r[len-1],beta_vp2_r[len-1],beta_vp3_r[len-1],beta_vp4_r[len-1],
                  beta_vp5_r[len-1],beta_vp6_r[len-1],beta_vp7_r[len-1],beta_vp8_r[len-1],beta_vp9_r[len-1],
                  beta_vp10_r[len-1])
    
    
    cur_change <- sum((ThetaNew - ThetaOld)^2)/(sum(ThetaOld^2)+0.00001)
 
    print(paste0("(iter finished) cur_iter: ",n_iter, " || cur_diff: ",round(cur_change,16)))
    
    n_iter <- n_iter + 1
    
    
  }
  
  return(list(w_1_r=w_1_r,true_group=true_group,yalld=yalld,
              ###
              fitgss1=fitgss1,fitgss2=fitgss2,
              ###
              sigma2=sigma2,sigma2_a_1=sigma2_a_1,sigma2_a_2=sigma2_a_2,
              sigma2_as_1=sigma2_as_1,sigma2_as_2 = sigma2_as_2,sigma2_s_1 = sigma2_s_1,
              sigma2_s_2 = sigma2_s_2,sigma2_np_1 = sigma2_np_1,sigma2_np_2 = sigma2_np_2,
              ###
              vp0hat=vp0hat,vp1hat=vp1hat,vp2hat=vp2hat,vp3hat=vp3hat,
              vp4hat=vp4hat,vp5hat=vp5hat,vp6hat=vp6hat,vp7hat=vp7hat,
              vp8hat=vp8hat,vp9hat=vp9hat,vp10hat=vp10hat,
              ###
              betashat=betashat,p1=p1,p2=p2,ids=ids))
  
}













# clean results
RESULTS <- function(hypers,DATA,EMresults,tm.all){
  ti=DATA$ti
  tiscaled=DATA$tiscaled
  y_all = DATA$y_all
  ids = EMresults$ids 
  D = DATA$D
  
  sigma2_1_true = hypers$sigma2_1_trueo
  sigma2_2_true = hypers$sigma2_2_trueo
  sigma2_a_1_true = hypers$sigma2_a_1_trueo
  sigma2_a_2_true = hypers$sigma2_a_2_trueo
  sigma2_as_1_true = hypers$sigma2_as_1_trueo
  sigma2_as_2_true = hypers$sigma2_as_2_trueo
  sigma2_s_1_true = hypers$sigma2_s_1_trueo
  sigma2_s_2_true = hypers$sigma2_s_2_trueo
  sigma2_np_1_true = hypers$sigma2_np_1_trueo
  sigma2_np_2_true = hypers$sigma2_np_2_trueo
  
  beta_vp0 <- hypers$beta_vp0o
  beta_vp1 <- hypers$beta_vp1o
  beta_vp2 <- hypers$beta_vp2o
  beta_vp3 <- hypers$beta_vp3o
  beta_vp4 <- hypers$beta_vp4o
  beta_vp5 <- hypers$beta_vp5o
  beta_vp6 <- hypers$beta_vp6o       
  beta_vp7 <- hypers$beta_vp7o        
  beta_vp8 <- hypers$beta_vp8o   
  beta_vp9 <- hypers$beta_vp9o
  beta_vp10 <- hypers$beta_vp10o

  w_1_r = w_1 = EMresults$w_1_r
  w_2_r = w_2 = 1 - w_1
  true_group = EMresults$true_group
  yalld = EMresults$yalld
  
  fitgss1=EMresults$fitgss1
  fitgss2=EMresults$fitgss2
  
  sigma2=EMresults$sigma2
  sigma2_a_1=EMresults$sigma2_a_1
  sigma2_a_2=EMresults$sigma2_a_2
  sigma2_as_1=EMresults$sigma2_as_1
  sigma2_as_2 = EMresults$sigma2_as_2
  sigma2_s_1 = EMresults$sigma2_s_1
  sigma2_s_2 = EMresults$sigma2_s_2
  sigma2_np_1 = EMresults$sigma2_np_1
  sigma2_np_2 = EMresults$sigma2_np_2
  
  vp0hat=EMresults$vp0hat
  vp1hat=EMresults$vp1hat
  vp2hat=EMresults$vp2hat
  vp3hat=EMresults$vp3hat
  vp4hat=EMresults$vp4hat
  vp5hat=EMresults$vp5hat
  vp6hat=EMresults$vp6hat
  vp7hat=EMresults$vp7hat
  vp8hat=EMresults$vp8hat
  vp9hat=EMresults$vp9hat
  vp10hat=EMresults$vp10hat
  betashat=EMresults$betashat
  
  p1=EMresults$p1
  p2=EMresults$p2
  
  mu1 = DATA$mu1
  mu2 = DATA$mu2
  
  n1 = DATA$n1
  n2 = DATA$n2
  
  ### clustering confusion matrix
  ### cluster according to threshold 0.5
  cluster <- ifelse(w_1_r[,ncol(w_1_r)]>=0.5,1,2)
  conf_mat <- table(true_group,cluster)
  clus_accuracy <- sum(diag(conf_mat))/sum(conf_mat)
  sensitivity <- conf_mat[1,1]/sum(conf_mat[1,])
  specificity <- conf_mat[2,2]/sum(conf_mat[2,])
  true_group_prob <- (yalld %>% group_by(id,group,prob) %>% dplyr::summarise(n=n()))[,1:3]
  
  pred_mu1 <- predict.ssanova(fitgss1,newdata=data.frame(ti=tiscaled),estSigma2in = sigma2)
  pred_mu2 <- predict.ssanova(fitgss2,newdata=data.frame(ti=tiscaled),estSigma2in = sigma2)
  
  if( mean(pred_mu1[29:31]) < mean(pred_mu2[29:31])){ # we used the wrong label, flipped the label
    cluster <- ifelse(w_1_r[,ncol(w_1_r)]>=0.5,2,1)
    conf_mat <- table(true_group,cluster)
    clus_accuracy <- sum(diag(conf_mat))/sum(conf_mat)
    sensitivity <- conf_mat[1,1]/sum(conf_mat[1,])
    specificity <- conf_mat[2,2]/sum(conf_mat[2,])
    
    # we also need to flip every estimates, i.e., sigma2_a_1 is actually sigma2_a_2
    estSigma2    <- sigma2
    estSigma2_a1 <- sigma2_a_2
    estSigma2_a2 <- sigma2_a_1
    estSigma2_as1 <- sigma2_as_2
    estSigma2_as2 <- sigma2_as_1
    estSigma2_s1 <- sigma2_s_2
    estSigma2_s2 <- sigma2_s_1
    estSigma2_np1 <- sigma2_np_2
    estSigma2_np2 <- sigma2_np_1
    
    estvp0        <- -1 * vp0hat
    estvp1        <- -1 * vp1hat
    estvp2        <- -1 * vp2hat
    estvp3        <- -1 * vp3hat
    estvp4        <- -1 * vp4hat
    estvp5        <- -1 * vp5hat
    estvp6        <- -1 * vp6hat
    estvp7        <- -1 * vp7hat
    estvp8        <- -1 * vp8hat
    estvp9        <- -1 * vp9hat
    estvp10       <- -1 * vp10hat
    #estvp11       <- -1 * vp11hat
    
    estimatebetas = -1 * betashat
    
    estP1        <- p2
    estP2        <- p1
    trueP1       <- true_group_prob$prob
    trueP2       <- 1 - true_group_prob$prob
    MSEP1        <- sum((trueP1 - estP1)^2)/length(trueP1)
    MSEP2        <- sum((trueP2 - estP2)^2)/length(trueP2)
    
    ### mean function error
    pred_mu1 <- predict.ssanova(fitgss2,newdata=data.frame(ti=tiscaled),se.fit = T,estSigma2in = estSigma2)
    pred_mu2 <- predict.ssanova(fitgss1,newdata=data.frame(ti=tiscaled),se.fit = T,estSigma2in = estSigma2)
    
    p_mu1 <- pred_mu1$fit
    p_mu2 <- pred_mu2$fit
    
    ### mean function MSE
    MSEmu1 <- sum((p_mu1-mu1(ti))^2)/length(ti)
    MSEmu2 <- sum((p_mu2-mu2(ti))^2)/length(ti)
    
    ### CI coverage for mean functions
    p_mu1L <- p_mu1 - qnorm(0.975)*pred_mu1$se.fit
    p_mu1U <- p_mu1 + qnorm(0.975)*pred_mu1$se.fit
    p_mu2L <- p_mu2 - qnorm(0.975)*pred_mu2$se.fit
    p_mu2U <- p_mu2 + qnorm(0.975)*pred_mu2$se.fit
    CI1 <- rbind(p_mu1L,p_mu1,p_mu1U,mu1(ti))
    CI2 <- rbind(p_mu2L,p_mu2,p_mu2U,mu2(ti))
    Fcoverage1 <- mean((CI1[4,]>CI1[1,]) & (CI1[4,]<CI1[3,]))
    Fcoverage2 <- mean((CI2[4,]>CI2[1,]) & (CI2[4,]<CI2[3,]))
    
    
    # used for CI of var. comp.
    ui1 = w_2
    ui2 = 1-ui1
    
    
    # flip fitgss1 and fitgss2, will be used later
    fitgsstemp = fitgss1
    fitgss1 = fitgss2
    fitgss2 = fitgsstemp
    
    
    
  }else{
    # no need to flip
    estSigma2    <- sigma2
    estSigma2_a1 <- sigma2_a_1
    estSigma2_a2 <- sigma2_a_2
    estSigma2_as1 <- sigma2_as_1
    estSigma2_as2 <- sigma2_as_2
    estSigma2_s1 <- sigma2_s_1
    estSigma2_s2 <- sigma2_s_2
    estSigma2_np1 <- sigma2_np_1
    estSigma2_np2 <- sigma2_np_2
    
    estvp0        <- vp0hat
    estvp1        <- vp1hat
    estvp2        <- vp2hat
    estvp3        <- vp3hat
    estvp4        <- vp4hat
    estvp5        <- vp5hat
    estvp6        <- vp6hat
    estvp7        <- vp7hat
    estvp8        <- vp8hat
    estvp9        <- vp9hat
    estvp10       <- vp10hat
    #estvp11       <- vp11hat
    
    estimatebetas <- betashat
    
    estP1        <- p1
    estP2        <- p2
    trueP1       <- true_group_prob$prob
    trueP2       <- 1 - true_group_prob$prob
    MSEP1        <- sum((trueP1 - estP1)^2)/length(trueP1)
    MSEP2        <- sum((trueP2 - estP2)^2)/length(trueP2)
    
    ### mean function error
    pred_mu1 <- predict.ssanova(fitgss1,newdata=data.frame(ti=tiscaled),se.fit = T,estSigma2in = estSigma2)
    pred_mu2 <- predict.ssanova(fitgss2,newdata=data.frame(ti=tiscaled),se.fit = T,estSigma2in = estSigma2)
    
    
    p_mu1 <- pred_mu1$fit
    p_mu2 <- pred_mu2$fit
    
    ### mean function MSE
    MSEmu1 <- sum((p_mu1-mu1(ti))^2)/length(ti)
    MSEmu2 <- sum((p_mu2-mu2(ti))^2)/length(ti)
    
    ### CI coverage
    p_mu1L <- p_mu1 - qnorm(0.975)*pred_mu1$se.fit
    p_mu1U <- p_mu1 + qnorm(0.975)*pred_mu1$se.fit
    p_mu2L <- p_mu2 - qnorm(0.975)*pred_mu2$se.fit
    p_mu2U <- p_mu2 + qnorm(0.975)*pred_mu2$se.fit
    CI1 <- rbind(p_mu1L,p_mu1,p_mu1U,mu1(ti))
    
    CI2 <- rbind(p_mu2L,p_mu2,p_mu2U,mu2(ti))
    Fcoverage1 <- mean((CI1[4,]>CI1[1,]) & (CI1[4,]<CI1[3,]))
    Fcoverage2 <- mean((CI2[4,]>CI2[1,]) & (CI2[4,]<CI2[3,]))
    
    
    # used for CI of var. comp.
    ui1 = w_1
    ui2 = 1-ui1
    
  }
  
  
  
  result <- matrix(c(sigma2_1_true,sigma2_a_1_true,sigma2_a_2_true,sigma2_as_1_true,sigma2_as_2_true,
                     sigma2_s_1_true,sigma2_s_2_true,sigma2_np_1_true,sigma2_np_2_true,
                     beta_vp0,beta_vp1,beta_vp2,beta_vp3,beta_vp4,beta_vp5,
                     beta_vp6,beta_vp7,beta_vp8,beta_vp9,beta_vp10,
                     estSigma2,estSigma2_a1,estSigma2_a2,estSigma2_as1,estSigma2_as2,
                     estSigma2_s1,estSigma2_s2,estSigma2_np1,estSigma2_np2,
                     estvp0,estvp1,estvp2,estvp3,estvp4,estvp5,estvp6,estvp7,estvp8,estvp9,estvp10),
                   nrow=2,byrow=T,
                   dimnames = list(c("true_val","est_val"),
                                   c("sigma2","sigma2_a1","sigma2_a2","sigma2_as1","sigma2_as2",
                                     "sigma2_s1","sigma2_s2","sigma2_np1","sigma2_np2",
                                     "b_p0","b_p1","b_p2","b_p3","b_p4","b_p5","b_p6",
                                     "b_p7","b_p8","b_p9","b_p10")))
  
  
  
  
  result <- cbind(result,c(NA,clus_accuracy))
  result <- cbind(result,c(NA,sensitivity))
  result <- cbind(result,c(NA,specificity))
  result <- cbind(result,c(NA,MSEP1))
  result <- cbind(result,c(NA,MSEP2))
  result <- cbind(result,c(NA,MSEmu1))
  result <- cbind(result,c(NA,MSEmu2))
  result <- cbind(result,c(n1,sum(cluster==1)))
  result <- cbind(result,c(n2,sum(cluster==2)))
  
  
  colnames(result)[21:29] <-
    c("cluster_accu","sensitivity","specificity","MSEp1","MSEp2",
      "MSEmu1","MSEmu2","n1","n2")
  
  
  tm.all <- proc.time() - tm.all
  
  result <- cbind(result,c(NA,tm.all[3]))
  
  colnames(result)[30] <- "elapsed CPU time"
  
  
  
  RESULTS = list()
  RESULTS$basics = result
  RESULTS$pred_mu1 = pred_mu1
  RESULTS$pred_mu2 = pred_mu2
  RESULTS$w_1 = w_1_r
  
  
  saveRDS(RESULTS,"results.rds")
  return(RESULTS)
}

















