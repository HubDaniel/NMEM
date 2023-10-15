
  
mkcov.known <- function(w) {
  ## check inputs
  if(class(w)[1]=="matrix"){
    if (nrow(w)-ncol(w)) stop("gss error in mkcov.known: matrix is not square")
  }else if(class(w)[1]=="list"){
    #print("passing a list")
  }
  # otherwise, w is a list of block diagonal elements
  env <- w
  init <- NULL
  fun <- function(env) env
  list(fun=fun,env=env,init=init)
}



ssanova9 <- function(formula,type=NULL,data=list(),subset,
                     offset,na.action=na.omit,partial=NULL,
                     method="v",alpha=1.4,varht=1,
                     id.basis=NULL,nbasis=NULL,seed=NULL,cov,
                     skip.iter=FALSE,tiidin,tiidcin)
{
  ## Obtain model frame and model terms
  mf <- match.call()
  

  mf$type <- mf$method <- mf$varht <- mf$partial <- NULL
  mf$alpha <- mf$id.basis <- mf$nbasis <- mf$seed <- NULL
  mf$cov <- mf$skip.iter <- mf$tiidin <- mf$tiidcin <- NULL
  

  mf[[1]] <- as.name("model.frame")
  

  
  
  mf <- eval(mf,parent.frame())

  
  
  
  ## Generate sub-basis
  nobs <- dim(mf)[1]
  if (is.null(id.basis)) {
    if (is.null(nbasis))  nbasis <- max(30,ceiling(10*nobs^(2/9)))
    if (nbasis>=nobs)  nbasis <- nobs
    if (!is.null(seed))  set.seed(seed)
    id.basis <- sample(nobs,nbasis)
  }
  else {
    if (max(id.basis)>nobs|min(id.basis)<1)
      stop("gss error in ssanova9: id.basis out of range")
    nbasis <- length(id.basis)
  }
  
 
  
  #tm.inner.cal.gtermSR <- proc.time()
  
  ## Generate terms
  term <- mkterm(mf,type)
  
  ## Generate cov
  if (is.null(cov$fun)) {
    type <- cov[[1]]
    if (type=="arma") {
      pq <- cov[[2]]
      cov <- mkcov.arma(pq[1],pq[2],nobs)
    }
    if (type=="long") {
      if (nobs<length(cov[[2]])) id <- cov[[2]][subset]
      else id <- cov[[2]]
      cov <- mkcov.long(id)
    }
    if (type=="known") {
      if(class(cov[[2]])[1]=="matrix"){
        #print("passing in a matrix as cov")
        cov <- mkcov.known(cov[[2]])
      }
      if(class(cov[[2]])[1]=="list"){
        #print("passing in a list as cov")
        cov <- mkcov.known(cov[[2]])
      }
      
    }
  }
  
  
  
  ## Generate s and r
  s <- r <- NULL
  nq <- 0
  for (label in term$labels) {
    if (label=="1") {
      s <- cbind(s,rep(1,len=nobs))
      next
    }
    x <- mf[,term[[label]]$vlist]
    x.basis <- mf[id.basis,term[[label]]$vlist]
    nphi <- term[[label]]$nphi # basis functions: 1 for ti
    nrk <- term[[label]]$nrk # reproducing kernel: 1 for ti
    
    # print(mf)
    # print("=======================")
    
    if (nphi) {
      phi <- term[[label]]$phi
      for (i in 1:nphi)
        s <- cbind(s,phi$fun(x,nu=i,env=phi$env)) # (x - env$min(-0.05))/(env$max(1.05) - env$min(-0.05)) - 0.5
    }
    if (nrk) {
      rk <- term[[label]]$rk
      
      # print(rk)
      # print("=======================")
      
       for (i in 1:nrk) {
        nq <- nq+1 # now nq=1
        r <- array(c(r,rk$fun(x,x.basis,nu=i,env=rk$env,out=TRUE)),c(nobs,nbasis,nq)) # dim=N*e*1
        # (1) scale x(ti) and x.basis(unique ti) by (x - env$min(-0.05))/(env$max(1.05) - env$min(-0.05))
        # (2) use the following function to evaluate rk
        # rk <- function(x, y) {
        #   k2 <- function(x) ((x - 0.5)^2 - 1/12)/2
        #   k4 <- function(x) ((x - 0.5)^4 - (x - 0.5)^2/2 + 7/240)/24
        #   k2(x) * k2(y) - k4(abs(x - y))
        # }
        # (3) outer(x,x.basis,rk), i.e., outer product using function rk defined above
        # (4) reshape result in (3) to be dim=N*e*1
      }
    }
  }
  if (is.null(r))
    stop("gss error in ssanova9: use lm for models with only unpenalized terms")
  ## Add the partial term
  if (!is.null(partial)) {
    mf.p <- model.frame(partial,data)
    for (lab in colnames(mf.p)) mf[,lab] <- mf.p[,lab]
    mt.p <- attr(mf.p,"terms")
    lab.p <- labels(mt.p)
    matx.p <- model.matrix(mt.p,data)[,-1,drop=FALSE]
    if (dim(matx.p)[1]!=dim(mf)[1])
      stop("gss error in ssanova9: partial data are of wrong size")
    matx.p <- scale(matx.p)
    center.p <- attr(matx.p,"scaled:center")
    scale.p <- attr(matx.p,"scaled:scale")
    s <- cbind(s,matx.p)
    part <- list(mt=mt.p,center=center.p,scale=scale.p)
  }
  else part <- lab.p <- NULL
  if (qr(s)$rank<dim(s)[2])
    stop("gss error in ssanova9: unpenalized terms are linearly dependent")
  ## Prepare the data
  y <- model.response(mf,"numeric")
  offset <- model.offset(mf)
  if (!is.null(offset)) {
    term$labels <- c(term$labels,"offset")
    term$offset <- list(nphi=0,nrk=0)
    y <- y - offset
  }
  
  
  #tm.inner.cal.gtermSR <-   proc.time() - tm.inner.cal.gtermSR
  
  #tm.inner.cal.gtermSR.total <<- tm.inner.cal.gtermSR.total + tm.inner.cal.gtermSR
  
  
  
  ## Fit the model
  if (nq==1) {
    r <- r[,,1]
    z <- sspreg91(s,r,r[id.basis,],y,cov,method,alpha,varht,tiidin,tiidcin)
  }
  else z <- mspreg91(s,r,id.basis,y,cov,method,alpha,varht,skip.iter)
  
  
  
  
  ## Brief description of model terms
  desc <- NULL
  for (label in term$labels)
    desc <- rbind(desc,as.numeric(c(term[[label]][c("nphi","nrk")])))
  if (!is.null(partial)) {
    desc <- rbind(desc,matrix(c(1,0),length(lab.p),2,byrow=TRUE))
  }
  desc <- rbind(desc,apply(desc,2,sum))
  if (is.null(partial)) rownames(desc) <- c(term$labels,"total")
  else rownames(desc) <- c(term$labels,lab.p,"total")
  colnames(desc) <- c("Unpenalized","Penalized")
  ## Return the results
  obj <- c(list(call=match.call(),mf=mf,terms=term,desc=desc,alpha=alpha,
                id.basis=id.basis,partial=part,lab.p=lab.p,cov=cov,
                skip.iter=skip.iter,r=r,s=s,q=r[id.basis,]),z)
  class(obj) <- c("ssanova9","ssanova")
  obj
}






## Fit Single Smoothing Parameter (Gaussian) REGression
sspreg91 <- function(s,r,q,y,cov,method,alpha,varht,tiidin,tiidcin)
{

  ## get dimensions
  nobs <- nrow(r)
  nxi <- ncol(r)
  if (!is.null(s)) {
    if (is.vector(s)) nnull <- 1
    else nnull <- ncol(s)
  }
  else nnull <- 0
  nn <- nxi + nnull
  
  
  
  ## cv function
  cv <- function(lambda) {
    q.wk <- 10^(lambda[1]+theta)*q
    if (length(lambda)-1) {
      stop("first if statement in cv is running, change the function!")
      ww <- cov$fun(lambda[-1],cov$env)
      ww <- chol(ww)
      y.wk <- forwardsolve(t(ww),y)
      s.wk <- forwardsolve(t(ww),s)
      r.wk <- forwardsolve(t(ww),r)
    }
    
    z <- .Fortran("reg",
                  as.double(cbind(s.wk,10^theta*r.wk)),
                  as.integer(nobs), as.integer(nnull),
                  as.double(q.wk), as.integer(nxi), as.double(y.wk),
                  as.integer(switch(method,"u"=1,"v"=2,"m"=3)),
                  as.double(alpha), varht=as.double(varht),
                  score=double(1), dc=double(nn),
                  as.double(.Machine$double.eps),
                  chol=double(nn*nn), double(nn),
                  jpvt=as.integer(c(rep(1,nnull),rep(0,nxi))),
                  wk=double(3*nobs+nnull), rkv=integer(1), info=integer(1),
                  PACKAGE="gss")[c("score","varht","dc","chol","jpvt","wk","rkv","info")]
    # print(z$jpvt)
    # print(z$rkv)
    if (z$info) stop("gss error in ssanova9: evaluation of GML score fails")
    if (!nnull|method%in%c("u","v")){
      
      if(class(ww)[1]=="matrix"){
        # ww now is the chol of original matrix ww
        detw <- 2*sum(log(diag(ww)))
      }else if(class(ww)[1]=="list"){
        # we need the chol of each element first and then sum the log of their diag
        detw <- 2 * sumLogD
        # check detw result
        # print("check detw result")
        # print(detw)
        # print(2*sum(log(diag(chol(WVW_all_g2)))))
      }
    } 
    else {
      
      # # qr(s) is fast since it only contains two columns (check by print out)
      # # we could re-essemble t(ww) here and see whether qr.qty could work
      # wk <- qr.qty(qr(s),t(ww))[-(1:nnull),]   # return: complete t(qr(s)) %*% t(ww), complete means complete=T in qr()
      # detw <- sum(log(eigen(wk%*%t(wk))$value))
      
      # cov matrix known in our case, detw not needed in the following calculation for score.
      detw <- 0
      
      
    }
    
    if (method=="m") score <- z$score*exp(detw/(nobs-nnull))
    if (method=="u") score <- z$wk[1]/varht+detw/nobs+2*alpha*z$wk[2]
    if (method=="v") score <- log(z$wk[1])+detw/nobs+2*alpha*z$wk[2]/(1-z$wk[2])
    alpha.wk <- max(0,log.la0-lambda[1]-5)*(3-alpha) + alpha
    alpha.wk <- min(alpha.wk,3)
    if (alpha.wk>alpha) {
      if (method=="u") score <- score + (alpha.wk-alpha)*2*z$wk[2]
      if (method=="v") score <- score + (alpha.wk-alpha)*2*z$wk[2]/(1-z$wk[2])
    }
    z$score <- score
    assign("fit",z[c(1:5,7)],inherits=TRUE)
    score
  }
  
  
  
  
  
  
  cv.wk <- function(lambda) cv.scale*cv(lambda)+cv.shift
  ## initialization
  #print(is.null(s))
  tmp <- sum(r^2)
  if (is.null(s)) theta <- 0
  else theta <- log10(sum(s^2)/nnull/tmp*nxi) / 2
  log.la0 <- log10(tmp/sum(diag(q))) + theta
  ## lambda search
  fit <- NULL
  la <- c(log.la0,cov$init)
  
  
  
  if (length(la)-1) {
    stop("first repeat: nlm is running, change the function!")
    counter <- 0
    ## scale and shift cv
    tmp <- abs(cv(la))
    cv.scale <- 1
    cv.shift <- 0
    if (tmp<1&tmp>10^(-4)) {
      cv.scale <- 10/tmp
      cv.shift <- 0
    }
    if (tmp<10^(-4)) {
      cv.scale <- 10^2
      cv.shift <- 10
    }
    repeat {
      zz <- nlm(cv.wk,la,stepmax=1,ndigit=7)
      if (zz$code<=3) break
      la <- zz$est
      counter <- counter + 1
      if (counter>=5) {
        warning("gss warning in ssanova9: iteration for model selection fails to converge")
        break
      }
    }
  }
  else {
    
    ww <- cov$fun(cov$env) # exactly the same as the list of cov we passed in

    countsub <- tiidin
    countsubs <- tiidcin
    
    
    if(class(ww)[1]=="matrix"){
      ww <- chol(ww)
      y.wk <- forwardsolve(t(ww),y)
      s.wk <- forwardsolve(t(ww),s)
      r.wk <- forwardsolve(t(ww),r)
    }else if(class(ww)[1]=="list"){

      
      cholfs <- cholfs.arma(countsubs,as.matrix(y),s,r,ww)
      sumLogD <- cholfs$sumLogD
      y.wk <- as.numeric(cholfs$y.wk)
      s.wk <- cholfs$s.wk
      r.wk <- cholfs$r.wk
      
    }
    
    
    
    mn0 <- log.la0-6
    mx0 <- log.la0+6
    
    # counting the number of iterations
    countingI = 1

    repeat {
      mn <- max(la-1,mn0)
      mx <- min(la+1,mx0)

      
      zz <- nlm0(cv,c(mn,mx))
      
      if ((min(zz$est-mn,mx-zz$est)>=1e-1)||
          (min(zz$est-mn0,mx0-zz$est)<1e-1)) break
      else la <- zz$est
      countingI = countingI + 1
      
    }

  }
  
  
  
  ## return
  lambda <- zz$est
  
  if(lambda < 0){
    lambda=0
  }
  
  jk1 <- cv(lambda)
  q.wk <- 10^(theta)*q
  if (length(lambda)-1) ww <- cov$fun(lambda[-1],cov$env)
  else ww <- cov$fun(cov$env)

  
  se.aux <- regaux(s.wk,10^theta*r.wk,q.wk,lambda[1],fit)
  c <- fit$dc[nnull+(1:nxi)]
  if (nnull) d <- fit$dc[1:nnull]
  else d <- NULL
  c(list(method=method,theta=theta,c=c,d=d,nlambda=lambda[1],zeta=lambda[-1],
         qout = q.wk,rout = 10^theta*r.wk,sout = s.wk),
    fit[-3],list(se.aux=se.aux))
}







## Auxiliary Quantities for Standard Error Calculation
regaux <- function(s,r,q,nlambda,fit)   #regaux(s.wk,10^theta*r.wk,q.wk,lambda[1],fit):  s.wk and r.wk are s and r pre-mult C^{-T};  q.wk=theta*Q
{
  nnull <- dim(s)[2]       # dimension of null space = length of vector d
  nn <- nnull +  dim(q)[1] # dimension of null space + dimension of approximated reproducing kernel space H1
  zzz <- eigen(q,symmetric=TRUE) 
  rkq <- min(fit$rkv-nnull,sum(zzz$val/zzz$val[1]>sqrt(.Machine$double.eps)))
  val <- zzz$val[1:rkq]             # first rkq eigenvalue of zzz -- element of D_Q page 90
  vec <- zzz$vec[,1:rkq,drop=FALSE] # first rkq eigenvector of zzz -- element of P_1 and P_2 page 90
  
  if (nnull) { # we run this
    wk1 <- qr(s) 
    wk1 <- (qr.qty(wk1,r%*%vec))[-(1:nnull),] 
  }
  else wk1 <- r%*%vec
  
  wk2 <- t(t(wk1)/sqrt(val))  
  wk2 <- t(wk2)%*%wk2         
  wk2 <- solve(wk2+diag(10^nlambda,dim(wk2)[1]),wk2)  # book page 89
  wk2 <- (wk2+t(wk2))/2
  wk2 <- t(wk2/sqrt(val))/sqrt(val)
  wk2 <- diag(1/val,dim(wk2)[1])-wk2
  z <- .Fortran("regaux",
                as.double(fit$chol), as.integer(nn),
                as.integer(fit$jpvt), as.integer(fit$rkv),
                drcr=as.double(t(cbind(s,r))%*%r%*%vec), as.integer(rkq),
                sms=double(nnull^2), as.integer(nnull), double(nn*nnull),
                PACKAGE="gss")[c("drcr","sms")]

  
  drcr <- matrix(z$drcr,nn,rkq)
  dr <- drcr[1:nnull,,drop=FALSE]
  sms <- 10^nlambda*matrix(z$sms,nnull,nnull)
  wk1 <- matrix(0,nnull+rkq,nnull+rkq)
  wk1[1:nnull,1:nnull] <- sms
  wk1[1:nnull,nnull+(1:rkq)] <- -t(t(dr)/val)
  wk1[nnull+(1:rkq),nnull+(1:rkq)] <- wk2
  suppressWarnings(z <- chol(wk1,pivot=TRUE))
  wk1 <- z
  rkw <- attr(z,"rank")
  while (wk1[rkw,rkw]<wk1[1,1]*sqrt(.Machine$double.eps)) rkw <- rkw-1
  wk1[row(wk1)>col(wk1)] <- 0
  if (rkw<nnull+rkq)
    wk1[(rkw+1):(nnull+rkq),(rkw+1):(nnull+rkq)] <- diag(0,nnull+rkq-rkw)
  hfac <- wk1
  hfac[,attr(z,"pivot")] <- wk1
  list(vec=vec,hfac=hfac)   # vec is first rkq eigen vector of matrix theta*Q
  # hfac is block matrix, left top: SMS, right top: d/eigenvalue(theta*Q)
  #                       left bot: 0  , right bot: ...
}






fitted.ssanova <- function(object,...)
{
  mf <- object$mf
  if (!is.null(object$random)) mf$random <- I(object$random$z)
  predict.ssanova(object,mf,estSigma2in = NULL)
}





predict.ssanova <- function(object,newdata,se.fit=FALSE,
                            include=c(object$terms$labels,object$lab.p),estSigma2in,...)
{


  nnew <- nrow(newdata)
  nbasis <- length(object$id.basis)
  nnull <- length(object$d) 
  nz <- length(object$b)
  nn <- nbasis + nnull + nz
  labels.p <- object$lab.p  

  term <- object$terms
  philist <- rklist <- NULL
  s <- r <- NULL
  nq <- 0

  
  
  for (label in include) {
    
    if (label=="1") {
      philist <- c(philist,term[[label]]$iphi)  # now: philist=1
      s <- cbind(s,rep(1,len=nnew))             # s= column of 1 with length nnew=nrow(newdata)
      next
    }
    
    if (label%in%labels.p) next # not run in our case
    if (label=="offset") next   # not run in our case
    
    xnew <- newdata[,term[[label]]$vlist] # data column: ti
    x <- object$mf[object$id.basis,term[[label]]$vlist] # mf is the dataframe.
    nphi <- term[[label]]$nphi  # 1
    nrk <- term[[label]]$nrk    # 1
    
    if (nphi) {
      iphi <- term[[label]]$iphi # 2
      phi <- term[[label]]$phi   # list of info. of phi for ti
      for (i in 1:nphi) {
        philist <- c(philist,iphi+(i-1)) # now: philist=c(1,2)
        s <- cbind(s,phi$fun(xnew,nu=i,env=phi$env))  # function is getting phi of var. ti 
        # s= column of 1 and column of phi of ti
      }
    }
    if (nrk) {
      irk <- term[[label]]$irk # 1
      rk <- term[[label]]$rk   # info. of reproducing kernel
      for (i in 1:nrk) {       
        rklist <- c(rklist,irk+(i-1)) # rklist=c(1)
        nq <- nq+1 # nq=1, number of Q matrix in the book
        r <- array(c(r,rk$fun(xnew,x,nu=i,env=rk$env,out=TRUE)),c(nnew,nbasis,nq)) # matrix R in the book
      }
    }
  }
  
  # below if, not run in our case
  if (!is.null(object$partial)) {
    vars.p <- as.character(attr(object$partial$mt,"variables"))[-1]
    facs.p <- attr(object$partial$mt,"factors")
    vlist <- vars.p[as.logical(apply(facs.p,1,sum))]
    for (lab in labels.p) {
      if (lab%in%include) {
        vlist.wk <- vars.p[as.logical(facs.p[,lab])]
        vlist <- vlist[!(vlist%in%vlist.wk)]
      }
    }
    if (length(vlist)) {
      for (lab in vlist) newdata[[lab]] <- 0
    }
    matx.p <- model.matrix(object$partial$mt,newdata)[,-1,drop=FALSE]
    matx.p <- sweep(matx.p,2,object$partial$center)
    matx.p <- sweep(matx.p,2,object$partial$scale,"/")
    nu <- nnull-dim(matx.p)[2]
    for (label in labels.p) {
      nu <- nu+1
      if (label%in%include) {
        philist <- c(philist,nu)
        s <- cbind(s,matx.p[,label])
      }
    }
  }
  
  
  r.wk <- matrix(0,nnew,nbasis)
  nq <- 0
  for (i in rklist) { # rklist=1
    nq <- nq + 1      # nq = 1
    r.wk <- r.wk + 10^object$theta[i]*r[,,nq]   # r.wk is now matrix R*theta in the book, why theta --> book p87
  }

  if (nz) {
    if (is.null(newdata$random)) z.wk <- matrix(0,nnew,nz)
    else z.wk <- newdata$random
    r.wk <- cbind(r.wk,z.wk)
  }
  ## Compute posterior mean
  nphi <- length(philist)  # philist=c(1,2); nphi=2
  pmean <- as.vector(r.wk%*%c(object$c,object$b))  # pmean = (theta * R %*% c)  since b is NULL
  if (nphi) pmean <- pmean + as.vector(s%*%object$d[philist]) # pmean = (theta * R %*% c) + S %*% d
  
  
  # not run below if in our case
  if (any(include=="offset")) {
    if (is.null(model.offset(object$mf)))
      stop("gss error: no offset in the fit")
    offset <- newdata$offset
    if (is.null(offset)) offset <- newdata$"(offset)"
    if (is.null(offset)) stop("gss error: missing offset")
    pmean <- pmean + offset
  }

  if (se.fit) {
    
    b <- estSigma2in/10^object$nlambda   # b in book page86
    ## Compute posterior variance
    ss <- matrix(0,nnull,nnew)  # nnull = 2
    if (!is.null(philist)) ss[philist,] <- t(s)  # matrix ss = S^T
    rr <- t(r.wk%*%object$se.aux$vec)            # matrix rr = (theta * R %*% vec)^T (vec is first rkq eigen vector of matrix theta*Q)
    wk <- object$se.aux$hfac%*%rbind(ss,rr)
    pse <- sqrt(b*apply(wk^2,2,sum))
    list(fit=pmean,se.fit=pse)
  }
  else
    pmean
}




