# Code related to the the paper "Functional static factor model" 
# Authors : G.Nisol, S.Tavakoli and M.Hallin.
# Responsability of the code : G.Nisol
# For the code related to the simulation, cf. file simu_FFM.R and FFM_all_functions.R
# For the code related to the real-data application, cf. "financial_app.R". It constitutes an example on how to use the FFM model
# In this code, we define a class "mits (for MIxed-nature Time Series, that is basically our panel dataset that can have scalar or functional units)

# Loading the required packages. 
library(methods)
library(LaplacesDemon)
require('fda')
require('pracma')
require('MASS')
require('MTS')
#require('PANICr')
require('xts')
require(gdata)
require('data.table')
require('doParallel')
require('foreach')
require('RSpectra')
require('tictoc')
require('VAR.etp')

#Definition of the class "MITS" that deals with mixed-nature FTS
# In this class, the scalar TS are stored in scalar while the functional ones (using package fda) in func. All fields and methods are described below
mits <- setRefClass('mits',
                    fields = list(
                      indexFunc="vector", #Index of every functional time series in the cross-section
                      indexScalar="vector",#Index of every scalar time series in the cross-section
                      indexType="vector", #Vector of boolean indicating whether indexes in the cross-section are functional or scalar
                      dimFunc="vector", #Dimension of every functional time series. Can be different
                      dimCum="vector", #Cumulative dimension : the i-th value of this vector is the dimension of the cross-section if we select the i first time series (functional or scalar)
                      nFunc="numeric", #Number of functional time series
                      nScalar="numeric", #Number of scalar time series
                      nTotal="numeric", #Number of time series (sum of nScalar and nFunc)
                      nDim="numeric", #Number of dimensions in the cross-section (last element of dimCum)
                      sampleSize="numeric", #=T
                      func="list", #List of FTS
                      scalar="list",  #List of scalar time series
                      coefs="matrix", #Cross-section of coefficients
                      inProd="matrix", #Cross-section of inner product between the cross-section and the basis functions
                      V="matrix"), #matrix of the inner product between the basis functions elements
                    methods= list(
                      initialize = function() {
                        indexFunc <<- vector()
                        indexScalar <<- vector()
                        indexType <<- vector()
                        dimFunc <<- vector()
                        dimCum <<- vector()
                        nFunc <<- 0
                        nScalar <<- 0
                        nTotal <<- 0
                        nDim <<- 0 
                        sampleSize <<- 0 
                        func <<- list()
                        scalar <<- list()
                        coefs <<- matrix()
                        inProd <<- matrix()
                        V <<- matrix()
                      },
                      loadFunc = function(funcList,position=(nTotal+1)) {
                        #This function adds "funcList" to the list of FTS stored in func. It adapts the indexes."position" refers to the position in the cross-section at which funcList should be added.
                        if(is.fd(funcList)) { funcList <- list(funcList)}
                        if(!is.list(funcList)) stop('The first argument is not a list nor a fd')
                        if(position<1 || position > nTotal +1) stop('The position that the user has specified is not in the correct range of pre-existing indexes')
                        size <- length(funcList)
                        saveType <- indexType[position:nTotal]
                        saveFunc <- indexFunc[position:nTotal]
                        saveScalar <- indexScalar[position:nTotal]
                        saveN <- nTotal+1
                        for(i in 1:size){
                          pushFunc(funcList[[i]])
                          indexType[position+i-1] <<- TRUE
                          indexFunc[position+i-1] <<- nFunc
                          indexScalar[position+i-1] <<- 0
                        }
                        if(position < saveN) {
                          indexType <<- c(indexType[1:(position+size-1)],saveType)
                          indexFunc <<- c(indexFunc[1:(position+size-1)],saveFunc)
                          indexScalar <<- c(indexScalar[1:(position+size-1)],saveScalar)
                        }
                        #computeV();computeCoefs(); #Rather use actualize()
                      },
                      pushFunc = function(funcA) {
                        #This function adds one single FTS to func. It is called by loadFunc(). It adaps the number of functions and dimensions as well.
                        if(!is.fd(funcA)) stop('The argument is not of class fd')
                        if(nTotal==0) sampleSize <<- dim(funcA$coefs)[2]
                        if (dim(funcA$coefs)[2] != sampleSize ) stop('The length of the functional time series introduced does not fit with previously introduced data')
                        nFunc <<- nFunc+1
                        nTotal <<- nTotal+1
                        nDim <<- nDim + dim(funcA$coefs)[1]
                        dimFunc[nFunc] <<- dim(funcA$coefs)[1]
                        func[[nFunc]] <<- funcA			
                      },
                      loadScalar = function(scalarList,position=(nTotal+1)) {
                        #Add the list (or matrix) of scalarList containing scalar time series to the class' fields scalar. "position" refers to the position in the cross-section at which scalarList should be added.
                        if(is.vector(scalarList)) { scalarList <- t(as.matrix(scalarList))}
                        if(!is.list(scalarList) && !is.matrix(scalarList)) stop('The first argument is not a matrix nor a ts')
                        if(position<1 || position > nTotal +1) stop('The position that the user has specified is not in the correct range of pre-existing indexes')
                        if(is.list(scalarList)) {
                          size <- length(scalarList)	
                        } else if(is.matrix(scalarList)) {
                          size <- dim(scalarList)[1]
                        }		
                        saveType <- indexType[position:nTotal]
                        saveScalar <- indexScalar[position:nTotal]
                        saveFunc <- indexFunc[position:nTotal]
                        saveN <- nTotal+1
                        for(i in 1:size){
                          if(is.list(scalarList)){
                            pushScalar(scalarList[[i]])
                          } else if(is.matrix(scalarList)) {
                            pushScalar(scalarList[i,])
                          }				
                          indexType[position+i-1] <<- FALSE
                          indexScalar[position+i-1] <<- nScalar
                          indexFunc[position+i-1] <<- 0
                        }
                        if(position < saveN) {
                          indexType <<- c(indexType[1:(position+size-1)],saveType)
                          indexScalar <<- c(indexScalar[1:(position+size-1)],saveScalar)
                          indexFunc <<- c(indexFunc[1:(position+size-1)],saveFunc)
                        }
                        #computeV();computeCoefs();
                      },
                      pushScalar = function(scalarA) {
                        #Add a single scalar time series to scalar. Adaps the dimension as well.
                        if(!is.vector(scalarA)) stop('The argument is not a vector')
                        if(nTotal==0) sampleSize <<- length(scalarA)
                        if (length(scalarA) != sampleSize ) stop('The length of the time series introduced does not fit with previously introduced data')
                        nScalar <<- nScalar+1
                        nTotal <<- nTotal+1
                        nDim <<- nDim + 1
                        scalar[[nScalar]] <<- scalarA			
                      },
                      computeV = function() {
                        #It computes the class' field V.
                        computeDimCum()
                        V <<- matrix(0,nDim,nDim)
                        for (i in 1:nTotal) {
                          if(!indexType[i]) {
                            V[(dimCum[i]+1),(dimCum[i]+1)] <<- 1
                          } else {
                            V[(dimCum[i]+1):(dimCum[i]+dimFunc[indexFunc[i]]),(dimCum[i]+1):(dimCum[i]+dimFunc[indexFunc[i]])] <<- inprod(func[[indexFunc[i]]]$basis,func[[indexFunc[i]]]$basis)
                          }
                        }
                      },
                      computeDimCum = function() {
                        # It computes the cumulative dimension based on indexFunc and indexType.
                        dimCum[1] <<- 0
                        for(i in 2:(nTotal+1)) {
                          if(!indexType[i-1]) {
                            dimCum[i] <<- dimCum[i-1]+1
                          } else {
                            dimCum[i] <<- dimCum[i-1]+ dimFunc[indexFunc[i-1]]
                          }
                        }
                        
                      },
                      computeCoefs = function() {
                        #It computes the coefficients of the mixed-nature panel
                        computeDimCum()
                        coefs <<- matrix(0, nDim, sampleSize)
                        for (j in 1:nTotal) {
                          if(!indexType[j]) {
                            coefs[(dimCum[j]+1),] <<- scalar[[indexScalar[j]]]
                          } else {
                            coefs[(dimCum[j]+1):(dimCum[j]+dimFunc[indexFunc[j]]),] <<- func[[indexFunc[j]]]$coefs
                          }
                        }
                      },
                      computeInProd = function() {
                        #It computes the inner product inProd
                        computeDimCum()
                        inProd <<- matrix(0, nDim, sampleSize)
                        for (j in 1:nTotal) {
                          if(!indexType[j]) {
                            inProd[(dimCum[j]+1),] <<- scalar[[indexScalar[j]]]
                          } else {
                            inProd[(dimCum[j]+1):(dimCum[j]+dimFunc[indexFunc[j]]),] <<- inprod(func[[indexFunc[j]]]$basis,func[[indexFunc[j]]])
                          }
                        }
                      },
                      actualize = function() {
                        #This function updates coefs, V, inprod and dimcum (through computeInProd())
                        computeV()
                        computeInProd()
                        computeCoefs()
                      },
                      trimSampleSize = function(N) {
                        #This function trims the sample size of the panel. It keeps only the first N time series
                        if(nFunc>= 1) {
                          for(i in 1:nFunc) {
                            func[[i]]$coefs <<- func[[i]]$coefs[,1:N]
                            
                          }
                        }
                        if(nScalar>=1) {
                          for(i in 1:nScalar) {
                            scalar[[i]] <<-scalar[[i]][1:N] 
                          }
                        }
                        sampleSize <<- N
                        coefs <<- coefs[,1:N]
                        inProd <<- inProd[,1:N]
                      },
                      shuffleAll =function() {
                        #randomly shuffles the indexes of the cross-section
                        perm <- sample.int(nFunc+nScalar)
                        indexAll <- pmax(indexScalar,indexFunc)
                        indexAllPerm <- indexAll[perm]
                        indexType <<- indexType[perm]
                        indexScalar <<- rep(0,nTotal)
                        indexScalar[indexType==0] <<- indexAllPerm[indexType==0]
                        indexFunc <<- rep(0,nTotal)
                        indexFunc[indexType==1] <<- indexAllPerm[indexType==1]
                        actualize()
                      },
                      reshapeFunc = function(index,trim) {
                        #Limit the function defined by its "index" to its first "trim" components
                        new_dim <- dimFunc[index]-trim
                        if(dimFunc[index]<=trim){
                          stop('number of dimension to remove smaller than number of basis functions of the function of the MTS to reshape')
                        }
                        else {
                          func[[index]]$basis$nbasis <<- new_dim
                          func[[index]]$basis$names <<- func[[1]]$basis$names[1:new_dim]
                          tmp <- func[[index]]$coefs[1:new_dim,]
                          func[[index]]$coefs <<- tmp
                        }
                      },
                      cut = function(n) {
                        #Remove the last (nTotal-n) components of the cross-section
                        indexType <<- indexType[1:n]
                        f <- rep(FALSE,nTotal-n)
                        nTotal <<- n
                        keep_func <- indexFunc[c(indexType,f)]
                        nFunc <<- length(keep_func)
                        indexFunc <<- indexFunc[1:nFunc]
                        func <<- func[sort(keep_func)]
                        indexFunc <<- rank(indexFunc)
                        dimFunc <<-dimFunc[keep_func]
                        if(nScalar>0){
                          keep_scalar <- indexScalar[c(!indexType,f)]
                          nScalar <<- length(keep_scalar)
                          indexScalar <<-indexScalar[1:nScalar]
                          scalar <<-scalar[sort(keep_scalar)]
                          indexScalar <<- rank(indexScalar)
                        }
                        
                        
                        dimCum <<- dimCum[1:(nFunc+1)]
                        nDim <<- dimCum[nFunc+1]
                        
                        coefs <<- coefs[1:nDim,]
                        inProd <<- inProd[1:nDim,]
                        V <<- V[1:nDim,1:nDim]
                      },
                      removeMean = function() {
                        #Standardize every time series by substracting the mean
                        mu <- NULL
                        if(nScalar>=1) {
                          for(i in 1:nScalar) {
                            mum <- mean(scalar[[i]])
                            mu <- c(mu,mum)
                            scalar[[i]] <<- scalar[[i]]-mum
                          }
                        }
                        muf <- mu
                        if(nFunc>= 1) {
                          for(i in 1:nFunc) {
                            mum <- rowMeans(func[[i]]$coefs)
                            mu <- c(mu,mum)
                            func[[i]]$coefs <<- func[[i]]$coefs - mum
                            muf <- c(muf,inprod(func[[i]],func[[i]]$basis))
                          }
                        }
                        #The next two lines assume that the cross-sectional order is such that the scalar have been loaded first then the functions. If it is not the case, the next two lines must be commented and replaced by "actualize()".
                        coefs <<- coefs-mu
                        inProd <<- inProd-muf
                        return(mu)
                      }
                    )
)

is.single_numeric <- function(x){
  # This function returns TRUE when p is a single numeric value
  # Input :
  #x : the value to be checked
  # Output : boolean
  if(is.numeric(x) & length(x)==1){
    return(TRUE)
  } else{
    return(FALSE)
  }
}
diagg <- function(a){
  if(length(a)==1){
    return(a)
  } else{
    return(diag(a))
  }
}
est_u <- function(ts,r){
  #Estimate the r first factors. Only for zero-mean process  
  #Input : -ts (mits class) : the mixed-nature panel
  #       - r (integer) : the number of factors to be estimated
  #Output : the multiple time series of factors
  Omega <- crossprod(ts$coefs,crossprod(t(ts$V),ts$coefs))
  T <- dim(Omega)[1]
  N <- dim(ts$coefs)[1]
  ei <- eigs(Omega,r)
  return(t(tcrossprod(ei$vectors,diagg(sqrt(T)*rep(1,r)))))
}

est_loadings <- function(ts,u) {
  #Estimate the factor loadings . Only for zero-mean process
  #Input : -ts (mits class) : the mixed-nature panel
  #-u (matrix) : the factors
  #Output : matrix consisting of the representation of the factor loadings in the basis specified by the ad hoc argument in ts
  T <- dim(ts$coefs)[2]
  return(ts$coefs%*%t(u)/T)
}

est_common <- function(ts,r){
  u_hat <- est_u(ts,r)
  B_hat <- est_loadings(ts,u_hat)
  return(B_hat%*%u_hat)
}

normVec <- function (x) {
  # Next function : Euclidian norm.
  # Args:
  #   x : vector
  # Returns:
  #  The Euclian norm of the vector x 
  return(sqrt(sum(x^2)))
} 


est_r <- function(ts,r_max,crit='IC1cor'){
  SSE <- rep(NA,r_max)
  IC <- rep(NA,r_max)
  
  T <- ts$sampleSize
  N <- ts$nTotal
  if(crit=='IC1cor'){
    penalty <- sqrt((N+T)/(N*T))*log((N*T)/(N+T))
  } else if(crit=='IC2cor') {
    penalty <- sqrt((N+T)/(N*T))*log(min(N,T))
  }
  
  for(r in 1:r_max){
    SSE[r] <- normVec(ts$coefs-est_common(ts,r))
    IC[r] <- log((SSE[r]^2)/(N*T))+r*penalty
  }
  r_min <- which(IC==min(IC))
  return(r_min)
}

est_r_abc <- function(ts,rmax,crit='IC2cor',N_perm=5,c=seq(0,10,by=0.05)) {
  # This function computes the estimate number of static factor r and the optimal dimension of the projection p_i of the functional time series ts based on our adaptation of the ABC method
  # Input :
  # ts (class : mits) : the MITS
  # rmax (single numeric) : we are looking for the optimal r in the range [1,rmax]
  # pmax (single numeric) : we are looking for the optimal p_i in the range [1,pmax]
  # c (vector, length: n_c) : we look for the optimal scale of the penalization function among the values contained in the vector c
  # n (vector, length: n_n) : the values of the cross-section dimension for which we will compute r_opt and p_opt (permitting to estimate the variance of these estimate with respect to c). The last value of this vector should be n
  # T (vector, length: n_n) : the values of the sample size for which we will compute r_opt and p_opt. The last value of this vector should be the sample size of ts
  # crit (either "FBIC" or "FIC") : the penalization function to use
  #method : 'VAR', 'HHNT' or 'red'
  #Output : a couple containing
  # the estimated r
  # the optimal p
  if(class(ts)!="mits"){
    stop('ts is not of class mits')
  }
  if(!is.single_numeric(rmax)){
    stop('rmax is not a single numeric')
  }
  if(!is.vector(c)){
    stop('c is not a vector')
  }
  if(crit!='IC1cor' & crit!='IC2cor'){
    stop('crit does not take a recognized value')
  }
  
  L_c <- length(c)
  T <- ts$sampleSize
  N <- ts$nTotal
  
  L <- 10
  
  
  T_all <- rep(T,L)
  N_all <- seq(round(0.8*N),N,length.out=L)
  
  
  criterion <- array(Inf,dim=c(L_c,L,N_perm,rmax))
  SSE <- array(0,dim=c(L,N_perm,rmax))
  penalty <- matrix(0,L,rmax)
  
  
  for(L_i in 1:L) {
    N_i <- N_all[L_i]
    T_i <- T_all[L_i]
    for(r_i in 1:rmax){
      if(crit=='IC1cor'){
        penalty[L_i,r_i] <- r_i*sqrt((N_i+T_i)/(N_i*T_i))*log((N_i*T_i)/(N_i+T_i))
      } else if(crit=='IC2cor') {
        penalty[L_i,r_i] <- r_i*sqrt((N_i+T_i)/(N_i*T_i))*log(min(N_i,T_i))
      }
    }
    
    for(perm_i in 1:N_perm){
      
      ts2 <- ts$copy()
      
      ts2$trimSampleSize(T_i)
      if(perm_i>1){
        ts2$shuffleAll()
      }
      
      ts2$cut(N_i)
      ts2$removeMean()
      
      
      
      for(r_i in 1:rmax){
        SSE[L_i,perm_i,r_i] <-  log((normVec(ts2$coefs-est_common(ts2,r_i))^2)/(N_i*T_i))
      }
      
    }
    
  }
  S_p <- rep(0,L_c)
  r_opt <- array(0,dim=c(L_c,L,N_perm))
  
  for(c_i in 1:L_c) {
    for(L_i in 1:L) {
      for(perm_i in 1:N_perm){
        criterion[c_i,L_i,perm_i,] <- SSE[L_i,perm_i,]+c[c_i]*penalty[L_i,]
        r_opt[c_i,L_i,perm_i]  <- which(criterion[c_i,L_i,perm_i,]==min(criterion[c_i,L_i,perm_i,]))[1]
      }
      
    }
    S_p[c_i] <- var(as.vector(r_opt[c_i,,]))
  }
  
  
  idx <- which(S_p>0)[1]
  if(is.na(idx)) {
    warning('idx second part has NA')
    idx <- floor(L_c/2)
  }
  c_min <- idx+which(S_p[(idx+1):L_c]==0)[1]
  if(is.na(c_min)) {
    warning('c_min has NA')
    c_min <- floor(L_c/2)
  }
  r_p <- median(r_opt[c_min,L,])
  return(r_p)
}


est_r_abc2 <- function(ts,rmax,crit='IC2cor',N_perm=5,c=seq(0,10,by=0.05)) {
  # This function computes the estimate number of static factor r and the optimal dimension of the projection p_i of the functional time series ts based on our adaptation of the ABC method
  # Input :
  # ts (class : mits) : the MITS
  # rmax (single numeric) : we are looking for the optimal r in the range [1,rmax]
  # pmax (single numeric) : we are looking for the optimal p_i in the range [1,pmax]
  # c (vector, length: n_c) : we look for the optimal scale of the penalization function among the values contained in the vector c
  # n (vector, length: n_n) : the values of the cross-section dimension for which we will compute r_opt and p_opt (permitting to estimate the variance of these estimate with respect to c). The last value of this vector should be n
  # T (vector, length: n_n) : the values of the sample size for which we will compute r_opt and p_opt. The last value of this vector should be the sample size of ts
  # crit (either "FBIC" or "FIC") : the penalization function to use
  #method : 'VAR', 'HHNT' or 'red'
  #Output : a couple containing
  # the estimated r
  # the optimal p
  if(class(ts)!="mits"){
    stop('ts is not of class mits')
  }
  if(!is.single_numeric(rmax)){
    stop('rmax is not a single numeric')
  }
  if(!is.vector(c)){
    stop('c is not a vector')
  }
  if(crit!='IC1cor' & crit!='IC2cor'){
    stop('crit does not take a recognized value')
  }
  
  L_c <- length(c)
  T <- ts$sampleSize
  N <- ts$nTotal
  
  L <- 10
  
  
  T_all <- rep(T,L)
  N_all <- seq(round(0.8*N),N,length.out=L)
  
  
  criterion <- array(Inf,dim=c(L_c,L,N_perm,rmax))
  SSE <- array(0,dim=c(L,N_perm,rmax))
  penalty <- matrix(0,L,rmax)
  
  
  for(L_i in 1:L) {
    N_i <- N_all[L_i]
    T_i <- T_all[L_i]
    for(r_i in 1:rmax){
      if(crit=='IC1cor'){
        penalty[L_i,r_i] <- r_i*sqrt((N_i+T_i)/(N_i*T_i))*log((N_i*T_i)/(N_i+T_i))
      } else if(crit=='IC2cor') {
        penalty[L_i,r_i] <- r_i*sqrt((N_i+T_i)/(N_i*T_i))*log(min(N_i,T_i))
      }
    }
    
    for(perm_i in 1:N_perm){
      
      ts2 <- ts$copy()
      
      ts2$trimSampleSize(T_i)
      if(perm_i>1){
        ts2$shuffleAll()
      }
      
      ts2$cut(N_i)
      
      
      
      for(r_i in 1:rmax){
        SSE[L_i,perm_i,r_i] <-  log((normVec(ts2$coefs-est_common(ts2,r_i))^2)/(N_i*T_i))
      }
      
    }
    
  }
  S_p <- matrix(0,L_c,N_perm)
  r_opt <- array(0,dim=c(L_c,L,N_perm))
  
  for(c_i in 1:L_c) {
    for(L_i in 1:L) {
      for(perm_i in 1:N_perm){
        criterion[c_i,L_i,perm_i,] <- SSE[L_i,perm_i,]+c[c_i]*penalty[L_i,]
        r_opt[c_i,L_i,perm_i]  <- which(criterion[c_i,L_i,perm_i,]==min(criterion[c_i,L_i,perm_i,]))[1]
      }
      
    }
    for(perm_i in 1:N_perm){
      S_p[c_i,perm_i] <- var(as.vector(r_opt[c_i,,perm_i]))
    }
  }
  
  r_p <- rep(0,N_perm)
  idx <- rep(0,N_perm)
  c_min <- rep(0,N_perm)
  for(perm_i in 1:N_perm){
    idx[perm_i] <- which(S_p[,perm_i]>0)[1]
    if(is.na(idx[perm_i])) {
      warning('idx second part has NA')
      idx[perm_i] <- floor(L_c/2)
    }
    c_min[perm_i] <- idx[perm_i]+which(S_p[(idx[perm_i]+1):L_c,perm_i]==0)[1]
    if(is.na(c_min[perm_i])) {
      warning('c_min has NA')
      c_min[perm_i] <- floor(L_c/2)
    }
    r_p[perm_i] <- median(r_opt[c_min,L,perm_i])
  }
  return(median(r_p))
}




predict_NTH_rknown <- function(ts,h,p,basis,args,r){
  
  T <- ts$sampleSize
  u <- est_u(ts,r)
  Y <- t(ts$coefs[,(1+h):T])
  u_reg <- t(u[,1:(T-h),drop=FALSE])
  reg <- lm(Y~u_reg)
  coef_predicted <- as.vector(reg$coefficients[1,]+t(u[,T])%*%reg$coefficients[2:(r+1),])
  pred_points <- matrix(0,ts$nTotal,length(args))
  for(i in 1:ts$nTotal){
    pred_fd <- fd(coef_predicted[((i-1)*p+1):(i*p)],basis)
    for(j in 1:length(args)){
      pred_points[i,j] <- eval.fd(pred_fd,args[j])
    }
  }
  return(pred_points)
}

predict_NTH_rknown2 <- function(ts_entry,h,p,basis,args,r){
  ts <- ts_entry$copy()
  T <- ts$sampleSize
  mu <- ts$removeMean()
  u <- est_u(ts,r)
  Y <- t(ts$coefs[,(1+h):T])
  u_reg <- t(u[,1:(T-h),drop=FALSE])
  reg <- lm(Y~u_reg)
  coef_predicted <- mu + as.vector(reg$coefficients[1,]+t(u[,T])%*%reg$coefficients[2:(r+1),])
  pred_points <- matrix(0,ts$nTotal,length(args))
  for(i in 1:ts$nTotal){
    pred_fd <- fd(coef_predicted[((i-1)*p+1):(i*p)],basis)
    for(j in 1:length(args)){
      pred_points[i,j] <- eval.fd(pred_fd,args[j])
    }
  }
  return(pred_points)
}

predict_NTH_rknown_v2 <- function(ts_entry,h,p,basis,args,r){
  ts <- ts_entry$copy()
  T <- ts$sampleSize
  N <- ts$nDim
  
  mu <- ts$removeMean()
  u <- est_u(ts,r)
  B <- est_loadings(ts,u)
  chi <-  B%*% u
  xi <- ts$coefs - chi
  u_pred <- rep(0,r)
  mod <- list()
  for (i in 1:r) {
    #mod[[i]] <- auto.arima(u[i,])
    mod[[i]] <- auto.arima(u[i,],  ic= "bic")
    u_pred[i] <-  forecast(mod[[i]], h)$mean[h]
  }
  chi_pred <- B %*% u_pred
  xi_pred  <- rep(0,N)
  
  mod <- list()
  for(i in 1:N){
    mod[[i]] <- auto.arima(xi[i,], ic= "bic")
    xi_pred[i] <-  forecast(mod[[i]], h)$mean[h]
  }
  coef_predicted <- mu + chi_pred + xi_pred
  pred_points <- matrix(0,ts$nTotal,length(args))
  for(i in 1:ts$nTotal){
    pred_fd <- fd(coef_predicted[((i-1)*p+1):(i*p)],basis)
    for(j in 1:length(args)){
      pred_points[i,j] <- eval.fd(pred_fd,args[j])
    }
  }
  return(pred_points)
}


predict_NTH_runknown <- function(ts,h,p,basis,args){
  r <- est_r_abc(ts,10)
  print('estimated r')
  print(r)
  T <- ts$sampleSize
  u <- est_u(ts,r)
  Y <- t(ts$coefs[,(1+h):T])
  u_reg <- t(u[,1:(T-h),drop=FALSE])
  reg <- lm(Y~u_reg)
  coef_predicted <- as.vector(reg$coefficients[1,]+t(u[,T])%*%reg$coefficients[2:(r+1),])
  pred_points <- matrix(0,ts$nTotal,length(args))
  for(i in 1:ts$nTotal){
    pred_fd <- fd(coef_predicted[((i-1)*p+1):(i*p)],basis)
    for(j in 1:length(args)){
      pred_points[i,j] <- eval.fd(pred_fd,args[j])
    }
  }
  return(pred_points)
}

compute_forecast_error <- function(real,estimated,h,type="MSE",expo=TRUE){
  N <- dim(estimated[[h]])[1]
  T <- dim(estimated[[h]])[2]
  Age <- dim(estimated[[h]])[3]
  error <- matrix(0,N,T)
  for(n in 1:N){
    for(t in 1:T){
      for(age in 1:Age){
        if(type=="MSE"){
          if(expo){
            error[n,t] <- error[n,t]+(exp(estimated[[h]][n,t,age])-exp(real[[n]][age,25+h+t]))^2
          } else{
            error[n,t] <- error[n,t]+(estimated[[h]][n,t,age]-real[[n]][age,25+h+t])^2
          }
        }
        else if(type=="MASE"){
          if(expo){
            error[n,t] <- error[n,t]+abs(exp(estimated[[h]][n,t,age])-exp(real[[n]][age,25+h+t]))
          } else{
            error[n,t] <- error[n,t]+abs(estimated[[h]][n,t,age]-real[[n]][age,25+h+t])
          }
        }
      }
    }
  }
  return(mean(error)/Age)#avant : on ne divisé pas par Age
}

compute_norm <- function(real,h,type="MSE"){
  N <- length(real)
  T <- 42-25-h
  Age <- dim(real[[1]])[1]
  error <- matrix(0,N,T)
  for(n in 1:N){
    for(t in 1:T){
      for(age in 1:Age){
        if(type=="MSE"){
          error[n,t] <- error[n,t]+(real[[n]][age,25+h+t])^2
        }
        else if(type=="MASE"){
          error[n,t] <- error[n,t]+abs(real[[n]][age,25+h+t])
        }
      }
    }
  }
  return(mean(error)/Age) #avant : on ne divisé pas par Age
}

path_male = '/home/gilles/Documents/ffm/codes/mortality data/results/mortality_Male_v1_Gaodata_nozeroref.Rdata'
path_female = '/home/gilles/Documents/ffm/codes/mortality data/results/mortality_Female_v1_Gaodata_nozeroref.Rdata'
make_Table_forecast <- function (path_male, path_female){
  require('xtable')
  load(path_male)
  err_GAO_male <- error_GAO_log
  err_NTH_male <- error_NTH_log
  err_UNI_male <- error_UNI_log
  load(path_female)
  err_GAO_female <- error_GAO_log
  err_NTH_female <- error_NTH_log
  err_UNI_female <- error_UNI_log
  table <- 1000*cbind(err_GAO_female[,2],err_UNI_female[,2],err_NTH_female[,2],err_GAO_female[,1],err_UNI_female[,1],err_NTH_female[,1],err_GAO_male[,2],err_UNI_male[,2],err_NTH_male[,2],err_GAO_male[,1],err_UNI_male[,1],err_NTH_male[,1])
  m = colMeans(table)
  med <- apply(table,2,median)
  table <- rbind(table,m,med)
  xtable(table, digits=3)
}

make_Table_forecast_nolog <- function (path_male, path_female){
  require('xtable')
  load(path_male)
  err_GAO_male <- error_GAO
  err_NTH_male <- error_NTH
  err_UNI_male <- error_UNI
  load(path_female)
  err_GAO_female <- error_GAO
  err_NTH_female <- error_NTH
  err_UNI_female <- error_UNI
  table <- cbind(1000*err_GAO_female[,2],1000*err_UNI_female[,2],1000*err_NTH_female[,2],10000*err_GAO_female[,1],10000*err_UNI_female[,1],10000*err_NTH_female[,1],1000*err_GAO_male[,2],1000*err_UNI_male[,2],1000*err_NTH_male[,2],10000*err_GAO_male[,1],10000*err_UNI_male[,1],10000*err_NTH_male[,1])
  m = colMeans(table)
  med <- apply(table,2,median)
  table <- rbind(table,m,med)
  xtable(table, digits=2)
}