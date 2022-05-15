

### Codes to implement the methods in the paper 
### Paul, S., Nath, S., Warren, K., (2022). Causal Network Influence with Latent Homophily and Measurement Error: An Application to Therapeutic Community. 
### Preprint at arXiv: 2203.14223

######### Homophily correction with SAR model #########

### Functions ###



#library(MASS)

library(igraph)


### Utility functions

### Compute Derivatives and Hessians of the objective function for bias corrected and uncorrected model

funcder<-function(rho,YL,Y,lambdas,M){
    n=length(Y)
    s2 = t(Y)%*%M%*%Y - 2*rho*t(Y)%*%M%*%YL + rho^2* t(YL)%*%M%*%YL
derval = (2/n)*sum(lambdas/(1-rho*lambdas)) + 2*(rho*t(YL)%*%M%*%YL-t(Y)%*%M%*%YL)/s2
return(derval)
}

funchess <-function(rho,YL,Y,lambdas,M){
    n=length(Y)
    s2 = t(Y)%*%M%*%Y - 2*rho*t(Y)%*%M%*%YL + rho^2* t(YL)%*%M%*%YL
    hessval = (2/n)*sum((lambdas/(1-rho*lambdas))^2) + 2*t(YL)%*%M%*%YL/s2 - 4*(rho*t(YL)%*%M%*%YL-t(Y)%*%M%*%YL)^2/s2^2
    return(hessval)
}



funcderc<-function(rho,YL,Y,lambdas,M){
    n=length(Y)
    s2 = t(Y)%*%(M)%*%Y - 2*rho*t(Y)%*%(M)%*%YL + rho^2* t(YL)%*%(M)%*%YL
derval = (2/n)*sum(lambdas/(1-rho*lambdas)) + 2*(rho*t(YL)%*%(M)%*%YL-t(Y)%*%(M)%*%YL)/s2
return(derval)
}

funchessc <-function(rho,YL,Y,lambdas,M){
    n=length(Y)
   s2 = t(Y)%*%(M)%*%Y - 2*rho*t(Y)%*%(M)%*%YL + rho^2* t(YL)%*%(M)%*%YL
    hessval = (2/n)*sum((lambdas/(1-rho*lambdas))^2) + 2*t(YL)%*%(M)%*%YL/s2 - 4*(rho*t(YL)%*%(M)%*%YL-t(Y)%*%(M)%*%YL)^2/s2^2
    return(hessval)
}


# A simple matrix trace function
tr<-function(X){
  return(sum(diag(X)))
}




### No latent variable, only predictor Z

## Y is univariate response
## X is binary adjacency matrix (undirected)
## L is the Laplacian or the row-normalized adjacency (can be weighted) matrix
## d is the number of dimensions - can either be given or automatically estimated
## Z is the set of additional predictors.


sarMLclassicZ<-function(Y,L,X,Z){
n=length(Y)
d2 = dim(as.matrix(Z))[2]
lambdas = eigen(L)$values
YL = L%*%Y
nhat=Z
M = diag(n) - nhat%*%solve(t(nhat)%*%(nhat))%*%t(nhat)
rholast = c((t(Y)%*%M%*%YL)/(t(Y)%*%M%*%Y))
rhodiff =1
t=0

### optimize for rho
while(rhodiff>0.001 & t<50){
    rho = rholast - c(funcder(rholast,YL,Y,lambdas,M)/funchess(rholast,YL,Y,lambdas,M))
    rhodiff = abs(rho-rholast)
    rholast = rho
    t=t+1
}

### compute estimates for beta and sigmasq
hatZ = (diag(n) - rho*L)%*%Y
beta = c(solve(t(nhat)%*%nhat)%*%t(nhat)%*%hatZ)
sigmasqhat = (t(Y)%*%M%*%Y - 2*rho*t(Y)%*%M%*%YL + rho^2* t(YL)%*%M%*%YL)/n

   G= L%*%solve(diag(n)-rho*L)
   H= G%*%(nhat%*%beta)

## compute the Fisher information matrix

 Infmatrix = matrix(NA,d2+2,d2+2)
  Infmatrix[1:d2,] = cbind(t(nhat)%*%nhat,t(nhat)%*%H,rep(0,d2))
  Infmatrix[d2+1,] = cbind(t(H)%*%nhat, t(H)%*%H + sigmasqhat*tr((G+t(G))%*%G),tr(G))
  Infmatrix[d2+2,] = cbind(t(rep(0,d2)), tr(G), n/(2*sigmasqhat))

  Infmatrix=Infmatrix/c(n*sigmasqhat) 
  
  ## compute the standard errors

  cov = solve(Infmatrix)
  serho = sqrt(diag(cov)[d2+1]/n)
  sebeta = sqrt(diag(cov)[1:d2]/n)

prho = 2*(1-pnorm(abs(rho)/serho))
pbeta = 2*(1-pnorm(abs(beta)/sebeta))

### Return as list the estimates, SEs, and p-values of influence (rho) and coefficeint of Z (beta), along with estimate of variance of error term

return(list(influence = rho, SEinfluence =serho, pvaluerho=prho, betacovariate=beta, SEbetacovariate=sebeta,  pvaluecovariate = pbeta, sigmasqhat=sigmasqhat))
}




### Latent homphily correction but no bias correction

## Y is univariate response
## X is binary adjacency matrix (undirected)
## L is the Laplacian or the row-normalized adjacency (can be weighted) matrix
## d is the number of dimensions - currently needs to be given but can also be  automatically estimated using "ase" package
## Z is the set of additional predictors.
## Uhat is a precomputed estimate of latent homophily vector. It is an optional argument and typically will not be provided by the user. 
## Latent homophily is extracted from spectral decomposition of the adjacency matrix


sarMLZ<-function(Y,L,X,Z,d,Uhat){

n=length(Y)
d2 = dim(as.matrix(Z))[2]


### If latent factors Uhat are already precomputed, then skip this step: relevant for out of sample prediction task

## Extract latent homophily vectors through a spectral decomposition

if(missing(Uhat)){
spectra<-svd(X)
Uhat<-spectra$u[,1:d]%*%diag(sqrt(spectra$d[1:d]))
}
nhat<-cbind(Uhat,Z)

lambdas = eigen(L)$values
YL = L%*%Y
M = diag(n) - nhat%*%solve(t(nhat)%*%(nhat))%*%t(nhat)
rholast = c((t(Y)%*%M%*%YL)/(t(Y)%*%M%*%Y))
rhodiff =1
t=0

### Estimate rho 
while(rhodiff>0.001 & t<50){
    rho = rholast - c(funcder(rholast,YL,Y,lambdas,M)/funchess(rholast,YL,Y,lambdas,M))
    rhodiff = abs(rho-rholast)
    rholast = rho
    t=t+1
}

### Estimate beta and gamma

hatZ = (diag(n) - rho*L)%*%Y
beta = solve(t(nhat)%*%nhat)%*%t(nhat)%*%hatZ
sigmasqhat = (t(Y)%*%M%*%Y - 2*rho*t(Y)%*%M%*%YL + rho^2* t(YL)%*%M%*%YL)/n

   G= L%*%solve(diag(n)-rho*L)
   H= G%*%(nhat%*%beta)

### Estimate information matrix

 Infmatrix = matrix(NA,d+d2+2,d+d2+2)
  Infmatrix[1:(d+d2),] = cbind(t(nhat)%*%nhat,t(nhat)%*%H,rep(0,d+d2))
  Infmatrix[d+d2+1,] = cbind(t(H)%*%nhat, t(H)%*%H + sigmasqhat*tr((G+t(G))%*%G),tr(G))
  Infmatrix[d+d2+2,] = cbind(t(rep(0,d+d2)), tr(G), n/(2*sigmasqhat))

  Infmatrix=Infmatrix/c(n*sigmasqhat) 

  ### compute standard errors

  cov = solve(Infmatrix)
  serho = sqrt(diag(cov)[d+d2+1]/n)
  sebeta = sqrt(diag(cov)[1:(d)]/n)
  segamma = sqrt(diag(cov)[(d+1):(d+d2)]/n)

  #### separate beta and gamma parameters

betalatent=beta[1:d]
gamma=beta[(d+1):(d+d2)]

## compute p values

prho = 2*(1-pnorm(abs(rho)/serho))
pbeta = 2*(1-pnorm(abs(betalatent)/sebeta))
pgamma = 2*(1-pnorm(abs(gamma)/segamma))


### Return as list the estimates, SEs, and p-values of influence (rho), coefficeints of Z (betacovariate), and coefficient of 
### latent variables (betalatent, not very meaningful) along with estimate of variance of error term


return(list(influence = rho, SEinfluence =serho, pvaluerho=prho, betalatent =betalatent, SEbetalatent = sebeta, 
  pvaluelatent = pbeta, betacovariate=gamma,   SEbetacovariate=segamma, pvaluecovariate = pgamma, sigmasqhat=sigmasqhat,d=d))
}



#### latent homophily with bias correction + additional covariate

### Returns both parameter estimates and standard errors.



## Y is univariate response
## X is binary adjacency matrix (undirected)
## L is the Laplacian or the row-normalized adjacency (can be weighted) matrix
## d is the number of dimensions - currently needs to be given but can also be  automatically estimated using "ase" package
## Z is the set of additional predictors.
## Uhat is a precomputed estimate of latent homophily vector. It is an optional argument and typically will not be provided by the user. 
## comms is an optional parameter for number of communities in the network SBM model. If not provided then it defaults to d.

### Function fits a RDPG SBM model to network to extract latent variables and estimate covariance of the estimated latent variables. 
## The network effect is estimated controlling for those latent variables


sarMLZc<-function(Y,L,X,Z,d,comms,Uhat){
  n=length(Y)
  d2 = dim(as.matrix(Z))[2]

  ### If number of communities are not given then default is it is equal to d.
  if(missing(comms)){
  comms=d
}
#}

  if(missing(Uhat)){
    spectra<-svd(X)
  Uhat<-spectra$u[,1:d]%*%diag(sqrt(spectra$d[1:d]))
}
  #lambdas = svd(L)$d
  lambdas = eigen(L)$values
  YL = L%*%Y

  ### Perform k means clustering on Uhat and estimate cluster centers and sizes
  clus = kmeans(Uhat, centers=comms, nstart=5)
  cents = clus$centers
  pihat = clus$size/n 
  

  ### Compute the covariance matrix of the estimated Uhat

  Delhat = Reduce("+",lapply(1:comms, function(m){
    return(pihat[m]*cents[m,]%*%t(cents[m,]))
  }))
  
  deltal = lapply(1:comms, function(l){
    expectl =  lapply(1:comms,function(m){
      expU = cents[m,]%*%t(cents[m,])
      prodU = t(cents[l,])%*%cents[m,]
      deltaq = c(prodU - prodU^2)*expU
      return(pihat[m]*deltaq)
    })
    explhat = Reduce("+",expectl)
    return(solve(Delhat)%*% explhat%*%solve(Delhat))
  })
  
  deltahat = Reduce("+",lapply(1:comms, function(l){
    return(deltal[[l]]*c(pihat[l]))
  }))
  

  ### Compute M
  
  MZW = t(cbind(Uhat,Z))%*%(cbind(Uhat,Z))
  omega  = matrix(0,(d+d2),(d+d2))
  omega[1:d,1:d] = deltahat
  MW = t(Uhat)%*%(Uhat) - deltahat
  M = diag(n) - (cbind(Uhat,Z))%*%solve(MZW - omega)%*%t(cbind(Uhat,Z))
 
  
  rholast = c((t(Y)%*%(M)%*%YL)/(t(Y)%*%(M)%*%Y))
  rhodiff =1
  t=0

  ## Estimate rho


  while(rhodiff>0.001 & t<50){
    rho = rholast - c(funcderc(rholast,YL,Y,lambdas,M)/funchessc(rholast,YL,Y,lambdas,M))
    rhodiff = abs(rho-rholast)
    rholast = rho
    t=t+1
  }
  
  hatY = (diag(n) - rho*L)%*%Y
  MY = t(cbind(Uhat,Z))%*%hatY 
  betac  = solve(MZW - omega)%*%MY
   

   ## estimate beta and gamma

   beta = betac[1:d]
   gamma = betac[(d+1):(d+d2)]
   sigmasqhat = c((t(Y)%*%(M)%*%Y - 2*rho*t(Y)%*%(M)%*%YL + rho^2* t(YL)%*%(M)%*%YL)/n)

   G= L%*%solve(diag(n)-rho*L)
   H= G%*%((Uhat%*%beta) + (Z%*%gamma))

## Estimate the corrected Fisher information matrix

  Infmatrix = matrix(NA,d+d2+2,d+d2+2)
  Infmatrix[1:d,] = cbind(MW,t(Uhat)%*%Z,t(Uhat)%*%H,t(deltahat)%*%beta/sigmasqhat)
  Infmatrix[(d+1):(d+d2),] = cbind(t(Z)%*%Uhat,t(Z)%*%Z,t(Z)%*%H,rep(0,d2))
  Infmatrix[d+d2+1,] = cbind(t(H)%*%Uhat, t(H)%*%Z, t(H)%*%H + sigmasqhat*tr((G+t(G))%*%G),tr(G))
  Infmatrix[d+d2+2,] = cbind(t(beta)%*%deltahat,t(rep(0,d2)), tr(G), n/(2*sigmasqhat)-t(beta)%*%deltahat%*%beta/(sigmasqhat)^2)

  Infmatrix=Infmatrix/(sigmasqhat) 


### Estimate S matrix from score function for one observation

V=hatY-Z%*%gamma-Uhat%*%beta
VZ = hatY-Z%*%gamma
VU = hatY-Uhat%*%beta
scorevec = matrix(NA,d+d2+2,n)
for (k in 1:n){
scorevec[1:d,k] = -(-Uhat[k,]%*%t(VZ[k,]) + Uhat[k,]%*%t(Uhat[k,])%*%beta -deltahat%*%beta/n) /sigmasqhat
scorevec[(d+1):(d+d2),k] = -(-Z[k,]%*%t(VU[k,]) + Z[k,]%*%t(Z[k,])%*%gamma) /sigmasqhat
scorevec[d+d2+1,k] = (YL[k,]%*%t(V[k,]))/sigmasqhat-tr(G)/n
scorevec[d+d2+2,k] = -1/(2*sigmasqhat) + (V[k,]%*%t(V[k,])-t(beta)%*%deltahat%*%beta/n)/(2*sigmasqhat^2)
 }

sigmamatrix = scorevec%*%t(scorevec)


### Compute the asymptotic covariance matrix as I^{-1} S I^{-1}

### Both infmatrix and sigmamatrix should be divided by n. Therefore cov is multiplied by n
### hence no more division by n for standard errors.

  cov = solve(Infmatrix)%*%sigmamatrix%*%solve(Infmatrix)
  serho = sqrt(diag(cov)[d+d2+1])
  sebeta = sqrt(diag(cov)[1:d])
  segamma = sqrt(diag(cov)[(d+1):(d+d2)])


prho = 2*(1-pnorm(abs(rho)/serho))
pbeta = 2*(1-pnorm(abs(beta)/sebeta))
pgamma = 2*(1-pnorm(abs(gamma)/segamma))

### Return as list the estimates, SEs, and p-values of influence (rho), coefficeints of Z (betacovariate), and coefficient of 
### latent variables (betalatent, not very meaningful) along with estimate of variance of error term


return(list(influence = rho, SEinfluence =serho, pvaluerho=prho, betalatent =beta, SEbetalatent = sebeta, 
  pvaluelatent = pbeta, betacovariate=gamma,   SEbetacovariate=segamma, pvaluecovariate = pgamma, sigmasqhat=sigmasqhat,d=d,k=comms))
}












################ Leave one out predictions    ###############


### The following function computes leave-one-out prediction for response - returns the predicted values and also computes 
### predicted R-squared (meaningful for continuous outcomes)

predRsq <- function(Y,L,X,Z,d)
{
  

### Pre compute Uhat for the whole network
n = length(Y)
  spectra<-svd(X)
  Uhatmain<-spectra$u[,1:d]%*%diag(sqrt(spectra$d[1:d]))


datapt=vector()
pred =vector()


### start leave one out loop
  for ( i in 1:n){

    ## every index becomes test index once
  testind = i

  ## create the test data: one element for Y, Z, Uhatmain
  Ytest = Y[testind]
  Ztest = Z[testind,]
  Uhattest=Uhatmain[testind,]

  ## create the vector of network links between test index and rest of data
  XoS =  X[testind,-testind]
  #LoS= XoS/sum(XoS)
  LoS=L[testind,-testind]

  ## Create the training dataset with the remaining n-1 points
  Ytrain=Y[-testind]
  Xtrain =X[-testind,-testind]
  Ltrain = L[-testind,-testind]

  ### Alternative way to compute Ltrain

  #Dtrain = rowSums(Xtrain)
  #Dtrain[Dtrain<1]=1
  #Ltrain = Xtrain/Dtrain
  
  Ztrain = Z[-testind,]
  Uhattrain=Uhatmain[-testind,]

  #### Estimates with the training data including training portion of pre computed Uhat
  ests = sarMLZ(Ytrain,Ltrain,Xtrain,Ztrain,d=d,Uhat=Uhattrain)

  ### Prediction using the SAR model on the test data
  predY = sum(Ztest*ests$betacovariate) + sum(Uhattest*ests$betalatent) + ests$influence*sum(LoS*Ytrain)

  

  pred = c(pred,predY)
  datapt = c(datapt,Ytest)
}

sqerr= (datapt-pred)^2

sst = sum((datapt-mean(datapt))^2)

predr2 = 1-(sum(sqerr)/sst)
return(list(predRsq = predr2, pred=pred, testdata=datapt))

}








