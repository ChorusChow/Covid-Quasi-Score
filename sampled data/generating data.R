# load packages
library(pracma)
library(ggplot2)
library(EpiEstim)
library(nlme)
library(dlnm)
library(tsModel)
library(corrplot)
library(mvtnorm)
library(stats)
library(gridExtra)
library(ggridges)
library(lubridate)

# set working directory
setwd("~/Dropbox/working on folder/Covid/Code_for_model/Code_reorganize/sampled data/")
## transmissibility rate of the epidemic (infectiousness profile)
Omega<-as.vector(c(data=NA,length = 25))
for (i in 1:25) {
  Omega[i]= (pgamma(i,shape = 2.5,scale = 3)-pgamma((i-1),shape = 2.5,scale = 3))/(pgamma(25,shape = 2.5,scale = 3))
}
## logistic function
logit<-function(x){
  result=log(x/(1-x))
  return(result)
}

# Overall Parameters

#####Model 1, corrctly specified model##########################################################

## tau_0, pre-specified time point where before date tau_0, the MLE will be applied. Default tau_0=5.
tau_0=5
## NoCov, number of covariates we choose
NoCov=2
## T, number of observations/days
T=120
### R[0], the instantaneous reproduction number at time 0.
R_0=3
## Parameter of interest
OracleBeta=as.matrix(c(-0.02,-0.125),nrow=NoCov,ncol=1)
OraclePhi=c(0.5,0.7)
## I_0, start cases
I_0=500
## rep, number of replications for bootstrap
rep=200
## tunningl, the block length of a fragment of the time series (daily incident cases) used in bootstrap
tunningl=45
## bias_corr_const, the bias correction constant, default is bias_corr_const=1
bias_corr_const=exp(-0.001/2)

### generate daily incident cases and covariates data
Z<-matrix(data=NA, nrow=T, ncol=NoCov)
for (t in 1:T) {Z[t,1]=5-(T/8)+((2*t)/8)+rnorm(1,mean=0,sd=3)}
Z[,2]=logit( runif(T,min=0.01,max=0.99) )+2
## R[t], the instantaneous reproduction number at time t.
R<-matrix(data = NA, nrow=1,ncol=T)
epsilon<-rmvnorm(1,mean=rep(0,T),sigma = diag( (-2*log(bias_corr_const))  ,nrow = T,ncol = T))
R[1]=exp( (OraclePhi[1]) + (OraclePhi[2]*log(R_0)) + (Z[1,]%*%OracleBeta) + (epsilon[1]) )
for (t in 2:T) { R[t]=exp( (OraclePhi[1]) + (OraclePhi[2]*log(R[(t-1)]))+ (Z[t,]%*%OracleBeta) + (epsilon[t]) ) }
## Generate I[t], number of incidence at time t.
I<-matrix(data = NA, nrow=1,ncol=T)
I[1]=rpois(1,lambda = (R[1]*I_0*Omega[1]))
for (t in 2:25) { I[t]=rpois(1,lambda = ( R[t]*((I_0*Omega[t])+ (I[1:(t-1)]%*%Omega[(t-1):1] )  )  ) ) }
for (t in 26:T) {
  if (( I[(t-1):(t-25)]%*%Omega[1:25]  )<=100000){
    I[t]=rpois(1,lambda = (R[t]*( I[(t-1):(t-25)]%*%Omega[1:25]  ) ))
  }else {
    multipconstant= (R[t]*( I[(t-1):(t-25)]%*%Omega[1:25]  ) )%/%100000  
    residueconstant= (R[t]*( I[(t-1):(t-25)]%*%Omega[1:25]  ) )%%100000  
    I[t]=sum( rpois(multipconstant,lambda = 100000   ) ) + rpois(1,lambda = residueconstant  )
  }
}



save(Z,I,R,file="sampled_data.rda")
















