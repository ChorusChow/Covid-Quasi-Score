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



###############  generate misspecified data trial one #########################
Ztrial1<-matrix(data=NA, nrow=T, ncol=(NoCov+1) )
for (t in 1:T) {  Ztrial1[t,1]=10-(T/8)+(t/4)+rnorm(1,mean=0,sd=3) }
Ztrial1[,2]=logit( runif(T,min=0.01,max=0.99) )+2
for (t in 1:T) {  Ztrial1[t,3]=-(T/18)+(t/9)+rnorm(1,mean=0,sd=5) }
## Oracle beta
OracleBeta<-as.matrix(c(-0.02,-0.125,-0.03),nrow=3,ncol=1)
## R[t], the instantaneous reproduction number at time t.
Rtrial1<-matrix(data = NA, nrow=1,ncol=T)
epsilon<-rmvnorm(1,mean=rep(0,T),sigma = diag( (-2*log(bias_corr_const))  ,nrow = T,ncol = T))
Rtrial1[1]=exp( (OraclePhi[1]) + (OraclePhi[2]*log(R_0)) + (Ztrial1[1,]%*%OracleBeta) + (epsilon[1]) )
for (t in 2:T) { Rtrial1[t]=exp( (OraclePhi[1]) + (OraclePhi[2]*log(Rtrial1[(t-1)]))+ (Ztrial1[t,]%*%OracleBeta) + (epsilon[t]) ) }
## Generate I[t], number of incidence at time t.
Itrial1<-matrix(data = NA, nrow=1,ncol=T)
Itrial1[1]=rpois(1,lambda = (Rtrial1[1]*I_0*Omega[1]))
for (t in 2:25) { Itrial1[t]=rpois(1,lambda = ( Rtrial1[t]*((I_0*Omega[t])+ (Itrial1[1:(t-1)]%*%Omega[(t-1):1] )  )  ) ) }
for (t in 26:T) {
  if (( Itrial1[(t-1):(t-25)]%*%Omega[1:25]  )<=100000){
    Itrial1[t]=rpois(1,lambda = (Rtrial1[t]*( Itrial1[(t-1):(t-25)]%*%Omega[1:25]  ) ))
  }else {
    multipconstant= (Rtrial1[t]*( Itrial1[(t-1):(t-25)]%*%Omega[1:25]  ) )%/%100000  
    residueconstant= (Rtrial1[t]*( Itrial1[(t-1):(t-25)]%*%Omega[1:25]  ) )%%100000  
    Itrial1[t]=sum( rpois(multipconstant,lambda = 100000   ) ) + rpois(1,lambda = residueconstant  )
  }
}



save(Ztrial1,Itrial1,Rtrial1,file="misspecified_data_trial1.rda")



###############  generate misspecified data trial one #########################
Ztrial2<-matrix(data=NA, nrow=T, ncol=NoCov )
for (t in 1:T) {  Ztrial2[t,1]=7.5-(T/8)+(t/4)+rnorm(1,mean=0,sd=3) }
Ztrial2[,2]=logit( runif(T,min=0.01,max=0.99) )+2
## Oracle phi
OraclePhi<-as.matrix(c(0.5,0.5,0.3),nrow=3,ncol=1)
## Oracle beta
OracleBeta=as.matrix(c(-0.02,-0.125),nrow=NoCov,ncol=1)
## R[t], the instantaneous reproduction number at time t.
Rmin1=2
Rtrial2<-matrix(data = NA, nrow=1,ncol=T)
epsilon<-rmvnorm(1,mean=rep(0,T),sigma = diag( (-2*log(bias_corr_const))  ,nrow = T,ncol = T))
Rtrial2[1]=exp( (OraclePhi[1]) + (OraclePhi[2]*log(R_0)) +(OraclePhi[3]*log(Rmin1)) + (Ztrial2[1,]%*%OracleBeta) + (epsilon[1]) )
Rtrial2[2]=exp( (OraclePhi[1]) + (OraclePhi[2]*log(Rtrial2[1])) +(OraclePhi[3]*log(R_0)) + (Ztrial2[2,]%*%OracleBeta) + (epsilon[2]) )
for (t in 3:T) { 
  Rtrial2[t]=exp( (OraclePhi[1]) + (OraclePhi[2]*log(Rtrial2[(t-1)])) +(OraclePhi[3]*log(Rtrial2[(t-2)])) + (Ztrial2[t,]%*%OracleBeta) + (epsilon[t]) )
  }
## Generate I[t], number of incidence at time t.
Itrial2<-matrix(data = NA, nrow=1,ncol=T)
Itrial2[1]=rpois(1,lambda = (Rtrial2[1]*I_0*Omega[1]))
for (t in 2:25) { Itrial2[t]=rpois(1,lambda = ( Rtrial2[t]*((I_0*Omega[t])+ (Itrial2[1:(t-1)]%*%Omega[(t-1):1] )  )  ) ) }
for (t in 26:T) {
  if (( Itrial2[(t-1):(t-25)]%*%Omega[1:25]  )<=100000){
    Itrial2[t]=rpois(1,lambda = (Rtrial2[t]*( Itrial2[(t-1):(t-25)]%*%Omega[1:25]  ) ))
  }else {
    multipconstant= (Rtrial2[t]*( Itrial2[(t-1):(t-25)]%*%Omega[1:25]  ) )%/%100000  
    residueconstant= (Rtrial2[t]*( Itrial2[(t-1):(t-25)]%*%Omega[1:25]  ) )%%100000  
    Itrial2[t]=sum( rpois(multipconstant,lambda = 100000   ) ) + rpois(1,lambda = residueconstant  )
  }
}



save(Ztrial2,Itrial2,Rtrial2,file="misspecified_data_trial2.rda")











