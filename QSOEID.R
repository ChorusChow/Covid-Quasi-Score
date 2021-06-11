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

# Overall Parameters
## transmissibility rate of the epidemic (infectiousness profile)
Omega<-as.vector(c(data=NA,length = 25))
for (i in 1:25) {
  Omega[i]= (pgamma(i,shape = 2.5,scale = 3)-pgamma((i-1),shape = 2.5,scale = 3))/(pgamma(25,shape = 2.5,scale = 3))
}
## tau_0, pre-specified time point where before date tau_0, the MLE will be applied.
tau_0=5 #default
## bias_corr_const, the bias correction constant, default is bias_corr_const=1
bias_corr_const=exp(-0.001/2)

# Functions
## logistic function
logit<-function(x){
  result=log(x/(1-x))
  return(result)
}
##
greater.or.na <- function(myobj, myconst, tocat){
  if(is.na(myobj)) cat(myobj, tocat, '\n')
  return(is.na(myobj) | myobj > myconst)
}


############ Main Function ##########################

QSOEID<-function(Z,I){

      ##################################################

      ### Lambda[t]
      Lambda<-matrix(data = NA, nrow=1,ncol=T)
      Lambda[1]=I_0*Omega[1]
      for (t in 2:25) { Lambda[t]=(I_0*Omega[t])+ (I[1:(t-1)]%*%Omega[(t-1):1] ) }
      for (t in 26:T) { Lambda[t]= I[(t-1):(t-25)]%*%Omega[1:25]   }

      ##################################################

      # Estimating EstR[1:tau_0] based on MLE
      EstR<-matrix(data = NA, nrow=1,ncol=T)
      EstR[1]=max(1,(I[1]/(I_0*Omega[1]) ) )
      for (t in 2:tau_0) { EstR[t]=max(1,  (I[t]/( (I_0*Omega[t])+ (I[1:(t-1)]%*%Omega[(t-1):1] ) ) ) )  }

      #################################################
      
      ### barZ[t], arverage of Z[i,]'s from date 1 to date t.
      barZ<-Z
      for (t in 2:T) { barZ[t,]=colSums(Z[1:t,])/t }
      
      ### WTilde, weight  \Big( \sum_{i=1}^t (Z_i-\bar{Z})(Z_i-\bar{Z})^T \Big)^{-1}
      WTilde<-list()
      WTilde[[1]]<-matrix(data=0, nrow=NoCov, ncol=NoCov)
      for (t in 2:tau_0) {
        WTilde[[t]]<-matrix(data=0, nrow=NoCov, ncol=NoCov)
        for (i in 1:t) {
          WTilde[[t]]=WTilde[[t]]+(Z[i,]-barZ[t,])%*%t((Z[i,]-barZ[t,]))
        }
        WTilde[[t]]=solve( WTilde[[t]] +diag(c(rep(0.1,NoCov)),nrow=NoCov ))
      }
      for (t in 3:T) {
        WTilde[[t]]<-matrix(data=0, nrow=NoCov, ncol=NoCov)
        for (i in 1:t) {
          WTilde[[t]]=WTilde[[t]]+(Z[i,]-barZ[t,])%*%t((Z[i,]-barZ[t,]))
        }
        WTilde[[t]]=solve( WTilde[[t]] )
      }

      ###
      YTilde<-matrix(data = NA, nrow=1,ncol=T)
      for (t in 1:tau_0) { YTilde[t]=log(EstR[t])}
      BetaTilde<-matrix(data=NA, nrow=T, ncol=NoCov)
      ZYTilde<-matrix(data = NA, nrow=T,ncol=NoCov)
      ### EstPhi, EstBeta, for instore the estimated regression parameter
      EstPhi<-matrix(data = NA,nrow=T,ncol = 2)
      EstBeta<-matrix(data = NA, nrow=T,ncol = NoCov)
      ### Intermediate variables, ZYHat
      ZYHat=ZYTilde
      
      ### Define the profile likelihood
      ell<- function(phi,k){
        ZYTilde[(k-1),]=( log(EstR[,1]) -phi[2]*log(R_0) -phi[1]  )%*%(Z[1,]-barZ[(k-1),])
        EstR[,1]=exp(phi[1]+phi[2]*log(R_0))
        for (i in 2:(k-1)) {
          ZYTilde[(k-1),]=ZYTilde[(k-1),]+ (log(EstR[,i])-phi[2]*log(EstR[,(i-1)]) - phi[1]  ) %*%(Z[i,]-barZ[(k-1),])
        }
        BetaTilde[k,]=ZYTilde[(k-1),]%*%WTilde[[(k-1)]]
        for (i in (tau_0+1):k){
          YTilde[,i]=phi[1]+phi[2]*log(EstR[,(i-1)])+t(Z[i,])%*%BetaTilde[k,]
        }
        result=0
        for (j in (tau_0+1):k) {
          result=result+I[j]*YTilde[,j]-bias_corr_const*(exp(YTilde[,j]))*Lambda[j]
        }
        return(-result)
      }
 
      ### Minimize over the minus profile log-likelihood
      for (t in (tau_0+1):T) {
        EstPhi[t,]=nlminb(c(0.05,0.7),ell,k=t,lower = c(-5,0.3),upper = c(5,0.95))$par
        
        ZYHat[(t-1),]=( log(EstR[,1]) -EstPhi[t,2]*log(R_0) -EstPhi[t,1]  )%*%(Z[1,]-barZ[(t-1),])
        for (i in 2:(t-1)) {
          ZYHat[(t-1),]=ZYHat[(t-1),]+ (log(EstR[,i])-EstPhi[t,2]*log(EstR[,(i-1)]) - EstPhi[t,1]  ) %*%(Z[i,]-barZ[(t-1),])
        }
        EstBeta[t,]=ZYHat[(t-1),]%*%WTilde[[(t-1)]]
        
        EstR[,t]=exp( EstPhi[t,1]+EstPhi[t,2]*log(EstR[,(t-1)])+t(Z[t,])%*%EstBeta[t,] ) 
        
        if(greater.or.na(abs(EstR[,t]-EstR[,(t-1)]), 5, 'L182?')){
          EstPhi[t,]=EstPhi[(t-1),]
          EstBeta[t,]=EstBeta[(t-1),]
          EstR[,t]=EstR[,(t-1)]
        }
      }	
      
      #############################################
      
      #Bootstrap Estimator for Different iterations
      IndexBootstrap<-sample(1:(T+1-tunningl),size=rep, replace = TRUE)
      
      ## Define outcome variables EstRBoot, EstPhiBoot, EstBetaBoot
      EstRBoot<-matrix(data = NA, nrow=rep,ncol=tunningl)
      EstPhiBoot<-list()
      EstBetaBoot<-list()
      
      J<-matrix(data = NA, nrow=1,ncol=tunningl)  #intermediate variable, J is the bootstrap version data I.
      
      for (jBoot in 1:rep){
        if (IndexBootstrap[jBoot]==1){
          for (k in 1:tunningl) { EstRBoot[jBoot,k]=EstR[k] }
          EstPhiBoot[[jBoot]]=EstPhi[1:tunningl,]
          EstBetaBoot[[jBoot]]=EstBeta[1:tunningl,]
        }else{
          J_0<-I[1,(IndexBootstrap[jBoot]-1)]
          JR_0<-EstR[1,(IndexBootstrap[jBoot]-1)]
          for (i in 1:tunningl) {  J[1,i]=I[1,(IndexBootstrap[jBoot]+i-1)] }

          #     ### JLambda[t]
          JLambda<-matrix(data = NA, nrow=1,ncol=tunningl)
          for (t in 1:tunningl) { JLambda[t]=Lambda[(IndexBootstrap[jBoot]+t-1)] }
          
          EstRBoot[jBoot,1]=max(1,J[1])/(JLambda[1])
          if (EstRBoot[jBoot,1]-JR_0 >3 ){
            EstRBoot[jBoot,1]=JR_0+3
          }else if (EstRBoot[jBoot,1]-JR_0< (-3) ){
            EstRBoot[jBoot,1]=JR_0-3
          }
          for (t in 2:tau_0) {
            EstRBoot[jBoot,t]=max(1,J[t])/( JLambda[t] )
            if(greater.or.na(EstRBoot[jBoot,t]-EstRBoot[jBoot,(t-1)], 3, 'L234?')){
              EstRBoot[jBoot,t]=EstRBoot[jBoot,(t-1)]+3
            }else if(greater.or.na(-EstRBoot[jBoot,t]+EstRBoot[jBoot,(t-1)], 3, 'L237?')){
              EstRBoot[jBoot,t]=EstRBoot[jBoot,(t-1)]-3
            }
          }
          
          # Estimating EstRBoot[jBoot,t] based on profile likelihood
          ## Define intermediate variables, JYTilde, JBetaTilde, JWTilde, JZYTilde, JbarZ, JZ
          JZ<- Z[IndexBootstrap[jBoot]:(IndexBootstrap[jBoot]+tunningl-1),]
          JbarZ<-JZ
          for (t in 2:length(JZ[,1])) { JbarZ[t,]=colSums(JZ[1:t,])/t }
          
          ### JWTilde, weight  \Big( \sum_{i=1}^t (Z_i-\bar{Z})(Z_i-\bar{Z})^T \Big)^{-1} for data Z[IndexBootstrap[jBoot]: (...+tunningl-1),]
          JWTilde<-list()
          JWTilde[[1]]<-matrix(data=0, nrow=NoCov, ncol=NoCov)
          JWTilde[[2]]=(JZ[1,]-JbarZ[2,])%*%t((JZ[1,]-JbarZ[2,]))+(JZ[2,]-JbarZ[2,])%*%t((JZ[2,]-JbarZ[2,]))
          JWTilde[[2]]=solve( JWTilde[[2]]+diag(c(rep(0.1,(NoCov-1)),0.1),nrow=NoCov ) )
          for (t in 3:tau_0) {
            JWTilde[[t]]<-matrix(data=0, nrow=NoCov, ncol=NoCov)
            for (i in 1:t) {
              JWTilde[[t]]=JWTilde[[t]]+(JZ[i,]-JbarZ[t,])%*%t((JZ[i,]-JbarZ[t,]))
            }
            JWTilde[[t]]=solve( JWTilde[[t]] +diag(c(rep(0.01,NoCov)),nrow=NoCov ))
          }
          for (t in (tau_0+1):length(JZ[,1])) {
            JWTilde[[t]]<-matrix(data=0, nrow=NoCov, ncol=NoCov)
            for (i in 1:t) {
              JWTilde[[t]]=JWTilde[[t]]+(JZ[i,]-JbarZ[t,])%*%t((JZ[i,]-JbarZ[t,]))
            }
            JWTilde[[t]]=solve( JWTilde[[t]] )
          }
          
          ###
          JYTilde<-matrix(data = NA, nrow=1,ncol=tunningl)
          for (t in 1:tau_0) {  JYTilde[,t]=log(EstRBoot[jBoot,t])  }
          
          JBetaTilde<-matrix(data=NA, nrow=tunningl, ncol=NoCov)
          JZYTilde<-matrix(data = NA, nrow=tunningl,ncol=NoCov)
          
          ### EstPhiBoot, EstBetaBoot, for instore the estimated regression parameter
          EstPhiBoot[[jBoot]]<-matrix(data = NA,nrow=tunningl,ncol = 2)
          EstBetaBoot[[jBoot]]<-matrix(data = NA, nrow=tunningl,ncol = NoCov)
          
          ### Intermediate variables, JZYHat
          JZYHat=JZYTilde
          
          ### Define the profile likelihood
          ellBoot<- function(phi,k){
            JZYTilde[(k-1),]=( log(EstRBoot[jBoot,1]) -phi[2]*log(JR_0) -phi[1]  )%*%(JZ[1,]-JbarZ[(k-1),])
            EstRBoot[jBoot,1]=exp(phi[1]+phi[2]*log(JR_0))
            for (i in 2:(k-1)) {
              JZYTilde[(k-1),]=JZYTilde[(k-1),]+ (log(EstRBoot[jBoot,i])-phi[2]*log(EstRBoot[jBoot,(i-1)]) - phi[1]  ) %*%(JZ[i,]-JbarZ[(k-1),])
            }
            JBetaTilde[k,]=JZYTilde[(k-1),]%*%JWTilde[[(k-1)]]
            for (i in (tau_0+1):k){
              JYTilde[,i]=phi[1]+phi[2]*log(EstRBoot[jBoot,(i-1)])+t(JZ[i,])%*%JBetaTilde[k,]
            }
            result=0
            for (j in (tau_0+1):k) {
              result=result+J[j]*JYTilde[,j]-bias_corr_const*(exp(JYTilde[,j]))*JLambda[j]
            }
            return(-result)
          }
          
          #Bootstrap Estimation procedure
          for (t in (tau_0+1):tunningl) {
            EstPhiBoot[[jBoot]][t,]=nlminb(c(0.05,0.7),ellBoot,k=t,lower = c(-5,0.05),upper = c(5,0.95))$par
            
            JZYHat[(t-1),]=( log(EstRBoot[jBoot,1]) -EstPhiBoot[[jBoot]][t,2]*log(JR_0) -EstPhiBoot[[jBoot]][t,1]  )%*%(JZ[1,]-JbarZ[(t-1),])
            for (i in 2:(t-1)) {
              JZYHat[(t-1),]=JZYHat[(t-1),]+ (log(EstRBoot[jBoot,i])-EstPhiBoot[[jBoot]][t,2]*log(EstRBoot[jBoot,(i-1)]) - EstPhiBoot[[jBoot]][t,1]  ) %*%(JZ[i,]-JbarZ[(t-1),])
            }
            EstBetaBoot[[jBoot]][t,]=JZYHat[(t-1),]%*%JWTilde[[(t-1)]]
            
            EstRBoot[jBoot,t]=exp( EstPhiBoot[[jBoot]][t,1]+EstPhiBoot[[jBoot]][t,2]*log(EstRBoot[jBoot,(t-1)])+t(JZ[t,])%*%EstBetaBoot[[jBoot]][t,] )
            if(greater.or.na(abs(EstRBoot[jBoot,t]-EstRBoot[jBoot,(t-1)]), 5, 'L329?')){
              EstPhiBoot[[jBoot]][t,]=EstPhiBoot[[jBoot]][(t-1),]
              EstBetaBoot[[jBoot]][t,]=EstBetaBoot[[jBoot]][(t-1),]
              EstRBoot[jBoot,t]=EstRBoot[jBoot,(t-1)]
            }
          }
          
        } ### end of the else loop
      } ### end of the bootstrap loop
      
      restEstPhiBoot<-matrix(data = NA, nrow=2,ncol = rep)
      restEstBetaBoot<-matrix(data = NA, nrow = NoCov, ncol = rep)
      
      for (jBoot in 1:rep) {
        restEstPhiBoot[1,jBoot]=EstPhiBoot[[jBoot]][tunningl,1]
        restEstPhiBoot[2,jBoot]=EstPhiBoot[[jBoot]][tunningl,2]
        for (k in 1:NoCov){
          restEstBetaBoot[k,jBoot]=EstBetaBoot[[jBoot]][tunningl,k]
        }
      }
      
      
      BootstrapEstR<-matrix(data=NA,nrow=rep,ncol = T)
      for (jBootrep in 1:rep) {
        BootstrapEstR[jBootrep,1]=exp( restEstPhiBoot[1,jBootrep]+
                                         restEstPhiBoot[2,jBootrep]*log(R_0)+
                                         Z[1,]%*%restEstBetaBoot[,jBootrep] )
        for (Rt in 2:T) {
          BootstrapEstR[jBootrep,Rt]= exp( restEstPhiBoot[1,jBootrep]+
                                             restEstPhiBoot[2,jBootrep]*log( BootstrapEstR[jBootrep,(Rt-1)] )+
                                             Z[Rt,]%*%restEstBetaBoot[,jBootrep] )
        }
      }
      
      ### Calculate the bootstrap band
      UpperCIBand<-matrix(data=NA,nrow=1,ncol = T)
      LowerCIBand<-matrix(data=NA,nrow=1,ncol = T)
      upperindex=round(rep*0.95)
      lowerindex=max(1,rep*0.05)
      
      for (Rt in 1:T) {
        a=sort(BootstrapEstR[,Rt])
        UpperCIBand[1,Rt]=a[upperindex]
        LowerCIBand[1,Rt]=a[lowerindex]
      }
      
      ### Calculate the bootstrap confidence interval for parameters.
      
      restPara<-matrix(data = NA,nrow=3,ncol = 4)
      restPara[2,1]=EstPhi[T,1]
      restPara[2,2]=EstPhi[T,2]
      restPara[2,3]=EstBeta[T,1]
      restPara[2,4]=EstBeta[T,2]
      restPara[1,1]=sort(restEstPhiBoot[1,])[lowerindex]
      restPara[1,2]=sort(restEstPhiBoot[2,])[lowerindex]
      restPara[1,3]=sort(restEstBetaBoot[1,])[lowerindex]
      restPara[1,4]=sort(restEstBetaBoot[2,])[lowerindex]
      restPara[3,1]=sort(restEstPhiBoot[1,])[upperindex]
      restPara[3,2]=sort(restEstPhiBoot[2,])[upperindex]
      restPara[3,3]=sort(restEstBetaBoot[1,])[upperindex]
      restPara[3,4]=sort(restEstBetaBoot[2,])[upperindex]
      
      
      
      EstR=rbind(LowerCIBand,EstR,UpperCIBand)
      
      
      resQSOEID=list(restPara=restPara,EstR=EstR)
      return(resQSOEID)
}












