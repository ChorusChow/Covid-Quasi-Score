library(ggplot2)
library(EpiEstim)
library(nlme)
library(dlnm)
library(dplyr)
library(tsModel)
library(sjstats)
library(corrplot)
library(pracma)
library(usmap)

# set working directory
setwd('/Users/shij/Desktop/Covid_code_real_data_organize/')
# source R function
source('./functions/covid-function_reorganized_0706.R')

load('covid_real_res_20200706_v2vari5.RData')

######################## Calculate EstR based on estimation #####################################

{
  ### Calculate EstR based on estimation
  for (i in 2:nrow(reg.dat.est2v5) ) {
    
    if ( is.na(reg.dat.est2v5$EstR[i]) & reg.dat.est2v5$FIPS[i]==
         reg.dat.est2v5$FIPS[(i-1)] ){
      phi0inter=EstPara[(as.numeric(reg.dat.est2v5$date[i]-pre.end.date)),1]
      phi1inter=EstPara[(as.numeric(reg.dat.est2v5$date[i]-pre.end.date)),2]
      betainter=EstPara[(as.numeric(reg.dat.est2v5$date[i]-pre.end.date)),3:7]
      
      reg.dat.est2v5$EstR[i] = exp( phi0inter + 
                                      phi1inter* log(reg.dat.est2v5$EstR[(i-1)])+
                                      betainter[1] * reg.dat.est2v5$PopO64[i]+
                                      betainter[2] * reg.dat.est2v5$popdens_s[i]+
                                      betainter[3] * reg.dat.est2v5$diabetes_s[i]+
                                      betainter[4] * reg.dat.est2v5$visit_mean[i]+
                                      betainter[5] * reg.dat.est2v5$temp_mean[i]  ) 
      
      if ( ( reg.dat.est2v5$EstR[i]-reg.dat.est2v5$EstR[(i-1)] )>5  ){
        reg.dat.est2v5$EstR[i]= reg.dat.est2v5$EstR[(i-1)] + 5
      }else if (( reg.dat.est2v5$EstR[i]-reg.dat.est2v5$EstR[(i-1)] ) < -5 ){
        reg.dat.est2v5$EstR[i]=reg.dat.est2v5$EstR[(i-1)] - 5
      }
      
    }
    
  }
  
}

######################## Bootstrap with block length ellboot=45 ##################################

{
  
  ellboot=45
  NumRep=Odays-ellboot
  EstParaBoot45=list()
  
  for (i in 1:NumRep) {
    Start.date.boot=Start.date.est+i
    End.date.boot=Start.date.boot+ellboot-1
    pre.end.boot=Start.date.boot+5
    
    reg.dat.boot2v5=reg.dat.est2v5[reg.dat.est2v5$date>=Start.date.boot & 
                                     reg.dat.est2v5$date<=End.date.boot, ]
    
    reg.dat.summarize.boot <- reg.dat.boot2v5%>% group_by(FIPS)%>%
      summarize(county = unique(county), state= unique(STATE_NAME), 
                max_cases = max(cases,na.rm = TRUE),ob_days= sum(cases>0),
                min_ob_day=min(date),max_ob_day=max(date))
    reg.dat.summarize.boot=reg.dat.summarize.boot[reg.dat.summarize.boot$ob_days>6,]
    
    reg.dat.boot2v5=reg.dat.boot2v5[reg.dat.boot2v5$FIPS %in% reg.dat.summarize.boot$FIPS,]
    
    reg.dat.boot2v5$EstRboot=NA
    
    NumCounty=nrow(reg.dat.summarize.boot)
    for (j in 1:NumCounty) {
      reg.dat.boot2v5[reg.dat.boot2v5$FIPS==reg.dat.summarize.boot$FIPS[j] & 
                        reg.dat.boot2v5$date==reg.dat.summarize.boot$min_ob_day[j],]$EstRboot =
        reg.dat.boot2v5[reg.dat.boot2v5$FIPS==reg.dat.summarize.boot$FIPS[j] & 
                          reg.dat.boot2v5$date==reg.dat.summarize.boot$min_ob_day[j],]$EstR
    }
    
    reg.dat.boot2v5[reg.dat.boot2v5$date<=pre.end.boot,]$EstRboot=
      reg.dat.boot2v5[reg.dat.boot2v5$date<=pre.end.boot,]$EstR
    
    EstParaBoot45[[i]]<-matrix(data = NA,nrow=as.numeric(End.date.boot-pre.end.boot),ncol = (2+NumCov))
    
    time0=Sys.time()
    for (j in 1:as.numeric(End.date.boot-pre.end.boot)) {
      #for (j in 1:3) {
      EstParaBoot45[[i]][j,]=estR.iterative.Boot(reg.dat.boot2v5,dayt=j)$para
      reg.dat.boot2v5[reg.dat.boot2v5$date <= (pre.end.boot+j),]$EstRboot=
        (estR.iterative.Boot(reg.dat.boot2v5,dayt=j)$my_dat_est)$EstRboot
      cat(j)
    }
    time1=Sys.time()-time0
    time1
    
    saveRDS(EstParaBoot[[i]], file=paste0('derived/result_ellboot45/covid_real_res_bs_20200721_v2vari5_ellboot', ellboot, '_i', i, '.RDS'))
    
  }

}


######################## Bootstrap with block length ellboot=60 ##################################

{
  
  ellbootv2=60
  NumRepv2=Odays-ellbootv2
  EstParaBoot60=list()
  
  for (i in 1:NumRepv2) {
    Start.date.boot=Start.date.est+i
    End.date.boot=Start.date.boot+ellbootv2-1
    pre.end.boot=Start.date.boot+5
    
    reg.dat.boot2v5=reg.dat.est2v5[reg.dat.est2v5$date>=Start.date.boot & 
                                     reg.dat.est2v5$date<=End.date.boot, ]
    
    reg.dat.summarize.boot <- reg.dat.boot2v5%>% group_by(FIPS)%>%
      summarize(county = unique(county), state= unique(STATE_NAME), 
                max_cases = max(cases,na.rm = TRUE),ob_days= sum(cases>0),
                min_ob_day=min(date),max_ob_day=max(date))
    reg.dat.summarize.boot=reg.dat.summarize.boot[reg.dat.summarize.boot$ob_days>6,]
    
    reg.dat.boot2v5=reg.dat.boot2v5[reg.dat.boot2v5$FIPS %in% reg.dat.summarize.boot$FIPS,]
    
    reg.dat.boot2v5$EstRboot=NA
    
    NumCounty=nrow(reg.dat.summarize.boot)
    for (j in 1:NumCounty) {
      reg.dat.boot2v5[reg.dat.boot2v5$FIPS==reg.dat.summarize.boot$FIPS[j] & 
                        reg.dat.boot2v5$date==reg.dat.summarize.boot$min_ob_day[j],]$EstRboot =
        reg.dat.boot2v5[reg.dat.boot2v5$FIPS==reg.dat.summarize.boot$FIPS[j] & 
                          reg.dat.boot2v5$date==reg.dat.summarize.boot$min_ob_day[j],]$EstR
    }
    
    reg.dat.boot2v5[reg.dat.boot2v5$date<=pre.end.boot,]$EstRboot=
      reg.dat.boot2v5[reg.dat.boot2v5$date<=pre.end.boot,]$EstR
    
    EstParaBoot60[[i]]<-matrix(data = NA,nrow=as.numeric(End.date.boot-pre.end.boot),ncol = (2+NumCov))
    
    time0=Sys.time()
    for (j in 1:as.numeric(End.date.boot-pre.end.boot)) {
      #for (j in 1:3) {
      EstParaBoot60[[i]][j,]=estR.iterative.Boot(reg.dat.boot2v5,dayt=j)$para
      reg.dat.boot2v5[reg.dat.boot2v5$date <= (pre.end.boot+j),]$EstRboot=
        (estR.iterative.Boot(reg.dat.boot2v5,dayt=j)$my_dat_est)$EstRboot
      cat(j)
    }
    time1=Sys.time()-time0
    time1
    
    saveRDS(EstParaBoot60[[i]], file=paste0('derived/result_ellboot60/covid_real_res_bs_20200706_v2vari5_ellboot', ellboot, '_i', i, '.RDS'))
  
  }

}

######################## Calculate Bootstrap Confidence Interval for each Parameter ##################################

{
  
  ### Take the case ellboot=45 as an example
  
  BootSampleSize=50
  EstBootInter=matrix(data = NA, nrow=BootSampleSize, ncol = (2+NumCov) )
  
  for (i in 1:BootSampleSize) {
    interk=sample(c(1:NumRep),1)
    interPara=readRDS(file = paste0('derived/result_ellboot45/covid_real_res_bs_20200706_v2vari5_ellboot45_i', interk, '.RDS'))
    EstBootInter[i,]=interPara[nrow(interPara),]
  }
  
  ### Each column of EstBCI is the 90% Bootstrap confidence interval (Lower bound, original estimator, upbound) for each parameter
  
  EstBCI=matrix(data = NA, nrow=3, ncol = (2+NumCov) )
  EstBCI[2,]=EstPara[nrow(EstPara),]
  
  upindex=floor(BootSampleSize*0.95)
  lowerindex=max( floor(BootSampleSize*0.05), 1)+1
  
  for (i in 1:(2+NumCov)) {
    EstBCI[1,i]=sort( EstBootInter[,i] )[lowerindex]
    EstBCI[3,i]=sort( EstBootInter[,i] )[upindex]
  }
  
}

######################## Calculate Bootstrap Confidence Interval for each county's instantaneous reproduction number #####################################

{
  
  reg.BCI=matrix(data = NA, nrow = nrow(reg.dat.est2v5), ncol = BootSampleSize)
  ########################################## Upbound for Bootstrap R0 #############
  
  for (Bootindex in 1:BootSampleSize) {
    
    reg.dat.test=reg.dat.est2v5
    
    for (i in 2:nrow(reg.dat.test) ) {
      
      if ( is.na(reg.dat.test$BCIupper[i]) & reg.dat.test$FIPS[i]==
           reg.dat.test$FIPS[(i-1)] ){
        phi0inter=EstBootInter[Bootindex,1]
        phi1inter=EstBootInter[Bootindex,2]
        betainter=EstBootInter[Bootindex,3:7]
        
        reg.dat.test$BCIupper[i] = exp( phi0inter + 
                                              phi1inter* log(reg.dat.test$BCIupper[(i-1)])+
                                              betainter[1] * reg.dat.test$PopO64[i]+
                                              betainter[2] * reg.dat.test$popdens_s[i]+
                                              betainter[3] * reg.dat.test$diabetes_s[i]+
                                              betainter[4] * reg.dat.test$visit_mean[i]+
                                              betainter[5] * reg.dat.test$temp_mean[i]  ) 
        
        if ( ( reg.dat.test$BCIupper[i]-reg.dat.test$BCIupper[(i-1)] )>5  ){
          reg.dat.test$BCIupper[i]= reg.dat.test$BCIupper[(i-1)] + 5
        }else if (( reg.dat.test$BCIupper[i]-reg.dat.test$BCIupper[(i-1)] ) < -5 ){
          reg.dat.test$BCIupper[i]=reg.dat.test$BCIupper[(i-1)] - 5
        }
        
      }
      
    }
    
    reg.BCI[,Bootindex]=reg.dat.test$BCIupper 
  
  }
  
  for (i in 1:nrow(reg.dat.est2v5)  ) {
    reg.dat.est2v5$BCIlower[i]=sort( reg.BCI[i,] )[lowerindex]
    reg.dat.est2v5$BCIupper[i]=sort( reg.BCI[i,] )[upindex]
  }
  
}

######################################

saveRDS(list(EstBCI,reg.dat.est2v5), file=paste0('covid_real_res_bs_20200706_v2vari5.RDS'))





