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
setwd('~/Desktop/Covid_code_real_data_organize')
# source R function
source('./functions/covid-function_reorganized_0706.R')


#################  Data Preparation  ####################
{
### specify assumptions
# Assumption 1: generation time
GT.ncov  <- discr_si(seq(1, 30), 7.5, 3.4)
GT.ncov=c(GT.ncov,(1-sum(GT.ncov))/2,(1-sum(GT.ncov))/2,rep(0,100))
# Assumptions 2: sliding window, 3-5 days
nd <- 3
# Assumptions 3: prior mean of R = 2.5
meanprior <- 2.5

### load the data
reg.dat<- read.csv('./rawdata/temp_casetest_soc_7_7.csv')
reg.dat$date <- as.Date(reg.dat$date,"%Y-%m-%d")
reg.dat$date <- as.Date(reg.dat$date,format = '%Y-%m-%d')

#drop Puerto Rico
#reg.dat <- reg.dat[reg.dat$FIPS!=72127,]
reg.dat <- reg.dat[reg.dat$STATE_NAME!="Rhode Island",]

names(reg.dat)

# standardize county level covariates
tt.dat <- distinct(reg.dat[,c('fips','DIABETES','OBESE','PopO64','Uninsured',
                              'PovRatioLT2','SMOKE',"Crowded",
                              'PopDens')])
tt.dat$popdens_s <- scale(log(tt.dat$PopDens))
tt.dat$popdens_s2 <- log(tt.dat$PopDens)


reg.dat <- left_join(reg.dat,tt.dat)

## scale population density after log
reg.dat$diabetes_s <- reg.dat$DIABETES/100
reg.dat$obese_s <- reg.dat$OBESE/100
reg.dat$smoke_s <- reg.dat$SMOKE/100

## generate lag socialing distancing measure (7 days ago)
reg.dat$visit_mean <- runMean(reg.dat$daily_visitation_diff,lags=4:14,group=reg.dat$fips)

## center wet_bulb_temp_C at 12, taking lag operator and stored in temp_mean
mean(na.omit(reg.dat$wet_bulb_temp_C))
reg.dat$temp <- reg.dat$wet_bulb_temp_C-12
reg.dat$temp_mean <- runMean(reg.dat$temp,lags=4:14,group=reg.dat$FIPS)

###Delete rows that are NA's
reg.dat <- reg.dat[rowSums(is.na(reg.dat)) != ncol(reg.dat), ]

### Check systemic error
reg.dat.est<- reg.dat[reg.dat$cases>=0 & !is.na(reg.dat$cases),]

}


#################  Missing Data Complement  ####################

{

### reg.dat.est with cases >5
reg.dat.est<- reg.dat[reg.dat$cases>5 & !is.na(reg.dat$cases),]
na_flag <- apply(is.na(reg.dat.est), 2, sum)
na_flag

### Keep covariate temp_mean and remove wet_bulb_temp_F, wet_bulb_temp_C, dry_bulb_temp
### rel_humid, abs_humid_est, abs_humid_R, pressure_in, temp_imp, temp
reg.dat.est = reg.dat.est[,!names(reg.dat.est) %in% 
                            c("wet_bulb_temp_F", "wet_bulb_temp_C","temp","dry_bulb_temp",
                              "rel_humid","abs_humid_est","abs_humid_R","pressure_in","temp_imp",
                              "DIABETES","SMOKE","OBESE","PopDens","daily_visitation_diff",
                              "daily_distance_diff","weekday","Positive","Negative",
                              "deaths","No_Masking_Requirement")]


{
  ######################################################################################
  #standardize county level covariates
  relative.dat.pic <- distinct(reg.dat[,c('DIABETES','OBESE','PopO64','Uninsured',
                                'PovRatioLT2','SMOKE',"Crowded","TotPop",
                                'PopDens',"temp_mean","visit_mean")])
  relative.dat.pic=relative.dat.pic[!is.na(relative.dat.pic$temp_mean)
                                      & !is.na(relative.dat.pic$visit_mean),]
  col3<-colorRampPalette(c("blue","gray90","violetred") )

  corrplot.mixed(cor(relative.dat.pic[,!names(relative.dat.pic) %in% c("FIPS")]),
                 is.corr=FALSE, lower = "color", upper = "number",order="FPC",
                tl.col="black",lower.col=col3(20), cl.lim=c(-0.5,1), upper.col = col3(20),
                tl.pos = "lt",diag = "l",tl.cex=0.9)
  ######################################################################################
}


### Based on the covariate heat pic, select visit_mean,popdens_s,diabetes_s,PopO64,PovRatioLT2
### temp_mean only and discard other covariates

reg.dat.est = reg.dat.est[,!names(reg.dat.est) %in% 
                            c( "Uninsured","Crowded","popdens_s2","obese_s","smoke_s","covid","newdeath",
                               "county_vmt","jan_avg_vmt","Localities_Permitted_to_Require_Masking","COUNTY_NAME",
                               "Statewide_Masking_Requirement","count","State","TotPop","PovRatioLT2","X")]

### Summarize data set
reg.dat.summarize.est <- reg.dat.est%>% group_by(fips)%>%
  summarize(county = unique(county), state= unique(STATE_NAME), 
            max_cases = max(cases,na.rm = TRUE),ob_days= sum(cases>0),
            min_ob_day=min(date),max_ob_day=max(date))
View(reg.dat.summarize.est)
na_flag <- apply(is.na(reg.dat.est), 2, sum)
#na_flag

### Dealing with visit_mean missing
reg.dat.est.missing.visit=reg.dat.est[is.na(reg.dat.est$visit_mean),]
reg.dat.summarize.est.missing.visit <- reg.dat.est.missing.visit%>% group_by(fips)%>%
  summarize(county = unique(county), state= unique(STATE_NAME), 
            max_cases = max(cases,na.rm = TRUE)-min(cases,na.rm = TRUE),
            miss_days= sum(is.na(visit_mean))  )
#View(reg.dat.summarize.est.missing.visit)

### indic.M being the true/false matrix indicting whether the visit_mean missing happens at
### the begining of that county's data.
indic.M<-matrix(data=NA,nrow=nrow(reg.dat.summarize.est.missing.visit),ncol=1)
  for (i in 1:nrow(reg.dat.summarize.est.missing.visit) ) {
    indic.M[i,1]=( reg.dat.est[reg.dat.est$fips == reg.dat.summarize.est.missing.visit$fips[i] &
                                 is.na(reg.dat.est$visit_mean), ]$date[1] ==
                     reg.dat.est[reg.dat.est$fips == reg.dat.summarize.est.missing.visit$fips[i],]$date[1] )
  }
indic.M
  for (i in 1:nrow(reg.dat.summarize.est.missing.visit) ) {
       indic.M[i,1]=( sum(  is.na(  reg.dat.est[reg.dat.est$fips == 
                                                reg.dat.summarize.est.missing.visit$fips[i],]$
                                     visit_mean[1:reg.dat.summarize.est.missing.visit$miss_days[i] ] ) ) ==
                        reg.dat.summarize.est.missing.visit$miss_days[i]   )
  }
indic.M

reg.dat.est=reg.dat.est[!is.na(reg.dat.est$visit_mean),]

na_flag <- apply(is.na(reg.dat.est), 2, sum)
na_flag



reg.dat.summarize.est <- reg.dat.est%>% group_by(fips)%>%
  summarize(county = unique(county), state= unique(STATE_NAME), 
            max_cases = max(cases,na.rm = TRUE),ob_days= sum(cases>0),
            min_ob_day=min(date),max_ob_day=max(date))

}



#################  Data Screening  ####################

{
  ### Check whether the data is consecutive
  
  reg.dat.est$date <- as.Date(reg.dat.est$date,"%Y-%m-%d")
  reg.dat.est$date <- as.Date(reg.dat.est$date,format = '%Y-%m-%d')  
  
  indic.M=matrix(data=NA,ncol=1,nrow = nrow(reg.dat.est))
  indic.M[1,1]=0
  for (i in 2:nrow(reg.dat.est)) {
    if ( reg.dat.est$fips[i] != reg.dat.est$fips[(i-1)]  ){
      indic.M[i,1]=0
    }else if ( as.numeric( reg.dat.est$date[i]) - as.numeric(reg.dat.est$date[(i-1)]) ==1   ){
      indic.M[i,1]=0
    }else if ( as.numeric( reg.dat.est$date[i]) - as.numeric(reg.dat.est$date[(i-1)]) >1  ){
      indic.M[i,1]=1
    }
  }

  which(indic.M==1) ### Line 42812, fips==42007, date= 2020-07-04
  a1=reg.dat.est[(which(indic.M==1)-1),]
  a1$date=a1$date+1
  a2=reg.dat.est[which(indic.M==1),]
  a1$cases=a2$cases-a2$newcase
  a1$newcase=a1$cases-reg.dat.est[(which(indic.M==1)-1),]$cases
  a1$visit_mean=(a1$visit_mean+a2$visit_mean)/2
  a1$temp_mean=(a1$temp_mean+a2$temp_mean)/2
  a3=reg.dat.est[1:(which(indic.M==1)-1),]
  a4=reg.dat.est[which(indic.M==1):nrow(reg.dat.est),]
  a5=rbind(a3,a1)
  a6=rbind(a5,a4)
  reg.dat.est=a6
  ### Special treatment for County 6025,8001,42003,51041
  #reg.dat.est=reg.dat.est[! ( reg.dat.est$FIPS==6025 & reg.dat.est$date<= "2020-03-12"),]
  #reg.dat.est=reg.dat.est[! ( reg.dat.est$FIPS==8001 & reg.dat.est$date<= "2020-03-14"),]
  #reg.dat.est=reg.dat.est[! ( reg.dat.est$FIPS==42003 & reg.dat.est$date<= "2020-03-15"),]
  #reg.dat.est=reg.dat.est[! ( reg.dat.est$FIPS==51041 & reg.dat.est$date<= "2020-03-19"),]
  
  reg.dat.est[reg.dat.est$newcase<0,]$newcase=0
  
  {
    ## Calculate incident cases / moving average of 4 days' newvases
    reg.dat.est$incident_cases=runMean(reg.dat.est$newcase,lags=0:3,group=reg.dat.est$fips)
    reg.dat.summarize.est <- reg.dat.est%>% group_by(fips)%>%
      summarize(county = unique(county), state= unique(STATE_NAME), 
                max_cases = max(cases,na.rm = TRUE),ob_days= sum(cases>0),
                min_ob_day=min(date),max_ob_day=max(date))
    
    for (i in 1:nrow(reg.dat.summarize.est)) {
      reg.dat.est[reg.dat.est$fips==reg.dat.summarize.est$fips[i],]$incident_cases[1]=
        reg.dat.est[reg.dat.est$fips==reg.dat.summarize.est$fips[i],]$cases[1]
      reg.dat.est[reg.dat.est$fips==reg.dat.summarize.est$fips[i],]$incident_cases[2]=
        reg.dat.est[reg.dat.est$fips==reg.dat.summarize.est$fips[i],]$incident_cases[4]
      reg.dat.est[reg.dat.est$fips==reg.dat.summarize.est$fips[i],]$incident_cases[3]=
        reg.dat.est[reg.dat.est$fips==reg.dat.summarize.est$fips[i],]$incident_cases[4]
    }
  }

  ### Calculate Lambda for reg.dat.edt
  #reg.dat.est$Lambda=reg.dat.est$incident_cases
  
  names(reg.dat.est)[1]="FIPS"
  
  #reg.dat.est=Lambda.gen.est(my.dat.est = reg.dat.est, Omega.nCov = GT.ncov, meanprior)
  
  #write.csv(reg.dat.est,paste('./derived/reg_dat_est_',Sys.Date(),'.csv',sep=""),row.names=F)
  
  ### Define EstR, BCIupper, BCIlower and calculate the I_0
  reg.dat.est$EstR=NA
  reg.dat.est$BCIupper=NA
  reg.dat.est$BCIlower=NA
  
  reg.dat.summarize.est <- reg.dat.est%>% group_by(FIPS)%>%
    summarize(county = unique(county), state= unique(STATE_NAME), 
              max_cases = max(cases,na.rm = TRUE),ob_days= sum(cases>0),
              min_ob_day=min(date),max_ob_day=max(date))
  
  NumCounty=nrow(reg.dat.summarize.est)
  for (i in 1:NumCounty) {
    reg.dat.est[reg.dat.est$FIPS==reg.dat.summarize.est$FIPS[i],]$EstR[1]=meanprior
  }
  
  reg.dat.est[!is.na(reg.dat.est$EstR),]$BCIupper=reg.dat.est[!is.na(reg.dat.est$EstR),]$EstR
  reg.dat.est[!is.na(reg.dat.est$EstR),]$BCIlower=reg.dat.est[!is.na(reg.dat.est$EstR),]$EstR
  
  ######################## Second Data Screening 
  
  ### select rows that have incident cases >=10 at least twice for every 5 days in a row..
  int.med.reg.dat.est=reg.dat.est[reg.dat.est$incident_cases >=10,]
  int.med.reg.dat.est$index=0
  for (i in 1: (nrow(int.med.reg.dat.est)-1) ) {
    if ( int.med.reg.dat.est$FIPS[i] != int.med.reg.dat.est$FIPS[(i+1)]  ){
      int.med.reg.dat.est$index[i]=2
    }else if ( as.numeric( int.med.reg.dat.est$date[(i+1)]) - 
               as.numeric(int.med.reg.dat.est$date[i]) <=4   ){
      int.med.reg.dat.est$index[i]=0
    }else {
      int.med.reg.dat.est$index[i]=1
    }
  }

  ### Summarize 
  indix=int.med.reg.dat.est[int.med.reg.dat.est$index==0,]
  int.med.reg.dat.summarize.est <- indix%>% group_by(FIPS)%>%
    summarize(county = unique(county), state= unique(STATE_NAME), 
              max_cases = max(cases),ob_days= sum(cases>0),
              min_ob_day=min(date),max_ob_day=max(date)
    )
  indix=int.med.reg.dat.summarize.est[int.med.reg.dat.summarize.est$ob_days>=20,]
  indix2=int.med.reg.dat.est[int.med.reg.dat.est$index==0 & int.med.reg.dat.est$FIPS %in%
                               indix$FIPS,]
  int.med.reg.dat.summarize.est <- indix2%>% group_by(FIPS)%>%
    summarize(county = unique(county), state= unique(STATE_NAME), 
              max_cases = max(cases),ob_days= sum(cases>0),
              min_ob_day=min(date),max_ob_day=max(date)
    )

  ### Version 2 of reg.dat.est
  reg.dat.est2v5=reg.dat.est[reg.dat.est$FIPS %in% int.med.reg.dat.summarize.est$FIPS,]
  
  for (i in 1:nrow(int.med.reg.dat.summarize.est)) {
    reg.dat.est2v5=reg.dat.est2v5[ !(reg.dat.est2v5$FIPS==int.med.reg.dat.summarize.est$FIPS[i] &
                                   reg.dat.est2v5$date<(int.med.reg.dat.summarize.est$min_ob_day[i]-1) ),] 
  }
  
  indix3=int.med.reg.dat.est[int.med.reg.dat.est$index==2,]
  indix3=indix3[indix3$FIPS %in% int.med.reg.dat.summarize.est$FIPS,]

  for (i in 1:nrow(indix3)) {
    if (indix3$FIPS[i]!=int.med.reg.dat.summarize.est$FIPS[i]){
      cat(m)
    }else if (indix3$date[i]<=( int.med.reg.dat.summarize.est$max_ob_day[i]  +4)){
      int.med.reg.dat.summarize.est$max_ob_day[i]<- indix3$date[i]
    }
  }
  
  for (i in 1:nrow(int.med.reg.dat.summarize.est)) {
    reg.dat.est2v5=reg.dat.est2v5[ !(reg.dat.est2v5$FIPS==int.med.reg.dat.summarize.est$FIPS[i] &
                                   reg.dat.est2v5$date>int.med.reg.dat.summarize.est$max_ob_day[i] ),] 
  }


  #########3### Check whether the data is consecutive again for reg.dat.est2v5
  indic.M=matrix(data=NA,ncol=1,nrow = nrow(reg.dat.est2v5))
  indic.M[1,1]=0
  for (i in 2:nrow(reg.dat.est2v5)) {
    if ( reg.dat.est2v5$FIPS[i] != reg.dat.est2v5$FIPS[(i-1)]  ){
      indic.M[i,1]=0
    }else if ( as.numeric( reg.dat.est2v5$date[i]) - as.numeric(reg.dat.est2v5$date[(i-1)]) ==1   ){
      indic.M[i,1]=0
    }else if ( as.numeric( reg.dat.est2v5$date[i]) - as.numeric(reg.dat.est2v5$date[(i-1)]) >1  ){
      indic.M[i,1]=1
    }
  }
  
  which(indic.M==1) 
  
}


#################  Estimation Setup  ####################

{
  ### Calculate incident cases for reg.dat.est2v5
  reg.dat.est2v5$incident_cases[1]<- reg.dat.est2v5$cases[1]
  for (i in 2:nrow(reg.dat.est2v5)) {
    if ( reg.dat.est2v5$FIPS[i] !=  reg.dat.est2v5$FIPS[(i-1)] ){
      reg.dat.est2v5$incident_cases[i]=reg.dat.est2v5$cases[i]
    }
  }
  
  reg.dat.est2v5[reg.dat.est2v5$incident_cases<0,]
  
  ### Calculate Lambda for reg.dat.edt1
  reg.dat.est2v5$Lambda=reg.dat.est2v5$incident_cases
  
  reg.dat.est2v5=Lambda.gen.est(my.dat.est = reg.dat.est2v5, Omega.nCov = GT.ncov, meanprior)
  
  write.csv(reg.dat.est2v5,paste('./derived/reg_dat_est_',Sys.Date(),'.csv',sep=""),row.names=F)
  
  ### Calculate the I_0
  
  reg.dat.summarize.est <- reg.dat.est2v5%>% group_by(FIPS)%>%
    summarize(county = unique(county), state= unique(STATE_NAME), 
              max_cases = max(cases,na.rm = TRUE),ob_days= sum(cases>0),
              min_ob_day=min(date),max_ob_day=max(date))
  
  NumCounty=nrow(reg.dat.summarize.est)
  for (i in 1:NumCounty) {
    reg.dat.est2v5[reg.dat.est2v5$FIPS==reg.dat.summarize.est$FIPS[i],]$EstR[1]=meanprior
  }
  
  reg.dat.est2v5[!is.na(reg.dat.est2v5$EstR),]$BCIupper=reg.dat.est2v5[!is.na(reg.dat.est2v5$EstR),]$EstR
  reg.dat.est2v5[!is.na(reg.dat.est2v5$EstR),]$BCIlower=reg.dat.est2v5[!is.na(reg.dat.est2v5$EstR),]$EstR

  ### Start date 
  Start.date.est = min(reg.dat.est2v5$date)
  ### End date
  End.date.est = max(reg.dat.est2v5$date)
  ### Days of instantaneous reproduction numbers to be estimated
  Odays=as.numeric(max(reg.dat.est2v5$date)-min(reg.dat.est2v5$date))+1
  
  ### Summarize
  
  County.Number<-reg.dat.est2v5[reg.dat.est2v5$FIPS==reg.dat.est2v5[reg.dat.est2v5$date==Start.date.est,]$FIPS[1],]
  County.Number$ob_county_number=NA
  County.Number=County.Number[,c('date','ob_county_number')]
  
  for (i in 1:Odays) {
    County.Number$ob_county_number[i]=nrow(reg.dat.est2v5[reg.dat.est2v5$date==County.Number$date[i],])
  }
  
  pre.end.date=min(reg.dat.est2v5$date)+5

}

#################  Estimation  ####################

{
  
  ### tau_0, or pre-estimation will be performed based on observations from 
  ### Start.date.est to pre.end.date, i.e, March 15 ~ Mar 20.
  
  iter.step=as.numeric(End.date.est-pre.end.date)  ## Iterative steps 
  NumCov=5 ### Number of covariates
  ### EstPhi, EstBeta, for instore the estimated regression parameter
  EstPara<-matrix(data = NA,nrow=iter.step,ncol = (2+NumCov) )
  EstPhi<-EstPara[,1:2]
  EstBeta<-EstPara[,3:(2+NumCov)]
  
  {  
    
    # Estimating EstR[1:tau_0] based on MLE
    
    reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$EstR = 
      reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$incident_cases/
      reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$Lambda
    
    NumCounty=nrow(reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ])
    
    for (i in 2:NumCounty ) {
      
      if ( reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$FIPS[i]==
           reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$FIPS[(i-1)] ){
        
        if ( ( reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$EstR[i]-
               reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$EstR[(i-1)] )>5  ){
          reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$EstR[i]= 
            reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$EstR[(i-1)] + 5
        }else if (( reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$EstR[i]-
                    reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$EstR[(i-1)] ) < -5 ){
          reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$EstR[i]=
            reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$EstR[(i-1)] - 5
        }else if( reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$EstR[i] ==0 ){
          reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$EstR[i]= 
            reg.dat.est2v5[reg.dat.est2v5$date<=pre.end.date, ]$EstR[(i-1)]/2
        }
        
      }
      
    }
    
  } 
  
  save.image('covid_real_0706_v2vari5.RData')
  
  ############### Iterations
  
  source('./functions/covid-function_reorganized_0706.R')
  
  time0=Sys.time()
  for (i in 1:iter.step) {
    inter.result=est.R.iterative.with.BCI(reg.dat.est2v5,dayt=i)
    EstPara[i,]=inter.result$para
    reg.dat.est2v5[reg.dat.est2v5$date <= (pre.end.date+i),]$EstR=
      (inter.result$my_dat_est)$EstR
    cat(i)
  }
  time1=Sys.time()-time0
  time1
  
  EstPara
  
  save.image('covid_real_res_20200706_v2vari5.RData')

}
