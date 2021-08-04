
########!My function############################################################

Lambda.gen.est <-function(my.dat.est = reg.dat.est2v5, Omega.nCov = GT.ncov, meanprior){
  
  my.dat.summarize.est <- my.dat.est%>% group_by(FIPS)%>%
    summarize(county = unique(county), state= unique(STATE_NAME), 
              max_cases = max(cases,na.rm = TRUE),ob_days= sum(cases>0))
  
  my.County.Number=nrow(my.dat.summarize.est)
  
  for (i in 1: my.County.Number ) {
    my.dat.est[my.dat.est$FIPS == my.dat.summarize.est$FIPS[i], ]$Lambda[1] =
      my.dat.est[my.dat.est$FIPS == my.dat.summarize.est$FIPS[i], ]$incident_cases[1]/meanprior
    for (j in 2:my.dat.summarize.est$ob_days[i]) {
      my.dat.est[my.dat.est$FIPS == my.dat.summarize.est$FIPS[i], ]$Lambda[j] =
        my.dat.est[my.dat.est$FIPS == my.dat.summarize.est$FIPS[i], ]$incident_cases[1:(j-1)]%*%
        Omega.nCov[(j-1):1]
    }
  }
  return(my.dat.est)
}

greater.or.na <- function(myobj, myconst, tocat){
  if(is.na(myobj)) cat(myobj, tocat, '\n')
  return(is.na(myobj) | myobj > myconst)
}

est.R.iterative.with.BCI <- function( my.dat.est = reg.dat.est2v5, dayt, meanprior, 
                                      Omega.nCov = GT.ncov, tau_0=5, NumCov=5
                                       #mystation = mystation, myst=myst,
                                       #opendate = as.Date('4000-01-01','%Y-%m-%d'), openpct = 0)
                                      ){
############################
 # my.dat.est = reg.dat.est2v5
 # Omega.nCov: corresponding to GT.ncov?
 # meanprior: I0, start reproduction number
 # my.dat.est$incident_cases = reg.dat.est2v5$incident_cases
############################
  my.dat.est=my.dat.est[my.dat.est$date<= (pre.end.date+ dayt),]
               
  
   my.dat.summarize.est <- my.dat.est%>% group_by(FIPS)%>%
     summarize(county = unique(county), state= unique(STATE_NAME), 
               max_cases = max(cases,na.rm = TRUE),ob_days= sum(cases>0),
               min_ob_day=min(date),max_ob_day=max(date))
  
   my.dat.est$label.ob1.day=0
  ### label each county's first observation day
   for (i in 1:nrow(my.dat.summarize.est)) {
     my.dat.est[my.dat.est$FIPS==my.dat.summarize.est$FIPS[i] & 
                    my.dat.est$date==my.dat.summarize.est$min_ob_day[i], ]$label.ob1.day=1
   }
   
   
   ### tau_0, or pre-estimation will be performed based on observations from 
  ### Start.date.est to pre.end.date, i.e, March 15 ~ Mar 20.
  
  ##################################################
  
  ## Covariates
  ## Define intermediate variables, YTilde, BetaTilde, WTilde, ZYTilde, barZ
  ### barZ[t], arverage of Z[i,]'s from date 1 to date t.
  ###
  my.dat.est$YTilde=NA
  my.dat.est[!is.na(my.dat.est$EstR),]$YTilde=
    log(my.dat.est[!is.na(my.dat.est$EstR),]$EstR)
  
  #BetaTilde<-matrix(data=NA, nrow=iter.step, ncol=NumCov)
  
  ### barZ coordinate corresponds to "PopO64","popdens_s","diabetes_s",
  ### "visit_mean","temp_mean" in order
  #barZ<-matrix(NA,nrow = iter.step,ncol = NumCov)
  #WTilde<-list()
  ### Calculate barZ and WTilde, weight  \Big( \sum_{i=1}^t (Z_i-\bar{Z})(Z_i-\bar{Z})^T \Big)^{-1}
  #ZZCov=my.dat.est
  
  
    ZCov=my.dat.est[my.dat.est$date<=(pre.end.date+ dayt -1) & my.dat.est$label.ob1.day != 1, ]
    #ZCov=ZCov[!(is.na(ZCov$dry_bulb_temp) | is.na(ZCov$rel_humid) | is.na(ZCov$visit_mean) ), ]
    ZCov=ZCov[,c("PopO64","popdens_s","diabetes_s","visit_mean","temp_mean")]
    ZCov=as.matrix(ZCov)
    barZ=colMeans(ZCov)
    WTilde<-matrix(NA,nrow = NumCov,ncol = NumCov)
    WTilde= solve( t(ZCov-as.matrix(rep(1,nrow(ZCov)))%*%barZ)%*%
                 (ZCov-as.matrix(rep(1,nrow(ZCov)))%*%barZ)  )
    
    ZYCov=my.dat.est[my.dat.est$date<=(pre.end.date+ dayt -1), ]
    #ZYCov=ZYCov[!(is.na(ZYCov$dry_bulb_temp) | is.na(ZYCov$rel_humid) | is.na(ZYCov$visit_mean) ), ]
    ### Define the profile likelihood
    ellk<- function(phi){
      phi_0=phi[1]
      phi_1=phi[2]
      
      ### Calculate ZYTilde, which is stored in ZYCov[,c("PopO64",
      ###                                 "popdens_s","diabetes_s","visit_mean","temp_mean")]
      for (i in 1:nrow(ZYCov)) {
        if (ZYCov$label.ob1.day[i]==1){
          ZYCov[i,c("PopO64","popdens_s","diabetes_s","visit_mean","temp_mean")]=0
        }else if(ZYCov$label.ob1.day[i]==0){
          ZYCov[i,c("PopO64","popdens_s","diabetes_s","visit_mean","temp_mean")]=  
          (ZYCov[i,c("PopO64","popdens_s",
                       "diabetes_s","visit_mean","temp_mean")] - barZ  )*
            ( log( ZYCov$EstR[i] )-phi_1*log( ZYCov$EstR[(i-1)] ) - phi_0 )
        }
      }
      
      ###Calculate BetaTilde[k,]
      BetaTilde=colSums(
          as.matrix(ZYCov[,c("PopO64","popdens_s",
                             "diabetes_s","visit_mean","temp_mean")]) %*% WTilde
      )
      
      ### Calculate YTilde in my.dat.est
      NumCounty=nrow(my.dat.est)
      for (i in 2:NumCounty ) {
        if (my.dat.est$date[i]> pre.end.date & my.dat.est$label.ob1.day[i]==0 ){
          my.dat.est$YTilde[i]= phi_0+ phi_1 * log( my.dat.est$EstR[(i-1)]  ) +
            as.matrix( my.dat.est[i,c("PopO64","popdens_s","diabetes_s",
                                      "visit_mean","temp_mean")] ) %*% BetaTilde
        }
      }
      
      ### Calculate Likelihood in my.dat.est
      Likelihood.est=my.dat.est[my.dat.est$date>pre.end.date & my.dat.est$date<= (pre.end.date+ dayt)
                                & my.dat.est$label.ob1.day != 1,] 
      Likelihood.est$prolikelihood= (Likelihood.est$incident_cases) * (Likelihood.est$YTilde ) - 
        (exp(Likelihood.est$YTilde))*(Likelihood.est$Lambda  )    
      
      result= - sum(Likelihood.est$prolikelihood)
      
      return(result)
    }     
    
    EstPhidayt=nlminb(c(0.15,0.7),ellk,lower = c((- Inf),0.01),upper = c(Inf,0.99))$par
    
    ### Calculate ZYHat, which is stored in ZYCov[,c("PopO64",
    ###                                 "popdens_s","diabetes_s","visit_mean","temp_mean")]
    ZYCov=my.dat.est[my.dat.est$date<=(pre.end.date+ dayt -1), ]
  
    for (i in 1:nrow(ZYCov)) {
      if (ZYCov$label.ob1.day[i]==1){
        ZYCov[i,c("PopO64","popdens_s","diabetes_s","visit_mean","temp_mean")]=0
      }else if(ZYCov$label.ob1.day[i]==0){
        ZYCov[i,c("PopO64","popdens_s","diabetes_s","visit_mean","temp_mean")]=  
          (ZYCov[i,c("PopO64","popdens_s",
                     "diabetes_s","visit_mean","temp_mean")] - barZ  )*
          ( log( ZYCov$EstR[i] )-EstPhidayt[2]*log( ZYCov$EstR[(i-1)] ) - EstPhidayt[1] )
      }
    }
    
    ###Calculate EstBeta[k,]
    EstBetadayt=colSums(
      as.matrix(ZYCov[,c("PopO64","popdens_s",
                         "diabetes_s","visit_mean","temp_mean")]) %*% WTilde
    )
    
    EstParadayt=c()
    EstParadayt[1:2]=EstPhidayt
    EstParadayt[3:(2+NumCov)]=EstBetadayt
    
    NumCounty= nrow(my.dat.est)
    for (i in 1:NumCounty) {
      
      if (is.na(my.dat.est$EstR[i]) ){
        my.dat.est$EstR[i]=
          exp( EstPhidayt[1]+
               EstPhidayt[2]*log(my.dat.est$EstR[(i-1)]) +
               as.matrix(my.dat.est[i,c("PopO64", "popdens_s","diabetes_s",
                                        "visit_mean","temp_mean")]) %*%EstBetadayt ) 
           if ( ( my.dat.est$EstR[i]- my.dat.est$EstR[(i-1)] )>5  ){
              my.dat.est$EstR[i]= my.dat.est$EstR[(i-1)] + 5
            }else if (( my.dat.est$EstR[i]- my.dat.est$EstR[(i-1)] ) < -5 ){
              my.dat.est$EstR[i]= my.dat.est$EstR[(i-1)] - 5
            }else if( my.dat.est$EstR[i] ==0 ){
              my.dat.est$EstR[i]=  my.dat.est$EstR[(i-1)]/2
            }
          
        }
        
    }
    
    my_fun_return_list=list("para" = EstParadayt, "my_dat_est" = my.dat.est )
    return(my_fun_return_list)
    
}


estR.iterative.Boot <- function( my.dat.est = reg.dat.boot2v5, dayt, meanprior, 
                                      Omega.nCov = GT.ncov, tau_0=5, NumCov=5){
  
  ######################
  my.dat.est=my.dat.est[my.dat.est$date<= (pre.end.boot+ dayt),]
  
  
  my.dat.summarize.est <- my.dat.est%>% group_by(FIPS)%>%
    summarize(county = unique(county), state= unique(STATE_NAME), 
              max_cases = max(cases,na.rm = TRUE),ob_days= sum(cases>0),
              min_ob_day=min(date),max_ob_day=max(date))
  
  my.dat.est$label.ob1.day=0
  ### label each county's first observation day
  for (i in 1:nrow(my.dat.summarize.est)) {
    my.dat.est[my.dat.est$FIPS==my.dat.summarize.est$FIPS[i] & 
                 my.dat.est$date==my.dat.summarize.est$min_ob_day[i], ]$label.ob1.day=1
  }
  
  
  ### tau_0, or pre-estimation will be performed based on observations from 
  ### Start.date.boot to pre.end.boot, i.e, March 15 ~ Mar 20.
  
  ##################################################
  
  ## Covariates
  ## Define intermediate variables, YTilde, BetaTilde, WTilde, ZYTilde, barZ
  ### barZ[t], arverage of Z[i,]'s from date 1 to date t.
  ###
  my.dat.est$YTilde=NA
  my.dat.est[!is.na(my.dat.est$EstRboot),]$YTilde=
    log(my.dat.est[!is.na(my.dat.est$EstRboot),]$EstRboot)
  
  #BetaTilde<-matrix(data=NA, nrow=iter.step, ncol=NumCov)
  
  ### barZ coordinate corresponds to "PopO64","popdens_s","diabetes_s",
  ### "visit_mean","temp_mean" in order
  #barZ<-matrix(NA,nrow = iter.step,ncol = NumCov)
  #WTilde<-list()
  ### Calculate barZ and WTilde, weight  \Big( \sum_{i=1}^t (Z_i-\bar{Z})(Z_i-\bar{Z})^T \Big)^{-1}
  #ZZCov=my.dat.est
  
  
  ZCov=my.dat.est[my.dat.est$date<=(pre.end.boot+ dayt -1) & my.dat.est$label.ob1.day != 1, ]
  #ZCov=ZCov[!(is.na(ZCov$dry_bulb_temp) | is.na(ZCov$rel_humid) | is.na(ZCov$visit_mean) ), ]
  ZCov=ZCov[,c("PopO64","popdens_s","diabetes_s","visit_mean","temp_mean")]
  ZCov=as.matrix(ZCov)
  barZ=colMeans(ZCov)
  WTilde<-matrix(NA,nrow = NumCov,ncol = NumCov)
  WTilde= solve( t(ZCov-as.matrix(rep(1,nrow(ZCov)))%*%barZ)%*%
                   (ZCov-as.matrix(rep(1,nrow(ZCov)))%*%barZ)  )
  
  ZYCov=my.dat.est[my.dat.est$date<=(pre.end.boot+ dayt -1), ]
  #ZYCov=ZYCov[!(is.na(ZYCov$dry_bulb_temp) | is.na(ZYCov$rel_humid) | is.na(ZYCov$visit_mean) ), ]
  ### Define the profile likelihood
  ellk<- function(phi){
    phi_0=phi[1]
    phi_1=phi[2]
    
    ### Calculate ZYTilde, which is stored in ZYCov[,c("PopO64",
    ###                                 "popdens_s","diabetes_s","visit_mean","temp_mean")]
    for (i in 1:nrow(ZYCov)) {
      if (ZYCov$label.ob1.day[i]==1){
        ZYCov[i,c("PopO64","popdens_s","diabetes_s","visit_mean","temp_mean")]=0
      }else if(ZYCov$label.ob1.day[i]==0){
        ZYCov[i,c("PopO64","popdens_s","diabetes_s","visit_mean","temp_mean")]=  
          (ZYCov[i,c("PopO64","popdens_s",
                     "diabetes_s","visit_mean","temp_mean")] - barZ  )*
          ( log( ZYCov$EstRboot[i] )-phi_1*log( ZYCov$EstRboot[(i-1)] ) - phi_0 )
      }
    }
    
    ###Calculate BetaTilde[k,]
    BetaTilde=colSums(
      as.matrix(ZYCov[,c("PopO64","popdens_s",
                         "diabetes_s","visit_mean","temp_mean")]) %*% WTilde
    )
    
    ### Calculate YTilde in my.dat.est
    NumCounty=nrow(my.dat.est)
    for (i in 2:NumCounty ) {
      if (my.dat.est$date[i]> pre.end.boot & my.dat.est$label.ob1.day[i]==0 ){
        my.dat.est$YTilde[i]= phi_0+ phi_1 * log( my.dat.est$EstRboot[(i-1)]  ) +
          as.matrix( my.dat.est[i,c("PopO64","popdens_s","diabetes_s",
                                    "visit_mean","temp_mean")] ) %*% BetaTilde
      }
    }
    
    ### Calculate Likelihood in my.dat.est
    Likelihood.est=my.dat.est[my.dat.est$date>pre.end.boot & my.dat.est$date<= (pre.end.boot+ dayt)
                              & my.dat.est$label.ob1.day != 1,] 
    Likelihood.est$prolikelihood= (Likelihood.est$incident_cases) * (Likelihood.est$YTilde ) - 
      (exp(Likelihood.est$YTilde))*(Likelihood.est$Lambda  )    
    
    result= - sum(Likelihood.est$prolikelihood)
    
    return(result)
  }     
  
  EstPhidayt=nlminb(c(0.15,0.7),ellk,lower = c((- Inf),0.01),upper = c(Inf,0.99))$par
  
  ### Calculate ZYHat, which is stored in ZYCov[,c("PopO64",
  ###                                 "popdens_s","diabetes_s","visit_mean","temp_mean")]
  ZYCov=my.dat.est[my.dat.est$date<=(pre.end.boot+ dayt -1), ]
  
  for (i in 1:nrow(ZYCov)) {
    if (ZYCov$label.ob1.day[i]==1){
      ZYCov[i,c("PopO64","popdens_s","diabetes_s","visit_mean","temp_mean")]=0
    }else if(ZYCov$label.ob1.day[i]==0){
      ZYCov[i,c("PopO64","popdens_s","diabetes_s","visit_mean","temp_mean")]=  
        (ZYCov[i,c("PopO64","popdens_s",
                   "diabetes_s","visit_mean","temp_mean")] - barZ  )*
        ( log( ZYCov$EstRboot[i] )-EstPhidayt[2]*log( ZYCov$EstRboot[(i-1)] ) - EstPhidayt[1] )
    }
  }
  
  ###Calculate EstBeta[k,]
  EstBetadayt=colSums(
    as.matrix(ZYCov[,c("PopO64","popdens_s",
                       "diabetes_s","visit_mean","temp_mean")]) %*% WTilde
  )
  
  EstParadayt=c()
  EstParadayt[1:2]=EstPhidayt
  EstParadayt[3:(2+NumCov)]=EstBetadayt
  
  NumCounty= nrow(my.dat.est)
  for (i in 1:NumCounty) {
    
    if (is.na(my.dat.est$EstRboot[i]) ){
      my.dat.est$EstRboot[i]=
        exp( EstPhidayt[1]+
               EstPhidayt[2]*log(my.dat.est$EstRboot[(i-1)]) +
               as.matrix(my.dat.est[i,c("PopO64", "popdens_s","diabetes_s",
                                        "visit_mean","temp_mean")]) %*%EstBetadayt ) 
      if ( ( my.dat.est$EstRboot[i]- my.dat.est$EstRboot[(i-1)] )>5  ){
        my.dat.est$EstRboot[i]= my.dat.est$EstRboot[(i-1)] + 5
      }else if (( my.dat.est$EstRboot[i]- my.dat.est$EstRboot[(i-1)] ) < -5 ){
        my.dat.est$EstRboot[i]= my.dat.est$EstRboot[(i-1)] - 5
      }else if( my.dat.est$EstRboot[i] ==0 ){
        my.dat.est$EstRboot[i]=  my.dat.est$EstRboot[(i-1)]/2
      }
      
    }
    
  }
  
  my_fun_return_list=list("para" = EstParadayt, "my_dat_est" = my.dat.est )
  return(my_fun_return_list)
  
}
  


