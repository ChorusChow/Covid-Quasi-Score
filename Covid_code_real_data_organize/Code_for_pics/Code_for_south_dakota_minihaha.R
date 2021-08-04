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

load('./covid_real_res_20200706_v2vari5.RData')

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


Minnehaha.R0=reg.dat.est2v5[reg.dat.est2v5$STATE_NAME=="South Dakota",]

####################################################################################


EstParaBoot=list()
EstParaEnds=matrix(data = NA, nrow=69,ncol=7)

for (i in 1:69) {
  EstParaBoot[[i]]=readRDS(file=paste0('derived/result_ellboot45/covid_real_res_bs_20200706_v2vari5_ellboot45_i', i, '.RDS'))
  
  EstParaEnds[i,]=EstParaBoot[[i]][39,]
}

repboot=200

set.seed(1115)
index2v5ellboot45=sample(1:69,size=repboot,replace = TRUE)

ResPara2v5ellboot45=matrix(data = NA, nrow = repboot,ncol = 7)

for (i  in 1:repboot) {
  ResPara2v5ellboot45[i,]=EstParaEnds[index2v5ellboot45[i],]
}

#########################################################################


BootstrapEstR<-matrix(data=NA,nrow=repboot,ncol = nrow(Minnehaha.R0))
for (jBootrep in 1:200) {
  BootstrapEstR[jBootrep,1]= meanprior
    #exp( res1rep1EstPhiBoot[1,jBootrep]+
    #                               res1rep1EstPhiBoot[2,jBootrep]*log(3)+
    #                               res1rep1Z[1,]%*%res1rep1EstBetaBoot[,jBootrep] )
  for (Rt in 2:53) {
    BootstrapEstR[jBootrep,Rt]= exp( ResPara2v5ellboot45[jBootrep,1]+
                                       ResPara2v5ellboot45[jBootrep,2]*log( BootstrapEstR[jBootrep,(Rt-1)] )+
                                       Minnehaha.R0$PopO64[Rt]*ResPara2v5ellboot45[jBootrep,3]+
                                       Minnehaha.R0$popdens_s[Rt]*ResPara2v5ellboot45[jBootrep,4]+
                                       Minnehaha.R0$diabetes_s[Rt]*ResPara2v5ellboot45[jBootrep,5]+
                                       Minnehaha.R0$visit_mean[Rt]*ResPara2v5ellboot45[jBootrep,6]+
                                       Minnehaha.R0$temp_mean[Rt]*ResPara2v5ellboot45[jBootrep,7])
  }
}

for (Rt in 1:nrow(Minnehaha.R0)) {
  a=sort(BootstrapEstR[,Rt])
  Minnehaha.R0$BCIupper[Rt] =a[190]
  Minnehaha.R0$BCIlower[Rt] =a[11]
}


Minnehaha.R0=Minnehaha.R0[,c("date","EstR","BCIupper","BCIlower")]

Minnehaha.R0$date=as.Date(Minnehaha.R0$date)
Minnehaha.R0=Minnehaha.R0[Minnehaha.R0$date<=(as.Date("2020-04-02") + 52), ]

inst.reproduction.nb=data.frame(day = as.Date("2020-04-02") + 0:52)
inst.reproduction.nb$Est.R0=Minnehaha.R0$EstR
inst.reproduction.nb$Est.R0.band.upp=Minnehaha.R0$BCIupper
inst.reproduction.nb$Est.R0.band.low=Minnehaha.R0$BCIlower

# ggplot(inst.reproduction.nb, aes(x=day)) + 
#   #geom_line(aes(y=oracle.R0), colour="blue") + 
#   geom_line(aes(y=Est.R0), colour="purple",lty=5)+
#   geom_line(aes(y=1), colour="black",lty=2)+
#   #labs(
#     #title = "",
#     #subtitle = "",
#     #tag = "",
#     #caption = "",
#   #  x = "Time",
#   #  y = "instantaneous reproduction number"
#   #)+
#   annotate(geom = "text",
#            x = as.Date("2020-05-30"), 
#            y = 2.5,
#            label = "Control Group")+
#   geom_ribbon(aes(ymin=Est.R0.band.low, ymax=Est.R0.band.upp), alpha=0.2)

inst.reproduction.nb$type.R0=c("Minnehaha, South Dakota")

ggplot(inst.reproduction.nb, aes(day)) + 
  geom_line(aes(y=Est.R0, colour="Daily R0"),colour="violetred") + 
  # geom_line(aes(y=Est.R0), colour="purple",lty=5)+
  geom_line(aes(y=1), colour="black",lty=2)+
  geom_vline(aes(xintercept=as.Date("2020-04-12")), colour="black",lty=4)+
  geom_vline(aes(xintercept=as.Date("2020-05-10")), colour="black",lty=4)+
  labs(
    #title = "",
    #subtitle = "",
    #tag = "",
    #caption = "",
    #x = "Time",
    y = "instantaneous reproduction number"
  )+
  annotate(geom = "text",
           x = as.Date("2020-04-18"), 
           y = 2.25,
           label = "April 12, First \n reported COVID-19 \n outbreak in Smith- \n field Foods, Inc.",
           size=3.45)+
  annotate(geom = "text",
           x = as.Date("2020-05-04"), 
           y = 1.75,
           label = "May 10, Second \n reported COVID-19 \n outbreak relates to \n testing site at Wash- \n ington High School \n for Smithfield emplo- \n yee and their family",
           size=3.45)+
  annotate(geom = "text",
           x = as.Date("2020-05-19"), 
           y = 0.69,
           label = "90% Interval Estimator",
           size=3)+
  geom_ribbon(aes(ymin=Est.R0.band.low, ymax=Est.R0.band.upp), alpha=0.2)+
  facet_grid(type.R0 ~ .)+
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0.75,0.9),
    axis.title.x=element_blank()
    # legend.background = 
  ) +
  scale_color_manual(values=c("dodgerblue", "orchid"))+
  scale_linetype_manual(values=c("twodash", "solid"))+
  ylim(0.45, 2.5)
  #scale_x_continuous(name='Time', breaks=c(30, 60, 90, 120), 
   #                  labels=c(30, 60, 90, 120), limits=c(0,120) ) 



