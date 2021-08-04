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
library(tools)
library(maps)
library(sf)
library(RColorBrewer)
library(cowplot) 
library(tidyverse)
library(urbnmapr)

# set working directory
setwd('/Users/shij/Desktop/Covid_code_real_data_organize/')
# source R function
source('./functions/covid-function_reorganized_0706.R')

load('./covid_real_res_20200706_v2vari5.RData')

{
  
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


R0summarize<-data.frame(date=c(rep(Start.date.est,Odays)))
R0summarize$R0RedCount<-0
for (i in 1:Odays) {
  R0summarize$date[i]=Start.date.est+i-1
  R0summarize$R0RedCount[i]=nrow(reg.dat.est2v5[reg.dat.est2v5$date==
                                                  R0summarize$date[i] & 
                                                  reg.dat.est2v5$EstR >=1 ,])
}

reg.dat.est.Mar23=reg.dat.est2v5[reg.dat.est2v5$date=="2020-03-23",]

reg.dat.est.Apr07=reg.dat.est2v5[reg.dat.est2v5$date=="2020-04-07",]
reg.dat.est.Apr20=reg.dat.est2v5[reg.dat.est2v5$date=="2020-04-20",]
reg.dat.est.May22=reg.dat.est2v5[reg.dat.est2v5$date=="2020-05-22",]
reg.dat.est.Jun13=reg.dat.est2v5[reg.dat.est2v5$date=="2020-06-13",]
reg.dat.est.Jun29=reg.dat.est2v5[reg.dat.est2v5$date=="2020-06-29",]
reg.dat.est.Jul06=reg.dat.est2v5[reg.dat.est2v5$date=="2020-07-06",]

for (i in 1:nrow(reg.dat.est.Mar23)) {
  idx=sprintf("%05d", as.numeric(reg.dat.est.Mar23$FIPS[i]) )
  reg.dat.est.Mar23$FIPS[i]=as.character(idx)
}


for (i in 1:nrow(reg.dat.est.Apr07)) {
  idx=sprintf("%05d", as.numeric(reg.dat.est.Apr07$FIPS[i]) )
  reg.dat.est.Apr07$FIPS[i]=as.character(idx)
}

for (i in 1:nrow(reg.dat.est.Apr20)) {
  idx=sprintf("%05d", as.numeric(reg.dat.est.Apr20$FIPS[i]) )
  reg.dat.est.Apr20$FIPS[i]=as.character(idx)
}

for (i in 1:nrow(reg.dat.est.May22)) {
  idx=sprintf("%05d", as.numeric(reg.dat.est.May22$FIPS[i]) )
  reg.dat.est.May22$FIPS[i]=as.character(idx)
}

for (i in 1:nrow(reg.dat.est.Jun13)) {
  idx=sprintf("%05d", as.numeric(reg.dat.est.Jun13$FIPS[i]) )
  reg.dat.est.Jun13$FIPS[i]=as.character(idx)
}

for (i in 1:nrow(reg.dat.est.Jun29)) {
  idx=sprintf("%05d", as.numeric(reg.dat.est.Jun29$FIPS[i]) )
  reg.dat.est.Jun29$FIPS[i]=as.character(idx)
}

for (i in 1:nrow(reg.dat.est.Jul06)) {
  idx=sprintf("%05d", as.numeric(reg.dat.est.Jul06$FIPS[i]) )
  reg.dat.est.Jul06$FIPS[i]=as.character(idx)
}


##############################################################
####### First Pic

household_data <- left_join(countydata, counties, by = "county_fips") 
a=household_data
a=a[,c("county_fips","long","lat","state_abbv")]
names(a)[1]="FIPS"
b1=reg.dat.est.Apr07[,c("FIPS","incident_cases","EstR")]
c1=left_join(a,b1,by="FIPS")
c1[is.na(c1$EstR),]$EstR=0

states <- plot_usmap("states", 
                     color = "black",
                     labels = TRUE,
                     fill = alpha(0.01))

names(c1)[1]<-"fips"
counties<-plot_usmap(data=c1,values = "EstR",color = "white",size =0.1,labels = TRUE)

fc<-rev(brewer.pal(11,"RdYlBu"))
fc=fc[-c(1,3,5)]
fc=c("gray96",fc)
fcbvc1<-c(0,0.1,0.5,0.75,0.9,1.1,1.25,1.5,2, max(c1$EstR)) / max(c1$EstR)

p1<-ggplot() +
  counties$layers[[1]] + #counties needs to be on top of states for this to work
  states$layers[[1]] +
  counties$theme + 
  coord_equal() +
  theme(legend.position="none")+
  labs(title = "April 07th, 2020      ")+
  scale_fill_gradientn(colours = fc, values =fcbvc1)

##############################################################
####### Second Pic

b2=reg.dat.est.Apr20[,c("FIPS","incident_cases","EstR")]
c2=left_join(a,b2,by="FIPS")
c2[is.na(c2$EstR),]$EstR=0

names(c2)[1]<-"fips"
counties<-plot_usmap(data=c2,values = "EstR",color = "white",size =0.1,labels = TRUE)

fcbvc2<-c(0,0.1,0.5,0.75,0.9,1.1,1.25,1.5,2, max(c2$EstR)) / max(c2$EstR)

p2<-ggplot() +
  counties$layers[[1]] + #counties needs to be on top of states for this to work
  states$layers[[1]] +
  counties$theme + 
  coord_equal() +
  theme(legend.position="top")+
  #labs(fill = "Instantaneous Reproduction Number         ")+
  labs(fill = "Instantaneous Reproduction Number         \n April 20th, 2020,")+
  #labs(title = "April 20th, 2020      ")+
  scale_fill_gradientn(colours = fc, values =fcbvc2)

##############################################################
####### Third Pic

b3=reg.dat.est.May22[,c("FIPS","incident_cases","EstR")]
c3=left_join(a,b3,by="FIPS")
c3[is.na(c3$EstR),]$EstR=0

names(c3)[1]<-"fips"
counties<-plot_usmap(data=c3,values = "EstR",color = "white",size =0.1,labels = TRUE)

fcbvc3<-c(0,0.1,0.5,0.75,0.9,1.1,1.25,1.5,2, max(c3$EstR)) / max(c3$EstR)

p3<-ggplot() +
  counties$layers[[1]] + #counties needs to be on top of states for this to work
  states$layers[[1]] +
  counties$theme + 
  coord_equal() +
  theme(legend.position="none")+
  #labs(fill = "Instantaneous Reproduction Number         ")+
  #labs(fill = "May 22nd, 2020, Instantaneous Reproduction Number         ")+
  labs(title = "May 22nd, 2020      ")+
  scale_fill_gradientn(colours = fc, values =fcbvc3)

##############################################################
####### Fourth Pic

b4=reg.dat.est.Jun13[,c("FIPS","incident_cases","EstR")]
c4=left_join(a,b4,by="FIPS")
c4[is.na(c4$EstR),]$EstR=0

names(c4)[1]<-"fips"
counties<-plot_usmap(data=c4,values = "EstR",color = "white",size =0.1,labels = TRUE)

fcbvc4<-c( 0,0.1,0.5,0.75,0.9,1.1,1.25,1.5,2,max(c4$EstR)) / max(c4$EstR)

p4<-ggplot() +
  counties$layers[[1]] + #counties needs to be on top of states for this to work
  states$layers[[1]] +
  counties$theme + 
  coord_equal() +
  theme(legend.position="none")+
  #labs(fill = "Instantaneous Reproduction Number         ")+
  labs(title = "June 13th, 2020      ")+
  scale_fill_gradientn(colours = fc, values =fcbvc4)

##############################################################
####### Fifth Pic

b5=reg.dat.est.Jun29[,c("FIPS","incident_cases","EstR")]
c5=left_join(a,b5,by="FIPS")
c5[is.na(c5$EstR),]$EstR=0

names(c5)[1]<-"fips"
counties<-plot_usmap(data=c5,values = "EstR",color = "white",size =0.1,labels = TRUE)

fcbvc5<-c( 0,0.1,0.5,0.75,0.9,1.1,1.25,1.5,2,max(c5$EstR)) / max(c5$EstR)

p5<-ggplot() +
  counties$layers[[1]] + #counties needs to be on top of states for this to work
  states$layers[[1]] +
  counties$theme + 
  coord_equal() +
  theme(legend.position="none")+
  #labs(fill = "Instantaneous Reproduction Number         ")+
  labs(title = "June 29th, 2020      ")+
  scale_fill_gradientn(colours = fc, values =fcbvc5)

##############################################################
####### Sixth Pic

b6=reg.dat.est.Jul06[,c("FIPS","incident_cases","EstR")]
c6=left_join(a,b6,by="FIPS")
c6[is.na(c6$EstR),]$EstR=0

names(c6)[1]<-"fips"
counties<-plot_usmap(data=c6,values = "EstR",color = "white",size =0.1,labels = TRUE)

fcbvc6<-c( 0,0.1,0.5,0.75,0.9,1.1,1.25,1.5,2,max(c6$EstR)) / max(c6$EstR)

p6<-ggplot() +
  counties$layers[[1]] + #counties needs to be on top of states for this to work
  states$layers[[1]] +
  counties$theme + 
  coord_equal() +
  theme(legend.position="none")+
  #labs(fill = "Instantaneous Reproduction Number         ")+
  labs(title = "July 06th, 2020      ")+
  scale_fill_gradientn(colours = fc, values =fcbvc6)

################################################################
################Seventh Pic


b7=reg.dat.est.Mar23[,c("FIPS","incident_cases","EstR")]
c7=left_join(a,b7,by="FIPS")
c7[is.na(c7$EstR),]$EstR=0

names(c7)[1]<-"fips"
counties<-plot_usmap(data=c7,values = "EstR",color = "white",size =0.1,labels = TRUE)

fcbvc7<-c(0,0.1,0.5,0.75,0.9,1.1,1.25,1.5,2, max(c7$EstR)) / max(c7$EstR)

p7<-ggplot() +
  counties$layers[[1]] + #counties needs to be on top of states for this to work
  states$layers[[1]] +
  counties$theme + 
  coord_equal() +
  theme(legend.position="none")+
  #labs(fill = "Instantaneous Reproduction Number         ")+
  labs(fill = "March 23rd, 2020      ")+
  #labs(title = "April 20th, 2020      ")+
  scale_fill_gradientn(colours = fc, values =fcbvc7)


###########################################################


###########Combine four pics

ggdraw() +
  draw_plot(p1, 0, .66, .5, .33) +
  draw_plot(p2, .5, .66, .5, .34) +
  draw_plot(p3, 0, 0.33, .5, .33) +
  draw_plot(p4, .5, .33, .5, .33)+
  draw_plot(p5, 0, 0, .5, .33)+
  draw_plot(p6, .5, 0, .5, .33)






