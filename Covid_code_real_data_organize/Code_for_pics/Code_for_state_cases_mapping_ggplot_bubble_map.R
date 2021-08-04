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
library(devtools)
library(ggrepel)
library(grid)

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


reg.dat.est.Apr07=reg.dat.est2v5[reg.dat.est2v5$date=="2020-04-07",]
reg.dat.est.Apr20=reg.dat.est2v5[reg.dat.est2v5$date=="2020-04-20",]
reg.dat.est.May22=reg.dat.est2v5[reg.dat.est2v5$date=="2020-05-22",]
reg.dat.est.Jun13=reg.dat.est2v5[reg.dat.est2v5$date=="2020-06-13",]
reg.dat.est.Jun29=reg.dat.est2v5[reg.dat.est2v5$date=="2020-06-29",]
reg.dat.est.Jul06=reg.dat.est2v5[reg.dat.est2v5$date=="2020-07-06",]


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

############# add 0 to county.fips if county.fips is only 4-digit
my.county.fips=county.fips
for (i in 1:nrow(my.county.fips)) {
  idx=sprintf("%05d", as.numeric(my.county.fips$fips[i]) )
  my.county.fips$fips[i]=as.character(idx)
}

##################################
##################################   First pic
household_data <- left_join(countydata, counties, by = "county_fips") 
a=household_data
a=a[,c("county_fips","long","lat","state_abbv")]
names(a)[1]="FIPS"
b1=reg.dat.est.Apr07[,c("FIPS","incident_cases","EstR")]
c1=left_join(a,b1,by="FIPS")
c1[is.na(c1$EstR),]$EstR=0
names(c1)[1]<-"fips"


c1=c1[c1$state_abbv %in% c("PA","ME","NY","NJ","MA","VT","RI","CT","NH"),]
c1summarize<- c1%>% group_by(fips)%>%
  summarize(state= unique(state_abbv), 
            long = mean(long,na.rm = TRUE),lat= mean(lat,na.rm = TRUE),
            incident_cases  =  unique(incident_cases))
c1summarize=c1summarize[!is.na(c1summarize$incident_cases),]
c1summarize<-c1summarize[c(3,4,1,2,5)]


# cities_t <- usmap_transform(citypop)
# cities_t1 <-cities_t[cities_t$abbr %in% c("PA","ME","NY","NJ","MA","VT","RI","CT","NH"),]

c1summarize <- usmap_transform(c1summarize)
c1summarize$county_name<-NA
most_incident_cases_city1=rev(order(c1summarize$incident_cases))[1:9]
city.for.label1=c1summarize[most_incident_cases_city1,]

############  Match city/county name

city.name.for.label1=my.county.fips[my.county.fips$fips %in% city.for.label1$fips,]
for (i in 1:nrow(city.name.for.label1)) {
  c1summarize[c1summarize$fips==city.name.for.label1$fips[i],]$county_name=
    city.name.for.label1$polyname[i]
}

############  Manual label
c1summarize$county_name[21]<-"NJ,Bergen,600+"
c1summarize$county_name[59]<-"NY,Westchester,600+"
c1summarize$county_name[28]<-"NJ,Hudson,500+"
c1summarize$county_name[41]<-"NY,Bronx,1000+"
c1summarize$county_name[44]<-"NY,Kings,1200+"
c1summarize$county_name[46]<-"NY,Nassau,1100+"
c1summarize$county_name[47]<-"NY,New York,600+"
c1summarize$county_name[53]<-"NY,Queens,1500+"
c1summarize$county_name[56]<-"NY,Suffolk,1000+"

####################################################


p1<-plot_usmap(fill="gray", alpha=0.3,include = c("PA","ME","NY","NJ","MA","VT","RI","CT","NH"),labels = TRUE) +
  ggrepel::geom_label_repel(data = c1summarize,
                           aes(x = long.1, y = lat.1, label = county_name),
                           size = 3, alpha = 0.8,
                           label.r = unit(0.5, "lines"), label.size = 0.5,
                           segment.color = "red", segment.size = 1,
                           seed = 1002) +
  geom_point(data = c1summarize, aes(x = long.1, y = lat.1, size = incident_cases),
             color = "violetred", alpha = 0.25) +
  scale_size_continuous(range=c(1,14))+
  labs(title = "Incident Cases \n Reported, April 07th",
       #subtitle = "Reported, April 07th",
       size = "Magnitude") +
  #theme(legend.position = "right")+
  theme(
    legend.position = "none",
    #text = element_text(color = "#22211d"),
    #plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    #panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    #legend.background = element_rect(fill = "#f5f5f2", color = NA),
    #plot.title = element_text(size= 16, hjust=0.1, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
  )

##################################
##################################   Third pic

b3=reg.dat.est.May22[,c("FIPS","incident_cases","EstR")]
c3=left_join(a,b3,by="FIPS")
c3[is.na(c3$EstR),]$EstR=0

names(c3)[1]<-"fips"

c3=c3[c3$state_abbv %in% c("PA","ME","NY","NJ","MA","VT","RI","CT","NH"),]
c3summarize<- c3%>% group_by(fips)%>%
  summarize(state= unique(state_abbv), 
            long = mean(long,na.rm = TRUE),lat= mean(lat,na.rm = TRUE),
            incident_cases  =  unique(incident_cases))
c3summarize=c3summarize[!is.na(c3summarize$incident_cases),]
c3summarize<-c3summarize[c(3,4,1,2,5)]

c3summarize <- usmap_transform(c3summarize)
c3summarize$county_name<-NA
most_incident_cases_city3=rev(order(c3summarize$incident_cases))[1:9]
city.for.label3=c3summarize[most_incident_cases_city3,]

############  Match city/county name

city.name.for.label3=my.county.fips[my.county.fips$fips %in% city.for.label3$fips,]
for (i in 1:nrow(city.name.for.label3)) {
  c3summarize[c3summarize$fips==city.name.for.label3$fips[i],]$county_name=
    city.name.for.label3$polyname[i]
}

############  Manual label
c3summarize$county_name[2]<-"CT,Hartford,120+"
c3summarize$county_name[11]<-"MA,Essex,150+"
c3summarize$county_name[13]<-"MA,Middlesex,180+"
c3summarize$county_name[16]<-"MA,Suffolk,120+"
c3summarize$county_name[17]<-"MA,Worcester,160+"
c3summarize$county_name[40]<-"NY,Bronx,150+"
c3summarize$county_name[43]<-"NY,Kings,280+"
c3summarize$county_name[52]<-"NY,Queens,220+"
c3summarize$county_name[74]<-"PA,Philadelphia,170+"


p2<-plot_usmap(fill="gray", alpha=0.3,include = c("PA","ME","NY","NJ","MA","VT","RI","CT","NH"),labels = TRUE) +
  ggrepel::geom_label_repel(data = c3summarize,
                            aes(x = long.1, y = lat.1, label = county_name),
                            size = 7.5, alpha = 0.8,
                            label.r = unit(0.5, "lines"), label.size = 0.5,
                            segment.color = "red", segment.size = 1,
                            seed = 1002) +
  geom_point(data = c3summarize, aes(x = long.1, y = lat.1, size = incident_cases),
             color = "violetred", alpha = 0.25) +
  scale_size_continuous(range=c(1,22.5))+
  labs(title = "Reported Cases, May 22nd",
       size = "Incident Cases") +
  theme(plot.title = element_text(size=19))+
  #theme(legend.position = "right")+
  theme(
    legend.position = "top",
    #text = element_text(color = "#22211d"),
    #plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    #panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    #legend.background = element_rect(fill = "#f5f5f2", color = NA),
    #plot.title = element_text(size= 16, hjust=0.1, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
  )


##################################
##################################   Forth pic

b3=reg.dat.est.May22[,c("FIPS","incident_cases","EstR")]
c3=left_join(a,b3,by="FIPS")
c3[is.na(c3$EstR),]$EstR=0

names(c3)[1]<-"fips"

c4=c3[c3$state_abbv %in% c("NC","TN","MS","SC","AL","GA","FL"),]
c4summarize<- c4%>% group_by(fips)%>%
  summarize(state= unique(state_abbv), 
            long = mean(long,na.rm = TRUE),lat= mean(lat,na.rm = TRUE),
            incident_cases  =  unique(incident_cases))
c4summarize=c4summarize[!is.na(c4summarize$incident_cases),]
c4summarize<-c4summarize[c(3,4,1,2,5)]

c4summarize <- usmap_transform(c4summarize)
c4summarize$county_name<-NA
most_incident_cases_city4=rev(order(c4summarize$incident_cases))[1:11]
city.for.label4=c4summarize[most_incident_cases_city4,]

############  Match city/county name

city.name.for.label4=my.county.fips[my.county.fips$fips %in% city.for.label4$fips,]
for (i in 1:nrow(city.name.for.label4)) {
  c4summarize[c4summarize$fips==city.name.for.label4$fips[i],]$county_name=
    city.name.for.label4$polyname[i]
}

############  Manual label
c4summarize$county_name[8]<-"AL,Montgomery,50+"
c4summarize$county_name[13]<-"FL,broward,60+"
c4summarize$county_name[23]<-"FL,Miami-dade,160+"
c4summarize$county_name[26]<-"FL,Palm Beach,100+"
c4summarize$county_name[41]<-"GA,De Kalb,50+"
c4summarize$county_name[46]<-"GA,Gwinnett,60+"
#c4summarize$county_name[60]<-"MS,Rankin,320+"
c4summarize$county_name[60]<-NA
c4summarize$county_name[74]<-"NC,Mecklenburg,70+"
#c4summarize$county_name[82]<-"NC,Wilson,270+"
c4summarize$county_name[82]<-NA
c4summarize$county_name[92]<-"TN,Davidson,80+"
c4summarize$county_name[98]<-"TN,Shelby,70+"


p4<-plot_usmap(fill="gray", alpha=0.3,include = c("NC","TN","MS","SC","AL","GA","FL"),labels = TRUE) +
  ggrepel::geom_label_repel(data = c4summarize,
                            aes(x = long.1, y = lat.1, label = county_name),
                            size = 7.5, alpha = 0.8,
                            label.r = unit(0.5, "lines"), label.size = 0.5,
                            segment.color = "red", segment.size = 1,
                            seed = 1002) +
  geom_point(data = c4summarize, aes(x = long.1, y = lat.1, size = incident_cases),
             color = "violetred", alpha = 0.25) +
  scale_size_continuous(range=c(1,21))+
  #labs(
    #title = "Reported, May 22nd",
  #     size = "Incident Cases") +
  #theme(legend.position = "right")+
  theme(
    legend.position = "none",
    #text = element_text(color = "#22211d"),
    #plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    #panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    #legend.background = element_rect(fill = "#f5f5f2", color = NA),
    #plot.title = element_text(size= 16, hjust=0.1, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
  )


##################################
##################################   Fifth pic

b5=reg.dat.est.Jun29[,c("FIPS","incident_cases","EstR")]
c5=left_join(a,b5,by="FIPS")
c5[is.na(c5$EstR),]$EstR=0

names(c5)[1]<-"fips"

c5=c5[c5$state_abbv %in% c("PA","ME","NY","NJ","MA","VT","RI","CT","NH"),]
c5summarize<- c5%>% group_by(fips)%>%
  summarize(state= unique(state_abbv), 
            long = mean(long,na.rm = TRUE),lat= mean(lat,na.rm = TRUE),
            incident_cases  =  unique(incident_cases))
c5summarize=c5summarize[!is.na(c5summarize$incident_cases),]
c5summarize<-c5summarize[c(3,4,1,2,5)]

c5summarize <- usmap_transform(c5summarize)
c5summarize$county_name<-NA
most_incident_cases_city5=rev(order(c5summarize$incident_cases))[1:9]
city.for.label5=c5summarize[most_incident_cases_city5,]

############  Match city/county name

city.name.for.label5=my.county.fips[my.county.fips$fips %in% city.for.label5$fips,]
for (i in 1:nrow(city.name.for.label5)) {
  c5summarize[c5summarize$fips==city.name.for.label5$fips[i],]$county_name=
    city.name.for.label5$polyname[i]
}

############  Manual label
c5summarize$county_name[8]<-"MA,Middlesex,50+"
c5summarize$county_name[15]<-"NJ,Bergen,50+"
c5summarize$county_name[31]<-"NY,Bronx,60+"
c5summarize$county_name[33]<-"NY,Kings,90+"
c5summarize$county_name[36]<-"NY,New York,50+"
c5summarize$county_name[41]<-"NY,Queens,80+"
#c5summarize$county_name[60]<-"MS,Rankin,320+"
c5summarize$county_name[47]<-"PA,Allegheny,80+"
#c5summarize$county_name[82]<-"NC,Wilson,270+"
c5summarize$county_name[54]<-"PA,Lancaster,50+"
c5summarize$county_name[60]<-"PA,Philadelphia,110+"


p3<-plot_usmap(fill="gray", alpha=0.3,include = c("PA","ME","NY","NJ","MA","VT","RI","CT","NH"),labels = TRUE) +
  ggrepel::geom_label_repel(data = c5summarize,
                            aes(x = long.1, y = lat.1, label = county_name),
                            size = 3, alpha = 0.8,
                            label.r = unit(0.5, "lines"), label.size = 0.5,
                            segment.color = "red", segment.size = 1,
                            seed = 1002) +
  geom_point(data = c5summarize, aes(x = long.1, y = lat.1, size = incident_cases),
             color = "violetred", alpha = 0.25) +
  scale_size_continuous(range=c(1,5))+
  labs(title = "Reported, June 29th",
       size = "Magnitude") +
  #theme(legend.position = "right")+
  theme(
    legend.position = "right",
    #text = element_text(color = "#22211d"),
    #plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    #panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    #legend.background = element_rect(fill = "#f5f5f2", color = NA),
    #plot.title = element_text(size= 16, hjust=0.1, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
  )


##################################
##################################   Sixth pic

b5=reg.dat.est.Jun29[,c("FIPS","incident_cases","EstR")]
c5=left_join(a,b5,by="FIPS")
c5[is.na(c5$EstR),]$EstR=0

names(c5)[1]<-"fips"

c6=c5[c5$state_abbv %in% c("NC","TN","MS","SC","AL","GA","FL"),]
c6summarize<- c6%>% group_by(fips)%>%
  summarize(state= unique(state_abbv), 
            long = mean(long,na.rm = TRUE),lat= mean(lat,na.rm = TRUE),
            incident_cases  =  unique(incident_cases))
c6summarize=c6summarize[!is.na(c6summarize$incident_cases),]
c6summarize<-c6summarize[c(3,4,1,2,5)]

c6summarize <- usmap_transform(c6summarize)
c6summarize$county_name<-NA
most_incident_cases_city6=rev(order(c6summarize$incident_cases))[1:9]
city.for.label6=c6summarize[most_incident_cases_city6,]

############  Match city/county name

city.name.for.label6=my.county.fips[my.county.fips$fips %in% city.for.label6$fips,]
for (i in 1:nrow(city.name.for.label6)) {
  c6summarize[c6summarize$fips==city.name.for.label6$fips[i],]$county_name=
    city.name.for.label6$polyname[i]
}

############  Manual label
c6summarize$county_name[14]<-"FL,Broward,600+"
c6summarize$county_name[16]<-"FL,Duval,500+"
c6summarize$county_name[19]<-"FL,Hillsborough,700+"
c6summarize$county_name[20]<-"FL,Lee,300+"
c6summarize$county_name[24]<-"FL,Miami-dade,1600+"
c6summarize$county_name[25]<-"FL,Orange,800+"
#c5summarize$county_name[60]<-"MS,Rankin,320+"
c6summarize$county_name[27]<-"FL,Palm Beach,400+"
#c5summarize$county_name[82]<-"NC,Wilson,270+"
c6summarize$county_name[29]<-"FL,Pinellas,300+"
c6summarize$county_name[73]<-"NC,Mecklenburg,300+"


p5<-plot_usmap(fill="gray", alpha=0.3,include = c("NC","TN","MS","SC","AL","GA","FL"),labels = TRUE) +
  ggrepel::geom_label_repel(data = c6summarize,
                            aes(x = long.1, y = lat.1, label = county_name),
                            size = 7.5, alpha = 0.8,
                            label.r = unit(0.5, "lines"), label.size = 0.5,
                            segment.color = "red", segment.size = 1,
                            seed = 1002) +
  geom_point(data = c6summarize, aes(x = long.1, y = lat.1, size = incident_cases),
             color = "violetred", alpha = 0.25) +
  scale_size_continuous(range=c(1,54))+
  labs(title = "Reported Cases, June 29th",
       size = "Magnitude") +
  theme(plot.title = element_text(size=19))+
  #theme(legend.position = "right")+
  theme(
    legend.position = "none",
    #text = element_text(color = "#22211d"),
    #plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    #panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    #legend.background = element_rect(fill = "#f5f5f2", color = NA),
    #plot.title = element_text(size= 16, hjust=0.1, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
  )


##################################
##################################   Eighth pic

b9=reg.dat.est.Jul06[,c("FIPS","incident_cases","EstR")]
c9=left_join(a,b9,by="FIPS")
c9[is.na(c9$EstR),]$EstR=0

names(c9)[1]<-"fips"

#c9=c5[c5$state_abbv %in% c("NC","TN","MS","SC","AL","GA","FL"),]
c9summarize<- c9%>% group_by(fips)%>%
  summarize(state= unique(state_abbv), 
            long = mean(long,na.rm = TRUE),lat= mean(lat,na.rm = TRUE),
            incident_cases  =  unique(incident_cases))
c9summarize=c9summarize[!is.na(c9summarize$incident_cases),]
c9summarize<-c9summarize[c(3,4,1,2,5)]

c9summarize <- usmap_transform(c9summarize)
c9summarize$county_name<-NA
most_incident_cases_city9=rev(order(c9summarize$incident_cases))[1:9]
city.for.label9=c9summarize[most_incident_cases_city9,]

############  Match city/county name

city.name.for.label9=my.county.fips[my.county.fips$fips %in% city.for.label9$fips,]
for (i in 1:nrow(city.name.for.label9)) {
  c9summarize[c9summarize$fips==city.name.for.label9$fips[i],]$county_name=
    city.name.for.label9$polyname[i]
}

############  Manual label
c9summarize$county_name[12]<-"AZ,Maricopa,2500+"
c9summarize$county_name[29]<-"CA,Los Angeles,2200+"
c9summarize$county_name[31]<-"CA,Orange,700+"
c9summarize$county_name[60]<-"FL,Broward,1100+"
c9summarize$county_name[70]<-"FL,Miami-dade,2100+"
c9summarize$county_name[71]<-"FL,Orange,700+"
c9summarize$county_name[193]<-"NV,Clark,600+"
c9summarize$county_name[304]<-"TX,Dallas,900+"
c9summarize$county_name[309]<-"TX,Harris,900+"


p6<-plot_usmap(fill="gray", alpha=0.3,labels = TRUE) +
  ggrepel::geom_label_repel(data = c9summarize,
                            aes(x = long.1, y = lat.1, label = county_name),
                            size = 3, alpha = 0.8,
                            label.r = unit(0.5, "lines"), label.size = 0.5,
                            segment.color = "red", segment.size = 1,
                            seed = 1002) +
  geom_point(data = c9summarize, aes(x = long.1, y = lat.1, size = incident_cases),
             color = "violetred", alpha = 0.25) +
  scale_size_continuous(range=c(1,17.5))+
  labs(title = "Reported, July 06th",
       size = "Magnitude") +
  #theme(legend.position = "right")+
  theme(
    legend.position = "right",
    #text = element_text(color = "#22211d"),
    #plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    #panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    #legend.background = element_rect(fill = "#f5f5f2", color = NA),
    #plot.title = element_text(size= 16, hjust=0.1, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
  )


###############################################################
###############################  Ninth Pic

b11=reg.dat.est.Apr07[,c("FIPS","cases","EstR")]
c11=left_join(a,b11,by="FIPS")
c11[is.na(c11$EstR),]$EstR=0
names(c11)[1]<-"fips"


#c1=c1[c1$state_abbv %in% c("PA","ME","NY","NJ","MA","VT","RI","CT","NH"),]
c11summarize<- c11%>% group_by(fips)%>%
  summarize(state= unique(state_abbv), 
            long = mean(long,na.rm = TRUE),lat= mean(lat,na.rm = TRUE),
            cases  =  unique(cases))
c11summarize=c11summarize[!is.na(c11summarize$cases),]
c11summarize<-c11summarize[c(3,4,1,2,5)]


# cities_t <- usmap_transform(citypop)
# cities_t1 <-cities_t[cities_t$abbr %in% c("PA","ME","NY","NJ","MA","VT","RI","CT","NH"),]

c11summarize <- usmap_transform(c11summarize)
c11summarize$county_name<-NA
most_cases_city11=rev(order(c11summarize$cases))[1:9]
city.for.label11=c11summarize[most_cases_city11,]

############  Match city/county name

city.name.for.label11=my.county.fips[my.county.fips$fips %in% city.for.label11$fips,]
for (i in 1:nrow(city.name.for.label11)) {
  c11summarize[c11summarize$fips==city.name.for.label11$fips[i],]$county_name=
    city.name.for.label11$polyname[i]
}

############  Manual label
c11summarize$county_name[79]<-"IL,Cook,9500+"
c11summarize$county_name[150]<-"MI,Wayne,9000+"
c11summarize$county_name[185]<-"NY,Bronx,15300+"
c11summarize$county_name[188]<-"NY,Kings,20800+"
c11summarize$county_name[190]<-"NY,Nassau,16600+"
c11summarize$county_name[191]<-"NY,New York,11000+"
c11summarize$county_name[197]<-"NY,Queens,24800+"
c11summarize$county_name[200]<-"NY,Suffolk,14500+"
c11summarize$county_name[203]<-"NY,WestChester,14800+"

####################################################


p11<-plot_usmap(fill="gray", alpha=0.3,labels = TRUE) +
  ggrepel::geom_label_repel(data = c11summarize,
                            aes(x = long.1, y = lat.1, label = county_name),
                            size = 3, alpha = 0.8,
                            label.r = unit(0.5, "lines"), label.size = 0.5,
                            segment.color = "red", segment.size = 1,
                            seed = 1002) +
  geom_point(data = c11summarize, aes(x = long.1, y = lat.1, size = cases),
             color = "violetred", alpha = 0.25) +
  scale_size_continuous(range=c(1,14))+
  labs(title = "Cumulated Cases, April 07th",
       #subtitle = "Reported, April 07th",
       size = "Magnitude") +
  #theme(legend.position = "right")+
  theme(
    legend.position = "top",
    #text = element_text(color = "#22211d"),
    #plot.background = element_rect(fill = "#f5f5f2", color = NA), 
    #panel.background = element_rect(fill = "#f5f5f2", color = NA), 
    #legend.background = element_rect(fill = "#f5f5f2", color = NA),
    #plot.title = element_text(size= 16, hjust=0.1, color = "#4e4d47", margin = margin(b = -0.1, t = 0.4, l = 2, unit = "cm")),
  )

##################################





################################################################
###########Combine four pics

ggdraw() +
  draw_plot(p1, 0, .5, .33, .5) +
  draw_plot(p2, .33, .5, .33, .5) +
  draw_plot(p3, .66, 0.5, .33, .5) +
  draw_plot(p4, 0, 0, .33, .5)+
  draw_plot(p5, .33, 0, .33, .5)+
  draw_plot(p6, .66, 0, .33, .5)


ggdraw() +
  draw_plot(p2, 0, 0, .3, 1) +
  draw_plot(p4, .32, 0, .33, .9)+
  draw_plot(p5, .66, 0, .34, .93)


p2$layers[[2]]$aes_params$size <- 7.5
p4$layers[[2]]$aes_params$size <- 7.5
p5$layers[[2]]$aes_params$size <- 7.5
ggdraw() +
  draw_plot(p2, 0, .64, 1, .36) +
  draw_plot(p4, 0, .32, 1, .32)+
  draw_plot(p5, 0, 0, 1, .32)





