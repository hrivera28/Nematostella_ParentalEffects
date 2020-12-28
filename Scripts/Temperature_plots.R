# Script used to visualize and analyze
# In situ logger data from Sachkova et al. 2020 
# https://bmcbiol.biomedcentral.com/articles/10.1186/s12915-020-00855-8#MOESM1
# Data S6 

##################
# Load Libraries 
library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(xts)
library(zoo)
library(TTR)
library(scales)
library(ggpubr)
library(signal)
library(lubridate)

###### 
# Load data
load("MA_NC_temps.Rdata")

#Weekly average data
MA_wk<-apply.monthly(MA, FUN=mean)
NC_wk<-apply.monthly(NC, FUN=mean)

plot(MA_wk)
plot(NC_wk)

######
## Create data frame with all MA-NC temperatures and plot 
# Name columns
MA_gg<-as.data.frame(MA)
row.names(MA_gg)<-as.character(seq(1,length(MA_gg[,1])))
# add in date time values as variable
MA_gg<-cbind(time(MA),MA_gg)
# fix column header 
colnames(MA_gg)[1]<-c("DateTime")
# convert to long format for use in ggplot 
# The datetime values are repeated over each site/temp combination 
colnames(MA_gg)<-c("DateTime", "MA_t")

NC_gg<-as.data.frame(NC)
row.names(NC_gg)<-as.character(seq(1,length(NC_gg[,1])))
# add in date time values as variable
NC_gg<-cbind(time(NC),NC_gg)
# fix column header 
colnames(NC_gg)[1]<-c("DateTime")
# convert to long format for use in ggplot 
# The datetime values are repeated over each site/temp combination 
colnames(NC_gg)<-c("DateTime", "NC_t")

all_temp<-left_join(NC_gg, MA_gg, by="DateTime")
pivot_longer(all_temp, cols=contains("_t"), names_to = "Site", values_to = "Temp")->all_temp

##### 
#Just summer months
ggplot(all_temp,aes(x=DateTime, y=Temp,colour=Site))+geom_line(size=.5, alpha=0.5)+theme_minimal()+
  xlab("")+ylab("Temperature (°C)")+theme(panel.grid.minor = element_blank(), 
                                          axis.text = element_text(face="bold", colour="black"), 
                                          axis.title = element_text(face="bold", colour="black"),
                                          legend.text = element_text(face="bold", colour="black"),
                                          legend.title =  element_text(face="bold", colour="black"),
                                          legend.position = "bottom", 
                                          legend.background = element_rect(fill=NA, colour="grey80"))+
  scale_y_continuous(limits=c(10,45), breaks=c(10,20,30,40), labels=c("10", "20", "30", "40"))+
  scale_x_datetime(limits=ymd_hms(c("2016-06-01 00:00:00", "2016-09-30 00:00:00")), 
                   breaks=ymd_hms(c("2016-06-01 00:00:00", "2016-07-01 00:00:00", 
                                    "2016-08-01 00:00:00", "2016-09-01 00:00:00", 
                                    "2016-09-30 00:00:00")), 
                   labels=c("June", "July", "August", "Sept", "Oct"))+
  scale_colour_manual(values=c("cornflowerblue", "coral3"), labels=c("MA", "NC"))

# Fig S1A 
ggsave("Summer_NC_MA.png", width = 6, height=4, dpi=300, units="in")



#### 
# Just one MA week in midsummer
ggplot(subset(all_temp,Site=="MA_t"), aes(x=DateTime, y=Temp,colour=Site))+geom_line(size=0.8)+theme_minimal()+
  xlab("Date")+ylab("Temperature (°C)")+theme(panel.grid.minor = element_blank(), 
                                          axis.text = element_text(face="bold", colour="black"), 
                                          axis.title = element_text(face="bold", colour="black"))+
  scale_y_continuous(limits=c(18,42), breaks=c(20,30,40), labels=c("20", "30", "40"))+
  scale_x_datetime(limits=ymd_hms(c("2016-07-15 00:00:00", "2016-07-22 00:00:00")), 
                   breaks=ymd_hms(c("2016-07-15 00:00:00", "2016-07-16 00:00:00", 
                                    "2016-07-17 00:00:00", "2016-07-18 00:00:00", 
                                    "2016-07-19 00:00:00", "2016-07-20 00:00:00", 
                                    "2016-07-21 00:00:00", "2016-07-22 00:00:00")), 
                   labels=c("7/15", "7/16", "7/17", "7/18", "7/19", "7/20", "7/21","7/22"))+
  scale_colour_manual(values="cornflowerblue")+guides(colour=F)+
  geom_hline(aes(yintercept = 33), colour="grey30", linetype="dashed")

# Fig 1A
ggsave("MA_summer_Week.png", width = 3, height=3.5, dpi=300, units="in")



######
# Daily mean summer temps
MA_mean<-apply.daily(MA, FUN=mean)
NC_mean<-apply.daily(NC, FUN=mean)

MA_mean_gg<-as.data.frame(MA_mean)
row.names(MA_mean_gg)<-as.character(seq(1,length(MA_mean_gg[,1])))
# add in date time values as variable
MA_mean_gg<-cbind(time(MA_mean),MA_mean_gg)
# fix column header 
colnames(MA_mean_gg)[1]<-c("DateTime")
# convert to long format for use in ggplot 
# The datetime values are repeated over each site/temp combination 
colnames(MA_mean_gg)<-c("DateTime", "MA_mean_t")

NC_mean_gg<-as.data.frame(NC_mean)
row.names(NC_mean_gg)<-as.character(seq(1,length(NC_mean_gg[,1])))
# add in date time values as variable
NC_mean_gg<-cbind(time(NC_mean),NC_mean_gg)
# fix column header 
colnames(NC_mean_gg)[1]<-c("DateTime")
# convert to long format for use in ggplot 
# The datetime values are repeated over each site/temp combination 
colnames(NC_mean_gg)<-c("DateTime", "NC_mean_t")

all_mean<-left_join(NC_mean_gg, MA_mean_gg, by="DateTime")
pivot_longer(all_mean, cols=contains("_t"), names_to = "Site", values_to = "Mean")->all_mean

ggplot(all_mean, aes(x=Site, y=Mean, colour=Site, fill=Site))+geom_boxplot(aes(alpha=0.5))+theme_minimal()+
  scale_x_discrete(labels=c("MA", "NC"))+
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(face="bold", colour="black"), 
        axis.title = element_text(face="bold", colour="black"))+
  ylab("Daily Mean Temperature (°C)")+annotate("text", x=1, y=32, label="*", size=7)+
  scale_colour_manual(values=c("cornflowerblue", "coral3"), labels=c("MA", "NC"))+
  scale_fill_manual(values=c("cornflowerblue", "coral3"), labels=c("MA", "NC"))+
  guides(colour=F, alpha=F, fill=F)

t.test(MA_mean$V1, NC_mean$V1)

#Fig S1B
ggsave("MA_NC_meanDTemp.png", width = 3, height=4, dpi=300, units="in")



#####
# Summer daily variability 
MA_dtv<-apply.daily(MA, FUN=max)-apply.daily(MA, FUN=min)  
NC_dtv<-apply.daily(NC, FUN=max)-apply.daily(NC, FUN=min)  

MA_dtv_summer<-apply.daily(MA["2016-06/2016-09"], FUN=max)-apply.daily(MA["2016-06/2016-09"], FUN=min)  
NC_dtv_summer<-apply.daily(NC["2016-06/2016-09"], FUN=max)-apply.daily(NC["2016-06/2016-09"], FUN=min)  


MA_dtv_gg<-as.data.frame(MA_dtv)
row.names(MA_dtv_gg)<-as.character(seq(1,length(MA_dtv_gg[,1])))
# add in date time values as variable
MA_dtv_gg<-cbind(time(MA_dtv),MA_dtv_gg)
# fix column header 
colnames(MA_dtv_gg)[1]<-c("DateTime")
# convert to long format for use in ggplot 
# The datetime values are repeated over each site/temp combination 
colnames(MA_dtv_gg)<-c("DateTime", "MA_dtv_t")

NC_dtv_gg<-as.data.frame(NC_dtv)
row.names(NC_dtv_gg)<-as.character(seq(1,length(NC_dtv_gg[,1])))
# add in date time values as variable
NC_dtv_gg<-cbind(time(NC_dtv),NC_dtv_gg)
# fix column header 
colnames(NC_dtv_gg)[1]<-c("DateTime")
# convert to long format for use in ggplot 
# The datetime values are repeated over each site/temp combination 
colnames(NC_dtv_gg)<-c("DateTime", "NC_dtv_t")

all_dtv<-left_join(NC_dtv_gg, MA_dtv_gg, by="DateTime")
pivot_longer(all_dtv, cols=contains("_t"), names_to = "Site", values_to = "DTV")->all_dtv

ggplot(all_dtv, aes(x=Site, y=DTV, colour=Site, fill=Site))+geom_boxplot(aes(alpha=0.5))+theme_minimal()+
  scale_x_discrete(labels=c("MA", "NC"))+
  theme(panel.grid.minor = element_blank(), 
        axis.text = element_text(face="bold", colour="black"), 
        axis.title = element_text(face="bold", colour="black"))+
  ylab("Daily Temperature Range (°C)")+
  scale_colour_manual(values=c("cornflowerblue", "coral3"), labels=c("MA", "NC"))+
  scale_fill_manual(values=c("cornflowerblue", "coral3"), labels=c("MA", "NC"))+
  guides(colour=F, alpha=F, fill=F)

#Fig S1C
ggsave("MA_NC_DTV.png", width = 3, height=4, dpi=300, units="in")

t.test(MA_dtv$V1, NC_dtv$V1) # No signficant difference in daily temperature range



