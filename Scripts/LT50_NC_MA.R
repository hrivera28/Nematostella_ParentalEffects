### Analysis of adaptive vs parental effects (NC v. MA)
### Script associated with manuscript Rivera et al.
### Plasticity in parental effects confers rapid thermal tolerance in N. vectensis larvae
### Copyright H.E. Rivera 

library(drc)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggpubr)
library(rcompanion)
library(multcomp)
load("MA-NC_data.RData")

################################## 
# Model different LT50s for each Treatment

for (fam in levels(All_data$Pop)){
  assign(paste(fam,".LTDiff", sep=""),drm(Alive/(Alive+Dead) ~ Temp,Treatment, weights = (Alive+Dead), data = subset(All_data, Pop==fam), fct = LL.2(), type = "binomial"))
  assign(paste(fam,".LTSame", sep=""),drm(Alive/(Alive+Dead) ~ Temp,Treatment, weights = (Alive+Dead), data = subset(All_data, Pop==fam), fct = LL.2(), type = "binomial", pmodels = list(~Treatment-1, ~1)))
}

####### Looking first at the NC Differences between control and short term heat stress 
anova(NC.LTDiff,NC.LTSame)
## p = 0.0011 so there is a difference between the two, there is an increase ~ 0.2 after NC goes through the STHS 
# MA Differences between control and short term heat stress 
anova(MA.LTDiff, MA.LTSame)
# p = 0.0003, increase of 0.27 

####### Testing for difference between NC and MA under control and STHS treatments
NCMAC.LTDiff <- drm(Alive/(Alive+Dead) ~ Temp, Pop, weights = (Alive+Dead), data = subset(subset(All_data, Pop == "NC" | Pop=="MA"),Treatment =="Ctrl"), fct = LL.2(), type = "binomial")
summary(NCMAC.LTDiff)
## One in which we model the LT50s as being the same for both Treatments
NCMAC.LTSame <- drm(Alive/(Alive+Dead) ~ Temp, Pop, weights = (Alive+Dead), data = subset(subset(All_data, Pop == "NC" | Pop=="MA"),Treatment =="Ctrl"), fct = LL.2(), type = "binomial",pmodels = list(~factor(Pop)-1, ~1))
summary(NCMAC.LTSame)
## Then we can test the difference between the models using a Chi-Squared test 
anova(NCMAC.LTDiff,NCMAC.LTSame)
# p = 0.0088 so there's a difference, NC is about 0.2 higher 


## After the short term heat stress 
NCMAS.LTDiff <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = subset(subset(All_data, Pop == "NC" | Pop=="MA"),Treatment =="STHS"), fct = LL.2(), type = "binomial")
summary(NCMAS.LTDiff)
## One in which we model the LT50s as being the same for both Treatments
NCMAS.LTSame <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = subset(subset(All_data, Pop == "NC" | Pop=="MA"),Treatment =="STHS"), fct = LL.2(), type = "binomial", pmodels = list(~factor(Pop)-1, ~1))
summary(NCMAS.LTSame)
## Then we can test the difference between the models using a Chi-Squared test 
anova(NCMAS.LTDiff,NCMAS.LTSame)
# p = 0.0911 so the difference is no longer significant. 


### Between MA STHS and NC Control 
NCCMAS.LTDiff <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = rbind(subset(All_data, Pop == "NC" & Treatment=="Ctrl"),subset(All_data,  Pop=="MA" & Treatment =="STHS")), fct = LL.2(), type = "binomial")
summary(NCMAS.LTDiff)
## One in which we model the LT50s as being the same for both Treatments
NCCMAS.LTSame <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = rbind(subset(All_data, Pop == "NC" & Treatment=="Ctrl"),subset(All_data,  Pop=="MA" & Treatment =="STHS")), fct = LL.2(), type = "binomial", pmodels = list(~factor(Pop)-1, ~1))
summary(NCMAS.LTSame)
## Then we can test the difference between the models using a Chi-Squared test 
anova(NCCMAS.LTDiff,NCCMAS.LTSame)
# p = 0.17 there is not a significant difference in larval LT50's now. 

####### Between MA and the hybrids
MAHy.LTDiff <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = subset(subset(All_data, Pop == "MA" | Pop=="NC-Fem" | Pop=="NC-Male"),Treatment =="Ctrl"), fct = LL.2(), type = "binomial")
summary(MAHy.LTDiff)
## One in which we model the LT50s as being the same for both Treatments
MAHy.LTSame <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = subset(subset(All_data, Pop == "MA" | Pop=="NC-Fem" | Pop=="NC-Male"),Treatment =="Ctrl"), fct = LL.2(), type = "binomial",  pmodels = list(~factor(Pop)-1, ~1))
summary(MAHy.LTSame)
## Then we can test the difference between the models using a Chi-Squared test 
anova(MAHy.LTDiff,MAHy.LTSame)
# p = 0.41 so it's not really that different. Seems like the MA genes are really dragging these guys down
      
# After STHS
MAHy2.LTDiff <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = subset(subset(All_data, Pop == "MA" | Pop=="NC-Male"| Pop=="NC-Fem"),Treatment =="STHS"), fct = LL.2(), type = "binomial")
summary(MAHy2.LTDiff)
## One in which we model the LT50s as being the same for both Treatments
MAHy2.LTSame <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = subset(subset(All_data, Pop == "MA" | Pop=="NC-Male" | Pop=="NC-Fem"),Treatment =="STHS"), fct = LL.2(), type = "binomial", pmodels = list(~factor(Pop)-1, ~1))
summary(MAHy2.LTSame)
## Then we can test the difference between the models using a Chi-Squared test 
anova(MAHy2.LTDiff,MAHy2.LTSame)
# p = 0.6 also the same then too.

# Between the two different hybrid sets
# Under control conditions
HyC.LTDiff <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = subset(subset(All_data, Pop == "NC-Fem" | Pop=="NC-Male"),Treatment =="Ctrl"), fct = LL.2(), type = "binomial")
summary(HyC.LTDiff)
## One in which we model the LT50s as being the same for both Treatments
HyC.LTSame <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = subset(subset(All_data, Pop == "NC-Fem" | Pop=="NC-Male"),Treatment =="Ctrl"), fct = LL.2(), type = "binomial", pmodels = list(~factor(Pop)-1, ~1))
summary(HyC.LTSame)
## Then we can test the difference between the models using a Chi-Squared test 
anova(HyC.LTDiff,HyC.LTSame)
# p = 0.7

# After STHS 
HyS.LTDiff <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = subset(subset(All_data, Pop == "NC-Fem" | Pop=="NC-Male"),Treatment =="STHS"), fct = LL.2(), type = "binomial")
summary(HyS.LTDiff)
## One in which we model the LT50s as being the same for both Treatments
HyS.LTSame <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = subset(subset(All_data, Pop == "NC-Fem" | Pop=="NC-Male"),Treatment =="STHS"), fct = LL.2(), type = "binomial", pmodels = list(~factor(Pop)-1, ~1))
summary(HyS.LTSame)
## Then we can test the difference between the models using a Chi-Squared test 
anova(HyS.LTDiff,HyS.LTSame)
# p = 0.65 also the same then too.


# NC female Hybrids pre and post hs
HySF.LTDiff <- drm(Alive/(Alive+Dead) ~ Temp,Treatment, weights = (Alive+Dead), data = subset(All_data, Pop == "NC-Fem"), fct = LL.2(), type = "binomial")
summary(HySF.LTDiff)
## One in which we model the LT50s as being the same for both Treatments
HySF.LTSame <- drm(Alive/(Alive+Dead) ~ Temp,Treatment, weights = (Alive+Dead), data = subset(All_data, Pop == "NC-Fem"), fct = LL.2(), type = "binomial", pmodels = list(~factor(Treatment)-1, ~1))
summary(HySF.LTSame)
## Then we can test the difference between the models using a Chi-Squared test 
anova(HySF.LTDiff,HySF.LTSame)
# p = 0.001 so there's a significant increase 

#NC Dad 
HySM.LTDiff <- drm(Alive/(Alive+Dead) ~ Temp,Treatment, weights = (Alive+Dead), data = subset(All_data, Pop == "NC-Male"), fct = LL.2(), type = "binomial")
summary(HySM.LTDiff)
## One in which we model the LT50s as being the same for both Treatments
HySM.LTSame <- drm(Alive/(Alive+Dead) ~ Temp,Treatment, weights = (Alive+Dead), data = subset(All_data, Pop == "NC-Male"), fct = LL.2(), type = "binomial", pmodels = list(~factor(Treatment)-1, ~1))
summary(HySM.LTSame)
## Then we can test the difference between the models using a Chi-Squared test 
anova(HySM.LTDiff,HySM.LTSame)
# p = 0.09 so there's not a significant increase 


#### Doing all the pops together by treatment for error estimations
# All Control
Ctrl_NM.drm <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = subset(All_data, Treatment =="Ctrl"), fct = LL.2(), type = "binomial")
summary(Ctrl_NM.drm)

# All STHS 
STHS_NM.drm <- drm(Alive/(Alive+Dead) ~ Temp,Pop, weights = (Alive+Dead), data = subset(All_data,Treatment =="STHS"), fct = LL.2(), type = "binomial")
summary(STHS_NM.drm)

# Getting estimates corrected for simultaneous inference 
Ctrl_NM_SI<-summary(glht(Ctrl_NM.drm))
STHS_NM_SI<-summary(glht(STHS_NM.drm))
# The ED estimates live in Ctrl_NM_SI$test$coefficients[5:8]
# The SE measures live in Ctrl_NM_SI$test$sigma[5:8]


## Use the robust estimates and standard errors for the difference plots
LT_mat<-as.data.frame(matrix(ncol=4, nrow=8))
colnames(LT_mat)<-c("LT50", "Cohort", "Treatment", "SE")

for (i in seq(1,4)){
  LT_mat$LT50[i]<-Ctrl_NM_SI$test$coefficients[i+4] #LT50 estimate for Control
  LT_mat$SE[i]<-Ctrl_NM_SI$test$sigma[i+4] #SE estimate for Control
  LT_mat$Treatment[i]<-"Ctrl"
}

for (i in seq(5,8)){
  LT_mat$LT50[i]<-STHS_NM_SI$test$coefficients[i] #LT50 estimate for STHS
  LT_mat$SE[i]<-STHS_NM_SI$test$sigma[i] #SE estimate for STHS
  LT_mat$Treatment[i]<-"STHS"
}

LT_mat$Cohort<-c("NC-Mom", "NC", "NC-Dad", "MA", "NC", "NC-Mom","NC-Dad", "MA")
# For some reason the order of the families in the Ctrl model is different from that of the STHS model

# For thesis chapter
(D3<-ggplot(LT_mat, aes(x=Treatment, y=LT50, group=Cohort, colour=Cohort))+
    geom_point(position=position_dodge(width = 0.15))+
    geom_line(position=position_dodge(width = 0.15))+theme_minimal()+
    geom_errorbar(aes(ymin=LT50-SE, ymax=LT50+SE),position=position_dodge(width = 0.15), width=0)+
    xlab("")+ggtitle("")+ylab("LT50 (°C)")+
    scale_y_continuous(breaks=c(41.2,41.4,41.6), limits=c(41.1,41.7))+
    scale_x_discrete(labels=c("Controls", "STHS"),expand=c(0.1,0.1))+
    theme(legend.position = "bottom", 
          legend.background = element_rect(colour="grey"),
          legend.text = element_text(size=8, face="bold",colour="black"),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=10, colour="black", face="bold"), 
          axis.text=element_text(size=8, face="bold",colour="black"))+
    labs(colour="")+scale_colour_manual(values=c("cornflowerblue", "coral3", "darkolivegreen4","plum4"))+
    annotate("text", label="*", x=c(0.93,0.895,1.1), y=c(41.38,41.2,41.32), size=5)+
    annotate("text", label="+", x=2.01, y=41.485, size=4))

ggsave("Fig3D.png",  width=4, height=4, units="in", dpi=300)


#### Plots of response curve by treatment and pop
# For Control Trial 
MA_CT <- drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(All_data, Pop=="MA" & Treatment=="Ctrl"), fct = LL.2(), type = "binomial")
predict_MA_CT<-predict(MA_CT, newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence")
predict_MA_CT<-cbind(predict_MA_CT,expand.grid(Temp=seq(39.8, 42.3, length=25)))
colnames(predict_MA_CT)[4]<-"Temp"

NC_CT <- drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(All_data, Pop=="NC" & Treatment=="Ctrl"), fct = LL.2(), type = "binomial")
predict_NC_CT<-predict(NC_CT, newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence")
predict_NC_CT<-cbind(predict_NC_CT,expand.grid(Temp=seq(39.8, 42.3, length=25)))
colnames(predict_NC_CT)[4]<-"Temp"

NC_fem_CT <- drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(All_data, Pop=="NC-Fem" & Treatment=="Ctrl"), fct = LL.2(), type = "binomial")
predict_NC_femCT<-predict(NC_fem_CT, newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence")
predict_NC_femCT<-cbind(predict_NC_femCT,expand.grid(Temp=seq(39.8, 42.3, length=25)))
predict_NC_femCT$Upper[predict_NC_femCT$Upper>1.05]<-1.05
colnames(predict_NC_femCT)[4]<-"Temp"

NC_male_CT <- drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(All_data, Pop=="NC-Male" & Treatment=="Ctrl"), fct = LL.2(), type = "binomial")
predict_NC_maleCT<-predict(NC_male_CT, newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence")
predict_NC_maleCT<-cbind(predict_NC_maleCT,expand.grid(Temp=seq(39.8, 42.3, length=25)))
predict_NC_maleCT$Upper[predict_NC_maleCT$Upper>1.05]<-1.05
colnames(predict_NC_maleCT)[4]<-"Temp"

(B3<-ggplot(subset(All_data, Treatment=="Ctrl" & Pop=="NC"), aes(x=Temp, y=Alive/(Alive+Dead)))+
    geom_point(colour="coral4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
    geom_ribbon(data=predict_NC_CT,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral3")+
    geom_line(data=predict_NC_CT, aes(x=Temp, y=Prediction), colour="coral4")+
    geom_point(data=subset(All_data, Treatment=="Ctrl" & Pop=="MA"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="lightblue4")+
    geom_ribbon(data=predict_MA_CT,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightblue3")+
    geom_line(data=predict_MA_CT, aes(x=Temp, y=Prediction),colour="lightblue4")+
    geom_point(data=subset(All_data, Treatment=="Ctrl" & Pop=="NC-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="plum4")+
    geom_ribbon(data=predict_NC_femCT,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="plum3")+
    geom_line(data=predict_NC_femCT, aes(x=Temp, y=Prediction),  colour="plum4")+
    geom_point(data=subset(All_data, Treatment=="Ctrl" & Pop=="NC-Male"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="darkolivegreen4")+
    geom_ribbon(data=predict_NC_maleCT,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="darkolivegreen3")+
    geom_line(data=predict_NC_maleCT, aes(x=Temp, y=Prediction),colour="darkolivegreen4")+
    ylab("Survival")+xlab("Temperature (°C)")+
    scale_y_continuous(breaks=c(0,0.5,1), labels=c("0%", "50%", "100%"))+
    coord_cartesian(ylim=c(0,1.1))+
    scale_x_continuous(breaks=c(40,40.5,41,41.5,42))+
    theme(axis.text=element_text(size=8, face="bold", colour="black"),
          axis.title=element_text(size=10, face="bold", colour="black")))


# For STHS Trial 
MA_HS <- drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(All_data, Pop=="MA" & Treatment=="STHS"), fct = LL.2(), type = "binomial")
predict_MA_HS<-predict(MA_HS, newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence")
predict_MA_HS<-cbind(predict_MA_HS,expand.grid(Temp=seq(39.8, 42.3, length=25)))
colnames(predict_MA_HS)[4]<-"Temp"

NC_HS <- drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(All_data, Pop=="NC" & Treatment=="STHS"), fct = LL.2(), type = "binomial")
predict_NC_HS<-predict(NC_HS, newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence")
predict_NC_HS<-cbind(predict_NC_HS,expand.grid(Temp=seq(39.8, 42.3, length=25)))
colnames(predict_NC_HS)[4]<-"Temp"

NC_fem_HS <- drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(All_data, Pop=="NC-Fem" & Treatment=="STHS"), fct = LL.2(), type = "binomial")
predict_NC_femHS<-predict(NC_fem_HS, newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence")
predict_NC_femHS<-cbind(predict_NC_femHS,expand.grid(Temp=seq(39.8, 42.3, length=25)))
predict_NC_femHS$Upper[predict_NC_femHS$Upper>1.05]<-1.05
colnames(predict_NC_femHS)[4]<-"Temp"

NC_male_HS <- drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(All_data, Pop=="NC-Male" & Treatment=="STHS"), fct = LL.2(), type = "binomial")
predict_NC_maleHS<-predict(NC_male_HS, newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence")
predict_NC_maleHS<-cbind(predict_NC_maleHS,expand.grid(Temp=seq(39.8, 42.3, length=25)))
colnames(predict_NC_maleHS)[4]<-"Temp"

(C3<-ggplot(subset(All_data, Treatment=="STHS" & Pop=="NC"), aes(x=Temp, y=Alive/(Alive+Dead)))+
    geom_point(colour="coral4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
    geom_ribbon(data=predict_NC_HS,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral3")+
    geom_line(data=predict_NC_HS, aes(x=Temp, y=Prediction), colour="coral4")+
    geom_point(data=subset(All_data, Treatment=="STHS" & Pop=="MA"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="lightblue4")+
    geom_ribbon(data=predict_MA_HS,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightblue3")+
    geom_line(data=predict_MA_HS, aes(x=Temp, y=Prediction),colour="lightblue4")+
    geom_point(data=subset(All_data, Treatment=="STHS" & Pop=="NC-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="plum4")+
    geom_ribbon(data=predict_NC_femHS,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="plum3")+
    geom_line(data=predict_NC_femHS, aes(x=Temp, y=Prediction),  colour="plum4")+
    geom_point(data=subset(All_data, Treatment=="STHS" & Pop=="NC-Male"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="darkolivegreen4")+
    geom_ribbon(data=predict_NC_maleHS,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="darkolivegreen3")+
    geom_line(data=predict_NC_maleHS, aes(x=Temp, y=Prediction),colour="darkolivegreen4")+
    ylab("")+xlab("Temperature (°C)")+
    scale_y_continuous(breaks=c(0,0.5,1))+
    coord_cartesian(ylim=c(0,1.1))+
    scale_x_continuous(breaks=c(40,40.5,41,41.5,42))+
    theme(axis.text=element_text(size=8, face="bold", colour="black"),
          axis.text.y=element_blank(),
          axis.title=element_text(size=10, face="bold", colour="black")))

ggarrange(B3, C3, D3, nrow=1, ncol=3, common.legend = T, align = "h",legend="bottom")
ggsave("Fig3BD.png", width=7, height=3, units="in", dpi=300)
