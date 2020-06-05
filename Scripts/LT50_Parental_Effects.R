### Analysis of parental effects 
### Script associated with manuscript Rivera et al.
### Plasticity in parental effects confers rapid thermal tolerance in N. vectensis larvae
### Copyright H.E. Rivera 

##### 
# Required dependencies:

library(drc)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggpubr)
library(rcompanion)
library(multcomp)
load("Main_parental.RData")

#####
# Run survival models for each cohort 

 for (fam in levels(all_data$Cohort)){
  # Model different LT 50s for each Treatment
  assign(paste(fam,".LTDiff", sep=""), drm(Alive/(Alive+Dead) ~ Temp,Treatment, weights = (Alive+Dead), data = subset(all_data, Cohort==fam), fct = LL.2(), type = "binomial"))
  # Model the LT50s as being the same for both Treatments
  assign(paste(fam,".LTSame", sep=""), drm(Alive/(Alive+Dead) ~ Temp,Treatment, weights = (Alive+Dead), data = subset(all_data, Cohort==fam), fct = LL.2(), type = "binomial", pmodels = list(~factor(Treatment)-1, ~1)))
}

# Look at results, then test the difference between the models
summary(P1.LTDiff)
summary(P1.LTSame)
anova(P1.LTDiff,P1.LTSame)
## For P1 we conclude the LT50 values are NOT the same between treatments (p=0.0021)

summary(P2.LTDiff)
summary(P2.LTSame)
anova(P2.LTDiff,P2.LTSame)
## For P2 we conclude the LT50 values are NOT the same between treatments (p=0.0005)

summary(P4.LTDiff)
summary(P4.LTSame)
anova(P4.LTDiff,P4.LTSame)
## For P4 we conclude the LT50 values are NOT the same between treatments (p=0.0174)

summary(G1.LTDiff)
summary(G1.LTSame)
anova(G1.LTDiff,G1.LTSame)
## For G1 we conclude the LT50 values are NOT the same between treatments (p=0.0013)

summary(G2.LTDiff)
summary(G2.LTSame)
anova(G2.LTDiff,G2.LTSame)
## For G2 we conclude the LT50 values are NOT the same between treatments (p=0.0)

summary(G3.LTDiff)
summary(G3.LTSame)
anova(G3.LTDiff,G3.LTSame)
##For G3 we conclude the LT50 values are NOT the same between treatments (p=0.0173)

summary(G4.LTDiff)
summary(G4.LTSame)
anova(G4.LTDiff,G4.LTSame)
## For G4 we conclude the LT50 values are ARE the same between treatments (p=0.003)

summary(G5.LTDiff)
summary(G5.LTSame)
anova(G5.LTDiff,G5.LTSame)
## For G5 we conclude the LT50 values are NOT the same between treatments (p=0.001)

summary(G6.LTDiff)
summary(G6.LTSame)
anova(G6.LTDiff,G6.LTSame)
## For G6 we conclude the LT50 values are NOT the same between treatments (p=0.0)
  
      
#######################################
# Test for differences between the groups using paired t-tests 
# Exctract the ED50 (LT50) estimate from each dose response model 
# and run a linear model by cohort and treatment to look at overall 
# trends on LT50 by family and parental treatment

# First let's make sure the differences follow a normal distribution 
#Getting the robust estimates
Ctrl_drm<-drm(Alive/(Alive+Dead) ~ Temp, Cohort, weights = (Alive+Dead), data = subset(all_data, Treatment=="Ctrl"), fct = LL.2(), type = "binomial")
STHS_drm<-drm(Alive/(Alive+Dead) ~ Temp, Cohort,weights = (Alive+Dead), data = subset(all_data, Treatment=="STHS"),  fct = LL.2(), type = "binomial")

#Get simulatneously inferred paremeters, standard errors and confidence intervals
Ctrl_SI<-summary(glht(Ctrl_drm))
STHS_SI<-summary(glht(STHS_drm))
# The ED estimates live in Ctrl_SI$test$coefficients[10:18]
# The SE measures live in Ctrl_SI$test$sigma[10:18]

# Make a dataframe with the LT50 Estimates and SE by cohort and treatment
LT_mat<-as.data.frame(matrix(ncol=4, nrow=18))
colnames(LT_mat)<-c("LT50", "Cohort", "Treatment", "SE")

for (i in seq(1,9)){
  LT_mat$LT50[i]<-Ctrl_SI$test$coefficients[i+9] #LT50 estimate for Control
  LT_mat$SE[i]<-Ctrl_SI$test$sigma[i+9] #SE estimate for Control
  LT_mat$Treatment[i]<-"Ctrl"
}

for (i in seq(10,18)){
  LT_mat$LT50[i]<-STHS_SI$test$coefficients[i] #LT50 estimate for STHS
  LT_mat$SE[i]<-STHS_SI$test$sigma[i] #SE estimate for STHS
  LT_mat$Treatment[i]<-"STHS"
}
LT_mat$Cohort<-rep(c("G1", "G2", "G3","G4", "G5", "G6", "P1", "P2", "P4"),2)

# Paired t-test of the data 
t.test(LT50~Treatment, data=LT_mat, paired=TRUE, alternative="less")
# We use less because the Ctrl treatment comes first so we want to know the probability that the LT50 is 
# lower in the Ctrls that the STHS 
# p = 2.25 e -7 so very significant, woo!!! 
# estimated mean difference is -0.34 C 


# Let's plot the differences by cohort all on the same plot
cbbPalette <- c("#000000", "#E69F00", "#5956A0", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")
# Panel E Figure 1
(E1<-ggplot(LT_mat, aes(x=Treatment, y=LT50, group=Cohort, colour=Cohort))+
    geom_point(position=position_dodge(width = 0.1))+
    geom_line(position=position_dodge(width = 0.1))+
    geom_errorbar(aes(ymin=LT50-SE, ymax=LT50+SE),width=0, position=position_dodge(width = 0.1))+
    theme_minimal()+theme(panel.grid.minor = element_blank(), legend.position = "bottom",
                          legend.background = element_rect(colour="gray"),
                          legend.text=element_text(size=10, colour="black", face="bold"),
                          legend.title=element_text(size=10, colour="black", face="bold"),
                          axis.text = element_text(size=8, colour="black", face="bold"), 
                          axis.title.y=element_text(size=10, colour="black", face="bold"))+
    xlab("")+ylab("LT50 (°C)")+
    scale_y_continuous(breaks=c(40.8,41,41.2,41.4,41.6,41.8), limits=c(40.75,41.8))+
    scale_color_manual(values = cbbPalette)+scale_x_discrete(labels=c("Controls", "STHS"), expand=c(0.1,0.1))+
    annotate("text", label="*", size=6, x=2, y=41.7)+labs(colour="Cohort")+guides(colour = guide_legend(nrow = 1)))

ggarrange(E1, F1, common.legend = TRUE, align="h", legend="bottom")    
ggsave("Fig1EF.png", width = 6, height = 4, units= "in", dpi=300)

#######################################
# Create individual plots for supplement

## For these graphs make sure the prediction is only on one of the treatments/trials 
## and then plot only those point results else you get weird looking graphs

for (fam in levels(all_data$Cohort)){
  # Same model as previous but use results to calucate interval for ribbon plot
  name1<-paste(fam,".Ctrl", sep="")
  assign(paste(fam,".Ctrl", sep=""), drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(all_data, Cohort==fam & Treatment=="Ctrl"), fct = LL.2(), type = "binomial"))
  assign(paste("predict_", fam,"Ctrl", sep=""), predict(get(name1), newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence"))
  name2<-paste("predict_", fam,"Ctrl", sep="")
  assign(paste("predict_", fam,"Ctrl", sep=""), cbind(get(name2),expand.grid(Temp=seq(39.8, 42.3, length=25))))
  
  # For the STHS treatment
  name3<-paste(fam,".STHS", sep="")
  assign(paste(fam,".STHS", sep=""),drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(all_data, Cohort==fam & Treatment=="STHS"),  fct = LL.2(), type = "binomial"))
  assign(paste("predict_", fam,"STHS", sep=""), predict(get(name3), newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence"))
  name4<-paste("predict_", fam,"STHS", sep="")
  assign(paste("predict_", fam,"STHS", sep=""),cbind(get(name4),expand.grid(Temp=seq(39.8, 42.3, length=25))))
}
  

A<-ggplot(subset(all_data, Cohort=="P1" & Treatment=="Ctrl"), aes(x=Temp, y=Alive/(Alive+Dead)))+
  geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
  geom_ribbon(data=predict_P1Ctrl,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightskyblue2")+
  geom_line(data=predict_P1Ctrl, aes(x=Temp, y=Prediction))+ggtitle("P1")+
  geom_point(data=subset(all_data, Cohort=="P1" & Treatment=="STHS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral3")+
  geom_ribbon(data=predict_P1STHS,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral2")+
  geom_line(data=predict_P1STHS, aes(x=Temp, y=Prediction))+ylab("")+xlab("")+
  scale_x_continuous(breaks=c(40,40.5,41,41.5,42))+
  scale_y_continuous(breaks=c(0,0.5,1), limits=c(-0.1,1.1), labels=c("0%", "50%", "100%"))+
  theme(axis.text.x =element_blank())


B<-ggplot(subset(all_data, Cohort=="P2" & Treatment=="Ctrl"), aes(x=Temp, y=Alive/(Alive+Dead)))+
  geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
  geom_ribbon(data=predict_P2Ctrl,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightskyblue2")+
  geom_line(data=predict_P2Ctrl, aes(x=Temp, y=Prediction))+ggtitle("P2")+
  geom_point(data=subset(all_data, Cohort=="P2" & Treatment=="STHS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral3")+
  geom_ribbon(data=predict_P2STHS,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral2")+
  geom_line(data=predict_P2STHS, aes(x=Temp, y=Prediction))+ylab("")+xlab("")+
  scale_x_continuous(breaks=c(40,40.5,41,41.5,42))+scale_y_continuous(breaks=c(0,0.5,1), limits=c(-0.1,1.1))+
  theme(axis.text=element_blank())


C<-ggplot(subset(all_data, Cohort=="P4" & Treatment=="Ctrl"), aes(x=Temp, y=Alive/(Alive+Dead)))+
  geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
  geom_ribbon(data=predict_P4Ctrl,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightskyblue2")+
  geom_line(data=predict_P4Ctrl, aes(x=Temp, y=Prediction))+ggtitle("P4")+
  geom_point(data=subset(all_data, Cohort=="P4" & Treatment=="STHS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral3")+
  geom_ribbon(data=predict_P4STHS,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral2")+
  geom_line(data=predict_P4STHS, aes(x=Temp, y=Prediction))+ylab("")+xlab("")+
  scale_x_continuous(breaks=c(40,40.5,41,41.5,42), labels=c("40°C","40.5°C","41°C","41.5°C","42°C"))+
  theme(axis.text =element_blank())+scale_y_continuous(breaks=c(0,0.5,1), limits=c(-0.1,1.1))

D<-ggplot(subset(all_data, Cohort=="G1" & Treatment=="Ctrl"), aes(x=Temp, y=Alive/(Alive+Dead)))+
  geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
  geom_ribbon(data=predict_G1Ctrl,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightskyblue2")+
  geom_line(data=predict_G1Ctrl, aes(x=Temp, y=Prediction))+ggtitle("G1")+
  geom_point(data=subset(all_data, Cohort=="G1" & Treatment=="STHS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral3")+
  geom_ribbon(data=predict_G1STHS,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral2")+
  geom_line(data=predict_G1STHS, aes(x=Temp, y=Prediction))+ylab("Survival")+xlab("")+
  scale_x_continuous(breaks=c(40,40.5,41,41.5,42), labels=c("40°C","40.5°C","41°C","41.5°C","42°C"))+
  scale_y_continuous(breaks=c(0,0.5,1), limits=c(-0.05,1.1), labels=c("0%", "50%", "100%"))+theme(axis.text.x =element_blank())


E<-ggplot(subset(all_data, Cohort=="G2" & Treatment=="Ctrl"), aes(x=Temp, y=Alive/(Alive+Dead)))+
  geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
  geom_ribbon(data=predict_G2Ctrl,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightskyblue2")+
  geom_line(data=predict_G2Ctrl, aes(x=Temp, y=Prediction))+ggtitle("G2")+
  geom_point(data=subset(all_data, Cohort=="G2" & Treatment=="STHS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral3")+
  geom_ribbon(data=predict_G2STHS,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral2")+
  geom_line(data=predict_G2STHS, aes(x=Temp, y=Prediction))+ylab("")+xlab("")+
  scale_x_continuous(breaks=c(40,40.5,41,41.5,42),labels=c("40°C","40.5°C","41°C","41.5°C","42°C"))+
  theme(axis.text.y= element_blank(),axis.text.x =element_blank())+scale_y_continuous(breaks=c(0,0.5,1), limits=c(-0.1,1.1))


F<-ggplot(subset(all_data, Cohort=="G3" & Treatment=="Ctrl"), aes(x=Temp, y=Alive/(Alive+Dead)))+
  geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
  geom_ribbon(data=predict_G3Ctrl,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightskyblue2")+
  geom_line(data=predict_G3Ctrl, aes(x=Temp, y=Prediction))+ggtitle("G3")+
  geom_point(data=subset(all_data, Cohort=="G3" & Treatment=="STHS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral3")+
  geom_ribbon(data=predict_G3STHS,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral2")+
  geom_line(data=predict_G3STHS, aes(x=Temp, y=Prediction))+ylab("")+xlab("")+
  scale_x_continuous(breaks=c(40,40.5,41,41.5,42), labels=c("40°C","40.5°C","41°C","41.5°C","42°C"))+
  theme(axis.text.y= element_blank(),axis.text.x =element_blank())+scale_y_continuous(breaks=c(0,0.5,1), limits=c(-0.1,1.1))

G<-ggplot(subset(all_data, Cohort=="G4" & Treatment=="Ctrl"), aes(x=Temp, y=Alive/(Alive+Dead)))+
  geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
  geom_ribbon(data=predict_G4Ctrl,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightskyblue2")+
  geom_line(data=predict_G4Ctrl, aes(x=Temp, y=Prediction))+ggtitle("G4")+
  geom_point(data=subset(all_data, Cohort=="G4" & Treatment=="STHS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral3")+
  geom_ribbon(data=predict_G4STHS,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral2")+
  geom_line(data=predict_G4STHS, aes(x=Temp, y=Prediction))+ylab("")+xlab("")+
  scale_x_continuous(breaks=c(40,40.5,41,41.5,42), labels=c("40°C","40.5°C","41°C","41.5°C","42°C"))+
  scale_y_continuous(breaks=c(0,0.5,1), limits=c(-0.1,1.1), labels=c("0%", "50%", "100%"))

H<-ggplot(subset(all_data, Cohort=="G5" & Treatment=="Ctrl"), aes(x=Temp, y=Alive/(Alive+Dead)))+
  geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
  geom_ribbon(data=predict_G5Ctrl,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightskyblue2")+
  geom_line(data=predict_G5Ctrl, aes(x=Temp, y=Prediction))+ggtitle("G5")+
  geom_point(data=subset(all_data, Cohort=="G5" & Treatment=="STHS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral3")+
  geom_ribbon(data=predict_G5STHS,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral2")+
  geom_line(data=predict_G5STHS, aes(x=Temp, y=Prediction))+ylab("")+xlab("Temperature (°C)")+
  scale_x_continuous(breaks=c(40,40.5,41,41.5,42), labels=c("40°C","40.5°C","41°C","41.5°C","42°C"))+
  theme(axis.text.y= element_blank())+scale_y_continuous(breaks=c(0,0.5,1), limits=c(-0.1,1.1))

I<-ggplot(subset(all_data, Cohort=="G6" & Treatment=="Ctrl"), aes(x=Temp, y=Alive/(Alive+Dead)))+
  geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
  geom_ribbon(data=predict_G6Ctrl,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightskyblue2")+
  geom_line(data=predict_G6Ctrl, aes(x=Temp, y=Prediction))+ggtitle("G6")+
  geom_point(data=subset(all_data, Cohort=="G6" & Treatment=="STHS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral3")+
  geom_ribbon(data=predict_G6STHS,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral2")+
  geom_line(data=predict_G6STHS, aes(x=Temp, y=Prediction))+ylab("")+xlab("")+
  scale_x_continuous(breaks=c(40,40.5,41,41.5,42), labels=c("40°C","40.5°C","41°C","41.5°C","42°C"))+
  theme(axis.text.y= element_blank())+scale_y_continuous(breaks=c(0,0.5,1), limits=c(-0.1,1.1))

ggarrange(A,B,C,D,E,F,G,H,I, nrow = 3, ncol = 3, align="h")
ggsave(filename="FigS1.png", width=8, height=6.5, units="in", dpi=300)
