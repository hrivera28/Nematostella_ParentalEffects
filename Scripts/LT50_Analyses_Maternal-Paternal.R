### Analysis of maternal/paternal effects  
### Script associated with manuscript Rivera et al.
### Plasticity in parental effects confers rapid thermal tolerance in N. vectensis larvae
### Copyright H.E. Rivera 

##### 
# Set up environment
library(drc)
library(ggplot2)
library(ggeffects)
library(multcomp)
library(lme4)
library(stargazer)
library(ggpubr)
library(plyr)
library(dplyr)
library(rcompanion)
load("Mat_pat_data.RData")

###########################################
# Running different survival models for each cohort 

for (fam in levels(All_data$Cohort)){
    assign(paste(fam,".LTDiff_all", sep=""),drm(Alive/(Alive+Dead) ~ Temp,Treatment, weights = (Alive+Dead), data = subset(All_data, Cohort==fam), fct = LL.2(), type = "binomial"))
    assign(paste(fam,"_SI", sep=""), summary(glht(get(paste(fam,".LTDiff_all", sep="")))))
}

# Making matrix with LT50 values and SE for linear model 
# Paramenters here have already been corrected for simulatneous inference
# The ED estimates live in "fam"_SI$test$coefficients[5:8]
# The SE measures live in "fam"_SI$test$sigma[5:8]

LT_mat_MD<-as.data.frame(matrix(ncol=4, nrow=24))
colnames(LT_mat_MD)<-c("LT50", "Cohort", "Treatment", "SE")

i=0
for (fam in levels(droplevels(subset(All_data, Cohort!="G2"))$Cohort)){
  for (x in seq(1,4)){
      LT_mat_MD$LT50[x+i]<-get(paste(fam,"_SI", sep=""))$test$coefficients[x+4] #LT50 estimate for Control
      LT_mat_MD$SE[x+i]<-get(paste(fam,"_SI", sep=""))$test$sigma[x+4] #SE estimate for Control
      LT_mat_MD$Cohort[x+i]<-fam
  }
  i=i+4
}
LT_mat_MD$Treatment<-rep(c("HS","CT", "HS-Mom", "HS-Dad"), 6)

# Manually coding in G2 since that one is missing the HS-mom data (bad spawning)
LT_mat_MD$LT50[21]<-G2_SI$test$coefficients[4]
LT_mat_MD$SE[21]<-G2_SI$test$sigma[4]

LT_mat_MD$LT50[22]<-G2_SI$test$coefficients[5]
LT_mat_MD$SE[22]<-G2_SI$test$sigma[5]

LT_mat_MD$LT50[24]<-G2_SI$test$coefficients[6]
LT_mat_MD$SE[24]<-G2_SI$test$sigma[6]

LT_mat_MD$Cohort[20:24]<-"G2"


#### Running linear model for the effect of treamtment with parental family (cohort) as a random effect
effects_model<-lmer(LT50~Treatment+(1|Cohort), data=LT_mat_MD, na.action=na.omit)
summary(effects_model)

#Get model results in table form 
stargazer(effects_model, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")

cbbPalette <- c("#000000", "#E69F00", "#5956A0", "#009E73", "#F0E442", "#0072B2")

# Adding in dummy variable for regression plotting
LT_mat_MD$Treatment2<-LT_mat_MD$Treatment
LT_mat_MD$Treatment2[LT_mat_MD$Treatment2=="CT"]<-1
LT_mat_MD$Treatment2[LT_mat_MD$Treatment2=="HS-Dad"]<-2
LT_mat_MD$Treatment2[LT_mat_MD$Treatment2=="HS-Mom"]<-3
LT_mat_MD$Treatment2[LT_mat_MD$Treatment2=="HS"]<-4
LT_mat_MD$Treatment2<-as.numeric(LT_mat_MD$Treatment2)

(C2<-ggplot(LT_mat_MD, aes(x=Treatment2, y=LT50))+
  geom_point(data=LT_mat_MD, aes(x=Treatment2, y=LT50, colour=Cohort))+
  stat_smooth(method="lm", colour="slategrey")+
  theme_minimal()+
  scale_y_continuous(breaks=c(40.8,41,41.2,41.4,41.6), limits=c(40.75,41.72), labels=c("40.8","41.0","41.2","41.4","41.6"))+
  theme(axis.text=element_text(size=8, face="bold", colour="black"),
        axis.title=element_text(size=10, face="bold", colour="black"),
        legend.position = "bottom",
        legend.background = element_rect(colour="grey"),
        legend.text=element_text(size=6, colour="black", face="bold"),
        legend.title=element_text(size=8, colour="black", face="bold"),
        panel.grid.minor = element_blank())+
  labs(colour="Cohort")+xlab("")+ylab("LT50 (°C)")+
  guides(colour = guide_legend(nrow = 1))+
  scale_color_manual(values=cbbPalette)+
  scale_x_continuous(breaks=c(1,2,3,4), labels=c("Controls", "Ctrl Mom x\nSTHS Dad","STHS Mom x\nCtrl Dad", "STHS"))+
  annotate("text", label="*", x=3, y=41.6, size=9)+
  annotate("text", label="**", x=4, y=41.7, size=9))

ggsave("Fig2C.png", width = 5, height = 4, units= "in", dpi=300)


# Fig 2
(B2<-ggplot(subset(All_data, Cohort=="G3" &  Treatment=="CT"), aes(x=Temp, y=Alive/(Alive+Dead)))+
  geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
  geom_ribbon(data=predict_G3_CT_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightblue3")+
  geom_line(data=predict_G3_CT_all, aes(x=Temp, y=Prediction),colour="lightblue4")+
  geom_point(data=subset(All_data, Cohort=="G3" &  Treatment=="HS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral4")+
  geom_ribbon(data=predict_G3_HS_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral3")+
  geom_line(data=predict_G3_HS_all, aes(x=Temp, y=Prediction), colour="coral4")+
  geom_point(data=subset(All_data, Cohort=="G3" &  Treatment=="HS-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="plum4")+
  geom_ribbon(data=predict_G3_HS_fem,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="plum3")+
  geom_line(data=predict_G3_HS_fem, aes(x=Temp, y=Prediction), colour="plum4")+
  geom_point(data=subset(All_data, Cohort=="G3" &  Treatment=="CT-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="darkolivegreen4")+
  geom_ribbon(data=predict_G3_CT_fem,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="darkolivegreen3")+
  geom_line(data=predict_G3_CT_fem, aes(x=Temp, y=Prediction), colour="darkolivegreen4")+
  ylab("Survival")+xlab("Temperature (°C)")+
  scale_y_continuous(breaks = c(0,0.5,1),limits=c(-0.1,1.05), labels=c("0%", "50%", "100%"))+
  scale_x_continuous(breaks=c(40,40.5,41,41.5,42))+
  theme(axis.text=element_text(size=8, face="bold", colour="black"),
        axis.title=element_text(size=10, face="bold", colour="black"),
        panel.grid.minor = element_blank()))

ggsave("Fig2B.png", width = 5, height =4, units= "in", dpi=300)


###########################################
# Plotting the survival surves for each cohort
## For these graphs make sure the prediction is only on one of the treatments/trials 
## and then plot only those point results else you get weird looking graphs


for (fam in levels(All_data$Cohort)){
  # Same model as previous but use results to calucate interval for ribbon plot
  name1<-paste(fam,".CT", sep="")
  assign(paste(fam,".CT", sep=""), drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(All_data, Cohort==fam & Treatment=="CT"), fct = LL.2(), type = "binomial"))
  assign(paste("predict_", fam,"_CT_all", sep=""), predict(get(name1), newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence"))
  name2<-paste("predict_", fam,"_CT_all", sep="")
  assign(paste("predict_", fam,"_CT_all", sep=""), cbind(get(name2),expand.grid(Temp=seq(39.8, 42.3, length=25))))
  
  # For the STHS treatment
  name3<-paste(fam,".HS", sep="")
  assign(paste(fam,".HS", sep=""),drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(All_data, Cohort==fam & Treatment=="HS"),  fct = LL.2(), type = "binomial"))
  assign(paste("predict_", fam,"_HS_all", sep=""), predict(get(name3), newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence"))
  name4<-paste("predict_", fam,"_HS_all", sep="")
  assign(paste("predict_", fam,"_HS_all", sep=""),cbind(get(name4),expand.grid(Temp=seq(39.8, 42.3, length=25))))

  # For the STHS Mother treatment
  if (fam!="G2"){
  name5<-paste(fam,".HSF", sep="")
  assign(paste(fam,".HSF", sep=""),drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(All_data, Cohort==fam & Treatment=="HS-Fem"),  fct = LL.2(), type = "binomial"))
  assign(paste("predict_", fam,"_HS_fem", sep=""), predict(get(name5), newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence"))
  name6<-paste("predict_", fam,"_HS_fem", sep="")
  assign(paste("predict_", fam,"_HS_fem", sep=""),cbind(get(name6),expand.grid(Temp=seq(39.8, 42.3, length=25))))
  }
  
  # For the STHS father treatment
  name7<-paste(fam,".CTF", sep="")
  assign(paste(fam,".CTF", sep=""),drm(Alive/(Alive+Dead) ~ Temp, weights = (Alive+Dead), data = subset(All_data, Cohort==fam & Treatment=="CT-Fem"),  fct = LL.2(), type = "binomial"))
  assign(paste("predict_", fam,"_CT_fem", sep=""), predict(get(name7), newdata=expand.grid(Temp=seq(39.8, 42.3, length=25)), interval="confidence"))
  name8<-paste("predict_", fam,"_CT_fem", sep="")
  assign(paste("predict_", fam,"_CT_fem", sep=""),cbind(get(name8),expand.grid(Temp=seq(39.8, 42.3, length=25))))
  
}


# G1
(G1_plot<-ggplot(subset(All_data, Cohort=="G1" & Treatment=="CT"), aes(x=Temp, y=Alive/(Alive+Dead)))+
    geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
    geom_ribbon(data=predict_G1_CT_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightblue3")+
    geom_line(data=predict_G1_CT_all, aes(x=Temp, y=Prediction),colour="lightblue4")+
    geom_point(data=subset(All_data, Cohort=="G1" &  Treatment=="HS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral4")+
    geom_ribbon(data=predict_G1_HS_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral3")+
    geom_line(data=predict_G1_HS_all, aes(x=Temp, y=Prediction), colour="coral4")+
    geom_point(data=subset(All_data, Cohort=="G1" &  Treatment=="HS-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="plum4")+
    geom_ribbon(data=predict_G1_HS_fem,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="plum3")+
    geom_line(data=predict_G1_HS_fem, aes(x=Temp, y=Prediction), colour="plum4")+
    geom_point(data=subset(All_data, Cohort=="G1" &  Treatment=="CT-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="darkolivegreen4")+
    geom_ribbon(data=predict_G1_CT_fem,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="darkolivegreen3")+
    geom_line(data=predict_G1_CT_fem, aes(x=Temp, y=Prediction), colour="darkolivegreen4")+
    ylab("Survival")+xlab("")+ggtitle("G1")+
    scale_y_continuous(breaks = c(0,0.5,1), limits=c(-0.1,1.25), labels=c("0%", "50%", "100%"))+
    theme(panel.grid.major.x = element_blank(),
          axis.text.x=element_blank(),
          axis.text = element_text(face="bold", colour="black", size=8), 
          axis.title.y = element_text(face="bold", colour="black", size=10),
          plot.title = element_text(face="bold", colour="black", size=12)))

# G2
(G2_plot<-ggplot(subset(All_data, Cohort=="G2" & Treatment=="CT"), aes(x=Temp, y=Alive/(Alive+Dead)))+
    geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
    geom_ribbon(data=predict_G2_CT_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightblue3")+
    geom_line(data=predict_G2_CT_all, aes(x=Temp, y=Prediction),colour="lightblue4")+
    geom_point(data=subset(All_data, Cohort=="G2" & Treatment=="HS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral4")+
    geom_ribbon(data=predict_G2_HS_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral3")+
    geom_line(data=predict_G2_HS_all, aes(x=Temp, y=Prediction), colour="coral4")+
    geom_point(data=subset(All_data, Cohort=="G2" & Treatment=="CT-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="darkolivegreen4")+
    geom_ribbon(data=predict_G2_CT_fem,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="darkolivegreen3")+
    geom_line(data=predict_G2_CT_fem, aes(x=Temp, y=Prediction), colour="darkolivegreen4")+
    ylab("")+xlab("")+ggtitle("G2")+
    scale_y_continuous(breaks = c(0,0.5,1),limits=c(-0.1,1.25))+
    theme(panel.grid.major.x = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(face="bold", colour="black", size=12)))

(G3_plot<-ggplot(subset(All_data, Cohort=="G3" &  Treatment=="CT"), aes(x=Temp, y=Alive/(Alive+Dead)))+
    geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
    geom_ribbon(data=predict_G3_CT_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightblue3")+
    geom_line(data=predict_G3_CT_all, aes(x=Temp, y=Prediction),colour="lightblue4")+
    geom_point(data=subset(All_data, Cohort=="G3" &  Treatment=="HS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral4")+
    geom_ribbon(data=predict_G3_HS_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral3")+
    geom_line(data=predict_G3_HS_all, aes(x=Temp, y=Prediction), colour="coral4")+
    geom_point(data=subset(All_data, Cohort=="G3" &  Treatment=="HS-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="plum4")+
    geom_ribbon(data=predict_G3_HS_fem,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="plum3")+
    geom_line(data=predict_G3_HS_fem, aes(x=Temp, y=Prediction), colour="plum4")+
    geom_point(data=subset(All_data, Cohort=="G3" &  Treatment=="CT-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="darkolivegreen4")+
    geom_ribbon(data=predict_G3_CT_fem,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="darkolivegreen3")+
    geom_line(data=predict_G3_CT_fem, aes(x=Temp, y=Prediction), colour="darkolivegreen4")+
    ylab("")+xlab("")+ggtitle("G3")+
    scale_y_continuous(breaks = c(0,0.5,1),limits=c(-0.1,1.25))+
    theme(panel.grid.major.x = element_blank(),
          axis.text = element_blank(),
          plot.title = element_text(face="bold", colour="black", size=12)))
    
# Trial 4 
(G4_plot<-ggplot(subset(All_data, Cohort=="G4" &  Treatment=="CT"), aes(x=Temp, y=Alive/(Alive+Dead)))+
  geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
  geom_ribbon(data=predict_G4_CT_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightblue3")+
  geom_line(data=predict_G4_CT_all, aes(x=Temp, y=Prediction),colour="lightblue4")+
  geom_point(data=subset(All_data, Cohort=="G4" &  Treatment=="HS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral4")+
  geom_ribbon(data=predict_G4_HS_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral3")+
  geom_line(data=predict_G4_HS_all, aes(x=Temp, y=Prediction), colour="coral4")+
  geom_point(data=subset(All_data, Cohort=="G4" &  Treatment=="HS-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="plum4")+
  geom_ribbon(data=predict_G4_HS_fem,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="plum3")+
  geom_line(data=predict_G4_HS_fem, aes(x=Temp, y=Prediction), colour="plum4")+
  geom_point(data=subset(All_data, Cohort=="G4" &  Treatment=="CT-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="darkolivegreen4")+
  geom_ribbon(data=predict_G4_CT_fem,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="darkolivegreen3")+
  geom_line(data=predict_G4_CT_fem, aes(x=Temp, y=Prediction), colour="darkolivegreen4")+
  ylab("Survival")+xlab("")+ggtitle("G4")+
  scale_y_continuous(breaks = c(0,0.5,1), limits=c(-0.1,1.25), labels=c("0%", "50%", "100%"))+
  scale_x_continuous(breaks=c(40,40.5,41,41.5,42), labels=c("40","40.5","41","41.5","42"))+
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(face="bold", colour="black", size=8),
        axis.title.y=element_text(face="bold", colour="black", size=10),
        plot.title = element_text(face="bold", colour="black", size=12)))


#Trial 5
(G5_plot<-ggplot(subset(All_data, Cohort=="G5" &  Treatment=="CT"), aes(x=Temp, y=Alive/(Alive+Dead)))+
  geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
  geom_ribbon(data=predict_G5_CT_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightblue3")+
  geom_line(data=predict_G5_CT_all, aes(x=Temp, y=Prediction),colour="lightblue4")+
  geom_point(data=subset(All_data, Cohort=="G5" & Treatment=="HS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral4")+
  geom_ribbon(data=predict_G5_HS_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral3")+
  geom_line(data=predict_G5_HS_all, aes(x=Temp, y=Prediction), colour="coral4")+
  geom_point(data=subset(All_data, Cohort=="G5" & Treatment=="HS-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="plum4")+
  geom_ribbon(data=predict_G5_HS_fem,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="plum3")+
  geom_line(data=predict_G5_HS_fem, aes(x=Temp, y=Prediction), colour="plum4")+
  geom_point(data=subset(All_data, Cohort=="G5" & Treatment=="CT-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="darkolivegreen4")+
  geom_ribbon(data=predict_G5_CT_fem,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="darkolivegreen3")+
  geom_line(data=predict_G5_CT_fem, aes(x=Temp, y=Prediction), colour="darkolivegreen3")+
  ylab("")+xlab("Temp (°C)")+ggtitle("G5")+
  scale_y_continuous(breaks = c(0,0.5,1),limits=c(-0.1,1.25))+
  scale_x_continuous(breaks=c(40,40.5,41,41.5,42), labels=c("40","40.5","41","41","42"))+
  theme(panel.grid.major.x = element_blank(),
        axis.text = element_text(face="bold", colour="black", size=8),
        axis.text.y= element_blank(),
        axis.title.x=element_text(face="bold", colour="black", size=10),
        plot.title = element_text(face="bold", colour="black", size=12)))

#Trial 6

(G6_plot<-ggplot(subset(All_data, Cohort=="G6" & Treatment=="CT"), aes(x=Temp, y=Alive/(Alive+Dead)))+
    geom_point(colour="lightblue4")+theme_minimal()+theme(panel.grid.minor = element_blank())+
    geom_ribbon(data=predict_G6_CT_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="lightblue3")+
    geom_line(data=predict_G6_CT_all, aes(x=Temp, y=Prediction),colour="lightblue4")+
    geom_point(data=subset(All_data, Cohort=="G6" & Treatment=="HS"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="coral4")+
    geom_ribbon(data=predict_G6_HS_all,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="coral3")+
    geom_line(data=predict_G6_HS_all, aes(x=Temp, y=Prediction), colour="coral4")+
    geom_point(data=subset(All_data, Cohort=="G6" & Treatment=="HS-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="plum4")+
    geom_ribbon(data=predict_G6_HS_fem,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="plum3")+
    geom_line(data=predict_G6_HS_fem, aes(x=Temp, y=Prediction), colour="plum4")+
    geom_point(data=subset(All_data, Cohort=="G6" & Treatment=="CT-Fem"), aes(x=Temp, y=Alive/(Alive+Dead)),colour="darkolivegreen4")+
    geom_ribbon(data=predict_G6_CT_fem,aes(x=Temp, y=Prediction, ymin=Lower, ymax=Upper), alpha=0.5, fill="darkolivegreen3")+
    geom_line(data=predict_G6_CT_fem, aes(x=Temp, y=Prediction), colour="darkolivegreen4")+
    ylab("")+xlab("")+ggtitle("G6")+
    scale_y_continuous(breaks = c(0,0.5,1), limits=c(-0.1,1.25))+
    scale_x_continuous(breaks=c(40,40.5,41,41.5,42))+
    theme(panel.grid.major.x = element_blank(),
          axis.text.y = element_blank(), 
          axis.text.x = element_text(face="bold", size=8),
          plot.title = element_text(face="bold", colour="black", size=12)))

ggarrange(G1_plot, G2_plot, G3_plot, G4_plot, G5_plot, G6_plot, nrow=2, ncol=3)

ggsave("FigS2.png", width=10, height = 6, units="in", dpi=300)

#Plot to grab legend from
ggplot(subset(All_data, Cohort=="G3"), aes(x=Temp, y=Alive/(Alive+Dead), colour=Treatment))+
  geom_point()+theme_minimal()+theme(panel.grid.minor = element_blank())+
  scale_colour_manual(values=c("lightblue3", "darkolivegreen3", "plum3", " coral3"),
                      limits=c("CT", "CT-Fem", "HS-Fem", "HS"),
                      labels=c("Controls", "Ctrl Mom x STHS Dad","STHS Mom x Ctrl Dad", "STHS"))+
  labs(colour="Treatment")+
  theme(legend.text = element_text(size=8, colour="black", face="bold"), 
        legend.position = "bottom", 
        legend.title=element_text(face="bold", colour="black", size=10),
        legend.background = element_rect(colour="grey"))

ggsave("Fig2_leg.png", width=7, height = 4, units="in", dpi=300)





