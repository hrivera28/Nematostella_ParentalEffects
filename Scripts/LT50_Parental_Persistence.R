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
load('Persistence_data.RData')

####################################### 
## Run two models for each cohort 

for (fam in levels(Persist_data$Cohort)){
  # One in which we model different LT 50s for each Treatment
  assign(paste(fam,"_P.LTDiff", sep=""), drm(Alive/(Alive+Dead) ~ Temp,Treatment, weights = (Alive+Dead), data = subset(Persist_data, Cohort==fam), fct = LL.2(), type = "binomial"))
  # One in which we model the LT50s as being the same for both Treatments
  assign(paste(fam,"_P.LTSame", sep=""), drm(Alive/(Alive+Dead) ~ Temp,Treatment, weights = (Alive+Dead), data = subset(Persist_data, Cohort==fam), fct = LL.2(), type = "binomial", pmodels = list(~factor(Treatment)-1, ~1)))
}

summary(P1_P.LTDiff)
summary(P1_P.LTSame)
anova(P1_P.LTDiff,P1_P.LTSame)
## For P1 we conclude the LT50 values ARE the same between treatments (p=0.24)

summary(P2_P.LTDiff)
summary(P2_P.LTSame)
anova(P2_P.LTDiff,P2_P.LTSame)
# For P2 we conclude the LT50 values are NOT the same between treatments (p=0.0496)
# Difference between treament values is due to the larvae from the STHS parents having a lower LT50 now

summary(P4_P.LTDiff)
summary(P4_P.LTSame)
anova(P4_P.LTDiff,P4_P.LTSame)
## For P4 we conclude the LT50 values ARE the same between treatments (p=0.9579)


# Test for differences between the groups using paired t-tests 
CtrlP_drm<-drm(Alive/(Alive+Dead) ~ Temp, Cohort, weights = (Alive+Dead), data = subset(Persist_data, Treatment=="Ctrl"), fct = LL.2(), type = "binomial")
STHSP_drm<-drm(Alive/(Alive+Dead) ~ Temp, Cohort,weights = (Alive+Dead), data = subset(Persist_data, Treatment=="STHS"),  fct = LL.2(), type = "binomial")

#Get simulatneous standard errors and confidence intervals
CtrlP_SI<-summary(glht(CtrlP_drm))
STHSP_SI<-summary(glht(STHSP_drm))
# The ED estimates live in CtrlP_SI$test$coefficients[4:6]
# The SE measures live in CtrlP_SI$test$sigma[4:6]

# Make a dataframe with the LT50 Estimates and SE by cohort and treatment
LT_mat_P<-as.data.frame(matrix(ncol=4, nrow=6))
colnames(LT_mat_P)<-c("LT50", "Cohort", "Treatment", "SE")

for (i in seq(1,3)){
  LT_mat_P$LT50[i]<-CtrlP_SI$test$coefficients[i+3] #LT50 estimate for Control
  LT_mat_P$LT50[i+3]<-STHSP_SI$test$coefficients[i+3] #LT50 estimate for STHS
  
  LT_mat_P$SE[i]<-CtrlP_SI$test$sigma[i+3] #SE estimate for Control
  LT_mat_P$SE[i+3]<-STHSP_SI$test$sigma[i+3] #SE estimate for Control
  
  LT_mat_P$Treatment[i]<-"Ctrl"
  LT_mat_P$Treatment[i+3]<-"STHS"
}

LT_mat_P$Cohort<-rep(c("P4","P2", "P1"),2)

# Paired t-test of the data 
t.test(LT50~Treatment, data=LT_mat_P, paired=TRUE)
# p = 0.8446 not sig different  
# estimated mean difference is 0.02 C 


### Plotting results
cons_palette<-c("#D55E00", "#CC79A7", "#999999") # keep colors consistent across parental pops

(F1<-ggplot(LT_mat_P, aes(x=Treatment, y=LT50, group=Cohort, colour=Cohort))+
  geom_point(position=position_dodge(width = 0.1))+
  geom_line(position=position_dodge(width = 0.1))+
  geom_errorbar(aes(ymin=LT50-SE, ymax=LT50+SE),width=0, position=position_dodge(width = 0.1))+
  theme_minimal()+theme(panel.grid.minor = element_blank(),
                        legend.position = "bottom", 
                        legend.background = element_rect(size=0.25), 
                        axis.text = element_text(size=8, colour="black", face="bold"))+
  scale_y_continuous(breaks=c(40.8,41,41.2,41.4,41.6,41.8), limits=c(40.75,41.8), labels=NULL)+
  xlab("")+ylab("")+
  scale_x_discrete(labels=c("Controls", "STHS"), expand=c(0.1,0.1))+
  scale_color_manual(values = cons_palette)+guides(colour=FALSE))

ggsave("Fig1F.png", width=3, height=3, units="in", dpi=300)
