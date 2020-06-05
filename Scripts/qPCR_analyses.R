# Analysis of qPCR data from Nematostella larval experiments 
# Data associated with manuscript XXXX 

########### 
# Set up environment
library(plyr)
library(dplyr)
library(ggplot2)
library(DescTools)
library(ggpubr)
load("qPCR.RData")
# The data frame has indvidual sample replicates by gene/plate and information on the 
# the parental cohort, parental treatment, larval treatment, and timepoint for each sample 
###########

### Subset to just reference genes and compute the geometric mean for each sample
subset(qPCR, Gene=="18S" | Gene=="Actin" | Gene=="L10")%>%group_by(Sample)%>%summarize(Norm_avg=Gmean(N0, method="classic"))->qPCR_ref_gmeans
# Add the reference geometric mean back to the main data frame
qPCR<-left_join(qPCR, qPCR_ref_gmeans, by='Sample')
# Normalize the expression by the reference genes 
qPCR<-mutate(qPCR, N0_ref_norm=N0/Norm_avg)

#Create data frame with the normalized values averaged by sample 
subset(qPCR, Gene=="HSP70" | Gene=="CitrateSynthase" | Gene=="MnSOD2")%>%group_by(Sample, Gene)%>%mutate(Mean_norm_exp=mean(N0_ref_norm))%>%distinct(Sample, Gene, .keep_all = TRUE)->GE_df

#Calculate mean expression by Family, Gene, and Long Larval Cat to then run paired T-test for the different comparisons 
GE_df%>%group_by(LongLarvalCategory,Family, Gene)%>%summarize(Mean_Fam_GE=mean(Mean_norm_exp))->GE_df_fam_mean


# run t-test for Baseline between larvae from STHS parents and Control parents 
#Citrate Synthase
t.test(subset(GE_df_fam_mean, LongLarvalCategory=="CCB" & Gene=="CitrateSynthase")$Mean_Fam_GE, subset(GE_df_fam_mean, LongLarvalCategory=="HCB" & Gene=="CitrateSynthase")$Mean_Fam_GE, paired = TRUE)
# p = 0.295

#HSP70
t.test(subset(GE_df_fam_mean, LongLarvalCategory=="CCB" & Gene=="HSP70")$Mean_Fam_GE, subset(GE_df_fam_mean, LongLarvalCategory=="HCB" & Gene=="HSP70")$Mean_Fam_GE, paired = TRUE)
#p = 0.31

#MnSOD2
t.test(subset(GE_df_fam_mean, LongLarvalCategory=="CCB" & Gene=="MnSOD2")$Mean_Fam_GE, subset(GE_df_fam_mean, LongLarvalCategory=="HCB" & Gene=="MnSOD2")$Mean_Fam_GE, paired = TRUE)
#p = 0.25


#Fig4B
(B4<-ggplot(subset(GE_df_fam_mean,LongLarvalCategory=="CCB" | LongLarvalCategory=="HCB"), 
       aes(x=LongLarvalCategory, y=Mean_Fam_GE))+geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
  geom_jitter(aes(colour=Family),height = 0, width=0.2)+
  facet_wrap(~Gene, nrow=1, ncol=3, shrink=TRUE)+
  scale_colour_manual(values=c("#D55E00", "#CC79A7","seagreen", "#999999"))+
  scale_fill_manual(values=c("lightskyblue2", "coral3"), labels=c("Controls", "STHS"))+
  theme_minimal()+guides(alpha=FALSE)+
  ylab("Starting Concentration (flourescence units)")+xlab("")+
  labs(fill="Parents")+
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        panel.grid.major.x = element_blank(),
              legend.background = element_rect(colour="gray"),
              legend.text=element_text(size=10, colour="black", face="bold"),
              legend.title=element_text(size=10, colour="black", face="bold"),
              axis.text.x = element_blank(),
              axis.text.y=element_text(size=8, colour="black", face="bold"),
              axis.title.y=element_text(size=10, colour="black", face="bold"),
        strip.background = element_rect(colour="grey", fill="grey96"),
        strip.text=element_text(size=10, colour="black", face="bold")))

ggsave(filename="/Users/hannyrivera/Documents/Stella_Paper/Manuscript/Figures/Fig4/GE_Baseline.png", 
       width=7, height=5, units="in", dpi=300)



#### Post Heat Stress Timepoint

# run T-test for Post stress between larvae from STHS parents and Control parents 
#Citrate Synthase
t.test(subset(GE_df_fam_mean, LongLarvalCategory=="CHP" & Gene=="CitrateSynthase")$Mean_Fam_GE, subset(GE_df_fam_mean, LongLarvalCategory=="HHP" & Gene=="CitrateSynthase")$Mean_Fam_GE, paired = TRUE)
# p = 0.68

t.test(subset(GE_df_fam_mean, LongLarvalCategory=="CHP" & Gene=="HSP70")$Mean_Fam_GE, subset(GE_df_fam_mean, LongLarvalCategory=="HHP" & Gene=="HSP70")$Mean_Fam_GE, paired = TRUE)
#p = 0.052

t.test(subset(GE_df_fam_mean, LongLarvalCategory=="CHP" & Gene=="MnSOD2")$Mean_Fam_GE, subset(GE_df_fam_mean, LongLarvalCategory=="HHP" & Gene=="MnSOD2")$Mean_Fam_GE, paired = TRUE)
#p = 0.12

# Make plot
# Fig4C
(C4<-ggplot(subset(GE_df_fam_mean,LongLarvalCategory=="CHP" | LongLarvalCategory=="HHP"), 
       aes(x=LongLarvalCategory, y=Mean_Fam_GE))+geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
  geom_jitter(aes(colour=Family),height = 0, width=0.2)+
  facet_wrap(~Gene, nrow=1, ncol=3)+
  scale_colour_manual(values=c("#D55E00", "#CC79A7","seagreen", "#999999"))+
  scale_fill_manual(values=c("lightskyblue2", "coral3"), labels=c("Controls", "STHS"))+
  theme_minimal()+guides(alpha=FALSE)+
  ylab("Starting Concentration (flourescence units)")+xlab("")+
  labs(fill="Parents")+
  coord_cartesian(ylim = c(0,0.05))+
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        panel.grid.major.x = element_blank(),
        legend.background = element_rect(colour="gray"),
        legend.text=element_text(size=10, colour="black", face="bold"),
        legend.title=element_text(size=10, colour="black", face="bold"),
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=8, colour="black", face="bold"),
        axis.title.y=element_text(size=10, colour="black", face="bold"),
        strip.background = element_rect(colour="grey", fill="grey96"),
        strip.text=element_text(size=10, colour="black", face="bold")))

ggsave(filename="/Users/hannyrivera/Documents/Stella_Paper/Manuscript/Figures/Fig4/GE_PostStress.png", 
       width=7, height=5, units="in", dpi=300)

#arrange both
ggarrange(B4, C4, common.legend = TRUE, legend="bottom", nrow=2, ncol=1)

ggsave(filename="/Users/hannyrivera/Documents/Stella_Paper/Manuscript/Figures/Fig4/Fig4BC.png", 
       width=7, height=10, units="in", dpi=300)



### Supplemental Figures 
# run t-test for Baseline between larvae from Ctrl parents and Immediately after heat stress from CTRL parents 
#Citrate Synthase
t.test(subset(GE_df_fam_mean, LongLarvalCategory=="CCB" & Gene=="CitrateSynthase")$Mean_Fam_GE, subset(GE_df_fam_mean, LongLarvalCategory=="CHI" & Gene=="CitrateSynthase")$Mean_Fam_GE, paired = TRUE)
# p = 0.03

#HSP70
t.test(subset(GE_df_fam_mean, LongLarvalCategory=="CCB" & Gene=="HSP70")$Mean_Fam_GE, subset(GE_df_fam_mean, LongLarvalCategory=="CHI" & Gene=="HSP70")$Mean_Fam_GE, paired = TRUE)
#p = 0.04

#MnSOD2
t.test(subset(GE_df_fam_mean, LongLarvalCategory=="CCB" & Gene=="MnSOD2")$Mean_Fam_GE, subset(GE_df_fam_mean, LongLarvalCategory=="CHI" & Gene=="MnSOD2")$Mean_Fam_GE, paired = TRUE)
#p = 0.84


(S4<-ggplot(subset(GE_df_fam_mean,LongLarvalCategory=="CCB" | LongLarvalCategory=="CHI"), 
            aes(x=LongLarvalCategory, y=Mean_Fam_GE))+geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
    geom_jitter(aes(colour=Family),height = 0, width=0.2)+
    facet_wrap(~Gene, nrow=1, ncol=3, scales="free")+
    scale_colour_manual(values=c("#D55E00", "#CC79A7","seagreen", "#999999"))+
    scale_fill_manual(values=c("lightskyblue2", "dodgerblue3"), labels=c("Baseline", "HS-Immediate"))+
    theme_minimal()+guides(alpha=FALSE)+
    ylab("Starting Concentration (flourescence units)")+xlab("")+
    labs(fill="Timepoint")+
    theme(panel.grid.minor = element_blank(), legend.position = "bottom",
          panel.grid.major.x = element_blank(),
          legend.background = element_rect(colour="gray"),
          legend.text=element_text(size=8, colour="black", face="bold"),
          legend.title=element_text(size=8, colour="black", face="bold"),
          axis.text.x = element_blank(),
          axis.text.y=element_text(size=6, colour="black", face="bold"),
          axis.title.y=element_text(size=8, colour="black", face="bold"),
          strip.background = element_rect(colour="grey", fill="grey96"),
          strip.text=element_text(size=8, colour="black", face="bold")))

ggsave(filename="/Users/hannyrivera/Documents/Stella_Paper/Manuscript/Figures/FigS4/FigS4.png", 
       width=6, height=3, units="in", dpi=300)

## Making each plot individually to control axes better 

(S4A<-ggplot(subset(GE_df_2, Family=="P1" & Gene=="CitrateSynthase"), aes(x=factor(LongLarvalCategory,level=c("CCB", "HCB","CHI", "HHI","CHP","HHP")),
                                                         y=Mean_norm_exp))+
    geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
    geom_jitter(aes(colour=LongLarvalCategory),height = 0, width=0.2)+
    scale_colour_manual(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
      values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"))+
    scale_x_discrete(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"))+
    scale_fill_manual(values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"), 
                      limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
                      labels=c("Ctrl-Base","STHS-Base","Ctrl-HS-Ac", "STHS-HS-Ac", "Ctrl-HS-Post", "STHS-HS-Post"))+
    theme_minimal()+
    xlab("")+ylab("")+guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    scale_y_continuous(limits=c(0.001,0.06), breaks=c(0,0.02,0.04,0.06))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.background = element_rect(colour="gray"),
          legend.text = element_text(size=10, colour="black", face="bold"),
          legend.title= element_text(size=10, colour="black", face="bold"),
          axis.text.y = element_text(size=8, colour="black", face="bold"),
          axis.text.x = element_blank(), 
          axis.title.y = element_text(size=10, colour="black", face="bold"),
          strip.background = element_rect(colour="grey", fill="grey96"),
          strip.text=element_blank()))


(S4B<-ggplot(subset(GE_df_2, Family=="P1" & Gene=="HSP70"), aes(x=factor(LongLarvalCategory,level=c("CCB", "HCB","CHI", "HHI","CHP","HHP")),
                                                                          y=Mean_norm_exp))+
    geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
    geom_jitter(aes(colour=LongLarvalCategory),height = 0, width=0.2)+
    scale_colour_manual(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
      values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"))+
    scale_x_discrete(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"))+
    scale_fill_manual(values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"), 
                      limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
                      labels=c("Ctrl-Base","STHS-Base","Ctrl-HS-Ac", "STHS-HS-Ac", "Ctrl-HS-Post", "STHS-HS-Post"))+
    theme_minimal()+
    xlab("")+ylab("")+guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    scale_y_continuous(limits=c(0.001,0.35), breaks=c(0,0.12,0.24,0.35))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.background = element_rect(colour="gray"),
          legend.text = element_text(size=10, colour="black", face="bold"),
          legend.title= element_text(size=10, colour="black", face="bold"),
          axis.text.y = element_text(size=8, colour="black", face="bold"),
          axis.text.x = element_blank(), 
          axis.title.y = element_text(size=10, colour="black", face="bold"),
          strip.background = element_rect(colour="grey", fill="grey96"),
          strip.text=element_blank()))

(S4C<-ggplot(subset(GE_df_2, Family=="P1" & Gene=="MnSOD2"), aes(x=factor(LongLarvalCategory,level=c("CCB", "HCB","CHI", "HHI","CHP","HHP")),
                                                                y=Mean_norm_exp))+
    geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
    geom_jitter(aes(colour=LongLarvalCategory),height = 0, width=0.2)+
    scale_colour_manual(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
      values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"))+
    scale_x_discrete(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"))+
    scale_fill_manual(values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"), 
                      limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
                      labels=c("Ctrl-Base","STHS-Base","Ctrl-HS-Ac", "STHS-HS-Ac", "Ctrl-HS-Post", "STHS-HS-Post"))+
    theme_minimal()+
    xlab("")+ylab("")+guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    scale_y_continuous(limits=c(0.001,0.2), breaks=c(0,0.06,0.13,0.2))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.background = element_rect(colour="gray"),
          legend.text = element_text(size=10, colour="black", face="bold"),
          legend.title= element_text(size=10, colour="black", face="bold"),
          axis.text.y = element_text(size=8, colour="black", face="bold"),
          axis.text.x = element_blank(), 
          axis.title.y = element_text(size=10, colour="black", face="bold"),
          strip.background = element_rect(colour="grey", fill="grey96"),
          strip.text=element_blank()))

(S4D<-ggplot(subset(GE_df_2, Family=="P2" & Gene=="CitrateSynthase"), aes(x=factor(LongLarvalCategory,level=c("CCB", "HCB","CHI", "HHI","CHP","HHP")),
                                                                 y=Mean_norm_exp))+
    geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
    geom_jitter(aes(colour=LongLarvalCategory),height = 0, width=0.2)+
    scale_colour_manual(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
      values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"))+
    scale_x_discrete(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"))+
    scale_fill_manual(values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"), 
                      limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
                      labels=c("Ctrl-Base","STHS-Base","Ctrl-HS-Ac", "STHS-HS-Ac", "Ctrl-HS-Post", "STHS-HS-Post"))+
    theme_minimal()+
    xlab("")+ylab("")+guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    scale_y_continuous(limits=c(0,0.08), breaks=c(0,0.026,0.052,0.08))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.background = element_rect(colour="gray"),
          legend.text = element_text(size=10, colour="black", face="bold"),
          legend.title= element_text(size=10, colour="black", face="bold"),
          axis.text.y = element_text(size=8, colour="black", face="bold"),
          axis.text.x = element_blank(), 
          axis.title.y = element_text(size=10, colour="black", face="bold"),
          strip.background = element_rect(colour="grey", fill="grey96"),
          strip.text=element_blank()))


(S4E<-ggplot(subset(GE_df_2, Family=="P2" & Gene=="HSP70"), aes(x=factor(LongLarvalCategory,level=c("CCB", "HCB","CHI", "HHI","CHP","HHP")),
                                                                          y=Mean_norm_exp))+
    geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
    geom_jitter(aes(colour=LongLarvalCategory),height = 0, width=0.2)+
    scale_colour_manual(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
      values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"))+
    scale_x_discrete(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"))+
    scale_fill_manual(values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"), 
                      limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
                      labels=c("Ctrl-Base","STHS-Base","Ctrl-HS-Ac", "STHS-HS-Ac", "Ctrl-HS-Post", "STHS-HS-Post"))+
    theme_minimal()+
    xlab("")+ylab("")+guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    scale_y_continuous(limits=c(0,0.3), breaks=c(0,0.1,0.2,0.3))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.background = element_rect(colour="gray"),
          legend.text = element_text(size=10, colour="black", face="bold"),
          legend.title= element_text(size=10, colour="black", face="bold"),
          axis.text.y = element_text(size=8, colour="black", face="bold"),
          axis.text.x = element_blank(), 
          axis.title.y = element_text(size=10, colour="black", face="bold"),
          strip.background = element_rect(colour="grey", fill="grey96"),
          strip.text=element_blank()))

(S4F<-ggplot(subset(GE_df_2, Family=="P2" & Gene=="MnSOD2"), aes(x=factor(LongLarvalCategory,level=c("CCB", "HCB","CHI", "HHI","CHP","HHP")),
                                                                y=Mean_norm_exp))+
    geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
    geom_jitter(aes(colour=LongLarvalCategory),height = 0, width=0.2)+
    scale_colour_manual(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
      values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"))+
    scale_x_discrete(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"))+
    scale_fill_manual(values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"), 
                      limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
                      labels=c("Ctrl-Base","STHS-Base","Ctrl-HS-Ac", "STHS-HS-Ac", "Ctrl-HS-Post", "STHS-HS-Post"))+
    theme_minimal()+
    xlab("")+ylab("")+guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    scale_y_continuous(limits=c(0,0.3), breaks=c(0,0.1,0.2,0.3))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.background = element_rect(colour="gray"),
          legend.text = element_text(size=10, colour="black", face="bold"),
          legend.title= element_text(size=10, colour="black", face="bold"),
          axis.text.y = element_text(size=8, colour="black", face="bold"),
          axis.text.x = element_blank(), 
          axis.title.y = element_text(size=10, colour="black", face="bold"),
          strip.background = element_rect(colour="grey", fill="grey96"),
          strip.text=element_blank()))

(S4G<-ggplot(subset(GE_df_2, Family=="P3" & Gene=="CitrateSynthase"), aes(x=factor(LongLarvalCategory,level=c("CCB", "HCB","CHI", "HHI","CHP","HHP")),
                                                                 y=Mean_norm_exp))+
    geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
    geom_jitter(aes(colour=LongLarvalCategory),height = 0, width=0.2)+
    scale_colour_manual(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
      values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"))+
    scale_x_discrete(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"))+
    scale_fill_manual(values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"), 
                      limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
                      labels=c("Ctrl-Base","STHS-Base","Ctrl-HS-Ac", "STHS-HS-Ac", "Ctrl-HS-Post", "STHS-HS-Post"))+
    theme_minimal()+
    xlab("")+ylab("")+guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    scale_y_continuous(limits=c(0.001,0.05), breaks=c(0,0.017,0.034,0.05))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.background = element_rect(colour="gray"),
          legend.text = element_text(size=10, colour="black", face="bold"),
          legend.title= element_text(size=10, colour="black", face="bold"),
          axis.text.y = element_text(size=8, colour="black", face="bold"),
          axis.text.x = element_blank(), 
          axis.title.y = element_text(size=10, colour="black", face="bold"),
          strip.background = element_rect(colour="grey", fill="grey96"),
          strip.text=element_blank()))

(S4H<-ggplot(subset(GE_df_2, Family=="P3" & Gene=="HSP70"), aes(x=factor(LongLarvalCategory,level=c("CCB", "HCB","CHI", "HHI","CHP","HHP")),
                                                                          y=Mean_norm_exp))+
    geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
    geom_jitter(aes(colour=LongLarvalCategory),height = 0, width=0.2)+
    scale_colour_manual(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
      values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"))+
    scale_x_discrete(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"))+
    scale_fill_manual(values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"), 
                      limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
                      labels=c("Ctrl-Base","STHS-Base","Ctrl-HS-Ac", "STHS-HS-Ac", "Ctrl-HS-Post", "STHS-HS-Post"))+
    theme_minimal()+
    xlab("")+ylab("")+guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    scale_y_continuous(limits=c(0.001,0.12), breaks=c(0,0.04,0.08,0.12))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.background = element_rect(colour="gray"),
          legend.text = element_text(size=10, colour="black", face="bold"),
          legend.title= element_text(size=10, colour="black", face="bold"),
          axis.text.y = element_text(size=8, colour="black", face="bold"),
          axis.text.x = element_blank(), 
          axis.title.y = element_text(size=10, colour="black", face="bold"),
          strip.background = element_rect(colour="grey", fill="grey96"),
          strip.text=element_blank()))

(S4I<-ggplot(subset(GE_df_2, Family=="P3" & Gene=="MnSOD2"), aes(x=factor(LongLarvalCategory,level=c("CCB", "HCB","CHI", "HHI","CHP","HHP")),
                                                                y=Mean_norm_exp))+
    geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
    geom_jitter(aes(colour=LongLarvalCategory),height = 0, width=0.2)+
    scale_colour_manual(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
      values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"))+
    scale_x_discrete(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"))+
    scale_fill_manual(values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"), 
                      limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
                      labels=c("Ctrl-Base","STHS-Base","Ctrl-HS-Ac", "STHS-HS-Ac", "Ctrl-HS-Post", "STHS-HS-Post"))+
    theme_minimal()+
    xlab("")+ylab("")+guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    scale_y_continuous(limits=c(0.001,0.05), breaks=c(0,0.017,0.034,0.05))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.background = element_rect(colour="gray"),
          legend.text = element_text(size=10, colour="black", face="bold"),
          legend.title= element_text(size=10, colour="black", face="bold"),
          axis.text.y = element_text(size=8, colour="black", face="bold"),
          axis.text.x = element_blank(), 
          axis.title.y = element_text(size=10, colour="black", face="bold"),
          strip.background = element_rect(colour="grey", fill="grey96"),
          strip.text=element_blank()))


(S4J<-ggplot(subset(GE_df_2, Family=="P4" & Gene=="CitrateSynthase"), aes(x=factor(LongLarvalCategory,level=c("CCB", "HCB","CHI", "HHI","CHP","HHP")),
                                                                 y=Mean_norm_exp))+
    geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
    geom_jitter(aes(colour=LongLarvalCategory),height = 0, width=0.2)+
    scale_colour_manual(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
                        values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"))+
    scale_x_discrete(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
                     labels=c("Ctrl-Base","STHS-Base","Ctrl-HS-Im", "STHS-HS-Im", "Ctrl-HS-Post", "STHS-HS-Post"))+
    scale_fill_manual(values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"), 
                      limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"))+
    theme_minimal()+
    xlab("")+ylab("")+guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    scale_y_continuous(limits=c(0.001,0.03), breaks=c(0,0.01,0.02,0.03))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.background = element_rect(colour="gray"),
          legend.text = element_text(size=10, colour="black", face="bold"),
          legend.title= element_text(size=10, colour="black", face="bold"),
          axis.text.y = element_text(size=8, colour="black", face="bold"),
          axis.text.x = element_text(size=8, colour="black", face="bold"),
          axis.title.y = element_text(size=10, colour="black", face="bold"),
          strip.background = element_rect(colour="grey", fill="grey96"),
          strip.text=element_blank()))

(S4K<-ggplot(subset(GE_df_2, Family=="P4" & Gene=="HSP70"), aes(x=factor(LongLarvalCategory,level=c("CCB", "HCB","CHI", "HHI","CHP","HHP")),
                                                                 y=Mean_norm_exp))+
    geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
    geom_jitter(aes(colour=LongLarvalCategory),height = 0, width=0.2)+
    scale_colour_manual(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
                        values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"))+
    scale_x_discrete(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"), 
                     labels=c("Ctrl-Base","STHS-Base","Ctrl-HS-Im", "STHS-HS-Im", "Ctrl-HS-Post", "STHS-HS-Post"))+
    scale_fill_manual(values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"), 
                      limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"))+    
    theme_minimal()+
    xlab("")+ylab("")+guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    scale_y_continuous(limits=c(0.001,0.075), breaks=c(0,0.025,0.05,0.075))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.background = element_rect(colour="gray"),
          legend.text = element_text(size=10, colour="black", face="bold"),
          legend.title= element_text(size=10, colour="black", face="bold"),
          axis.text.y = element_text(size=8, colour="black", face="bold"),
          axis.text.x = element_text(size=8, colour="black", face="bold"), 
          axis.title.y = element_text(size=10, colour="black", face="bold"),
          strip.background = element_rect(colour="grey", fill="grey96"),
          strip.text=element_blank()))

(S4L<-ggplot(subset(GE_df_2, Family=="P4" & Gene=="MnSOD2"), aes(x=factor(LongLarvalCategory,level=c("CCB", "HCB","CHI", "HHI","CHP","HHP")),
                                                                          y=Mean_norm_exp))+
    geom_boxplot(aes(fill=LongLarvalCategory,alpha=0.2), outlier.shape = NA)+
    geom_jitter(aes(colour=LongLarvalCategory),height = 0, width=0.2)+
    scale_colour_manual(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"),
                        values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"))+
    scale_x_discrete(limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"), 
                     labels=c("Ctrl-Base","STHS-Base","Ctrl-HS-Im", "STHS-HS-Im", "Ctrl-HS-Post", "STHS-HS-Post"))+
    scale_fill_manual(values=c("lightskyblue2","coral3","lightskyblue2","coral3", "lightskyblue2","coral3"), 
                      limits=c("CCB", "HCB", "CHI", "HHI", "CHP", "HHP"))+
    theme_minimal()+
    xlab("")+ylab("")+guides(fill=FALSE, colour=FALSE, alpha=FALSE)+
    scale_y_continuous(limits=c(0.001,0.03), breaks=c(0,0.01,0.02,0.03))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.background = element_rect(colour="gray"),
          legend.text = element_text(size=10, colour="black", face="bold"),
          legend.title= element_text(size=10, colour="black", face="bold"),
          axis.text.y = element_text(size=8, colour="black", face="bold"),
          axis.text.x = element_text(size=8, colour="black", face="bold"),
          axis.title.y = element_text(size=10, colour="black", face="bold"),
          strip.background = element_rect(colour="grey", fill="grey96"),
          strip.text=element_blank()))

ggarrange(S4A, S4B,S4C,ncol=3, nrow=1, align="h")
ggarrange(S4D, S4E,S4F,ncol=3, nrow=1, align="h")
ggarrange(S4G, S4H,S4I,ncol=3, nrow=1, align="h")
ggarrange(S4J, S4K,S4L,ncol=3, nrow=1, align="hv")

ggarrange(S4A, S4B,S4C, S4D, S4E, S4F, S4G, S4H, S4I, S4J, S4K, S4L, ncol=3, nrow=4, align="hv")

ggsave(filename="/Users/hannyrivera/Documents/Stella_Paper/Manuscript/Figures/FigS5/FigS5.png", 
       width=20, height=10, units="in", dpi=300)

