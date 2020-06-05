##### Analyses of Adult Short Term Heat Stress Mortality Experiment 

### Data
data<-matrix(c(64,73,16,7), nrow = 2)
dimnames(data) <- list("Treatment" = c("CTRL","STHS"), "Mortality" = c("Alive","Dead"))

fisher.test(data, alternative="less")
