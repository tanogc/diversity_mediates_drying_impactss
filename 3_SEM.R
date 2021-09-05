################################################################################################################
#                                                                                                               
# R script to run SEM to explore direct and indirect impacts of drying aspects on invertebrate density                                                                   
#                                                                                                               
# Diversity mediates the responses of invertebrate density to annual drying duration and frequency              
#                                                                                                               
# Rebeca Arias-Real, Cayetano Gutiérrez-Cánovas, Margarita Menéndez, Verónica Granados 
# & Isabel Muñoz
#                                                                                                               
# Code written by Rebeca Arias-Real and Cayetano Gutiérrez-Cánovas                                              
# email for queries: rebeca.arias.real@gmail.com or cayeguti@um.es                                              
################################################################################################################

# working folder
setwd("/your_folder")
# Setting results plots
plot_folder<-paste(getwd(),"/plots/",sep="")
assum_folder<-paste(getwd(),"/plots/assum/",sep="")

# loading packages
library(MuMIn)
library(plyr)
library(vegan)
library(usdm)
library(viridis)
library(piecewiseSEM)
library(reshape)

source("FD_functions.R")

#loading data
metric<-read.table("final_set.txt",h=T,sep="\t", dec=".")
env<-read.table("env.txt",h=T,sep="\t", dec=".")

# Arranging dataframes
hydro<-env[,2:5]
env<-env[,-c(2:5)]

# Arranging matrices
row.names(metric)<-metric$site

# Merging all matrices
dat2<-dat<-data.frame(metric, hydro, alt=env$alt, do=env$do)

# Variable transformation
q<-sapply(dat,class)=="numeric" | sapply(dat,class)=="integer"# selecting quantitative variables

# Assessing variable distributions

par(mfrow=c(3,3))
for (i in which(q==T)) hist(dat[,i], main=names(dat)[i])
par(mfrow=c(1,1))

sqrt(dat$ric)->dat$ric
dat$ZFP<-sqrt(dat$ZFP)
dat$RE<-(dat$ZFL^(1/4))
dat$RE<-(dat$RE^(1/4))
log(dat$abun)->dat$abun
log(dat$FR)->dat$FR
log(dat$FR.PRE+1)->dat$FR.PRE
log(dat$FR.GRA)->dat$FR.GRA
log(dat$FR.GAT)->dat$FR.GAT
log(dat$FR.OMNI)->dat$FR.OMNI
log(dat$FR.SHR+1)->dat$FR.SHR
log(dat$FR.FIL+1)->dat$FR.FIL

# Checking problems after transformation
apply(dat,2,function(x) any(is.infinite(x)==T))

# Standardising hydrological and non-hydrological env. metrics
dat<-data.frame(scale(dat))

# Variance Inflation Factor (VIF)
vifstep(dat[, c("alt","do","ZFT", "ric")])
vifstep(dat[, c("ZFT", "ric")])
vifstep(dat[, c("alt","do","ZFT", "meanRD")])
vifstep(dat[, c("ZFT", "meanRD")])

# Performing SEM using piecewiseSEM

# Simple linear regression models
summary(lm(abun~ric, dat))
summary(lm(abun~meanRD, dat))

# Competing SEM
sem1 <- summary(mod1<-psem(lm(meanRD ~ ZFT + do, dat), lm(abun ~ ZFT, dat)))
sem2 <- summary(mod2<-psem(lm(meanRD ~ ZFT + do, dat), lm(abun ~ meanRD, dat)))
sem3 <- summary(mod3<-psem(lm(meanRD ~ ZFT + do, dat), lm(abun ~ meanRD + ZFT, dat)))
sem4 <- summary(mod4<-psem(lm(meanRD ~ ZFP + do, dat), lm(abun ~ ZFT, dat)))
sem5 <- summary(mod5<-psem(lm(meanRD ~ ZFP + do, dat), lm(abun ~ meanRD, dat)))
sem6 <- summary(mod6<-psem(lm(meanRD ~ ZFP + do, dat), lm(abun ~ meanRD + ZFT, dat)))
sem7 <- summary(mod7<-psem(lm(meanRD ~ ZFP + ZFT + do, dat), lm(abun ~ meanRD + ZFT, dat)))

sem8 <- summary(mod8<-psem(lm(ric ~ ZFT + do, dat), lm(abun ~ ZFT, dat)))
sem9 <- summary(mod9<-psem(lm(ric ~ ZFT + do, dat), lm(abun ~ ric, dat)))
sem10 <- summary(mod10<-psem(lm(ric ~ ZFT + do, dat), lm(abun ~ ric + ZFT, dat)))
sem11 <- summary(mod11<-psem(lm(ric ~ ZFP + do, dat), lm(abun ~ ZFT, dat)))
sem12 <- summary(mod12<-psem(lm(ric ~ ZFP + do, dat), lm(abun ~ ric, dat)))
sem13 <- summary(mod13<-psem(lm(ric ~ ZFP + do, dat), lm(abun ~ ric + ZFT, dat)))
sem14 <- summary(mod14<-psem(lm(ric ~ ZFP + ZFT + do, dat), lm(abun ~ ric + ZFT, dat)))

# AICc values
AIC_sem<-c(sem1$IC$AICc, sem2$IC$AICc, sem3$IC$AICc, sem4$IC$AICc, sem5$IC$AICc,
           sem6$IC$AICc, sem7$IC$AICc, sem8$IC$AICc, sem9$IC$AICc, sem10$IC$AICc,
           sem11$IC$AICc, sem12$IC$AICc, sem13$IC$AICc, sem14$IC$AICc)

# Model results
sem_res<-data.frame(mod=1:14, rbind(sem1$C, sem2$C, sem3$C, sem4$C, sem5$C,sem6$C, sem7$C, sem8$C, sem9$C, sem10$C,
                 sem11$C, sem12$C, sem13$C, sem14$C), AICc=round(AIC_sem,2),w=round(Weights(AIC_sem),2),
           r2=round(rbind(rsquared(mod1)$R.squared,rsquared(mod2)$R.squared,rsquared(mod3)$R.squared,
                 rsquared(mod4)$R.squared,rsquared(mod5)$R.squared,rsquared(mod6)$R.squared,
                 rsquared(mod7)$R.squared,rsquared(mod8)$R.squared,rsquared(mod9)$R.squared,
                 rsquared(mod10)$R.squared,rsquared(mod11)$R.squared,rsquared(mod12)$R.squared,
                 rsquared(mod13)$R.squared,rsquared(mod14)$R.squared),2))

# Model names
models=c(rep("mod1", 3), rep("mod2", 3), rep("mod3", 4), rep("mod4",3), rep("mod5", 3), rep("mod6",4), rep("mod7",5),
         rep("mod8", 3), rep("mod9", 3), rep("mod10", 4), rep("mod11",3), rep("mod12", 3), rep("mod13",4), rep("mod14",5))

# Model coefficients
coef_sem<-data.frame(mod=models, 
           rbind(sem1$coef, sem2$coef, sem3$coef, sem4$coef, sem5$coef,
           sem6$coef, sem7$coef, sem8$coef, sem9$coef, sem10$coef,
           sem11$coef, sem12$coef, sem13$coef, sem14$coef))

element <- as.numeric(lapply(list(sem3$coef, sem5$coef,sem10$coef, sem12$coef),nrow))

w2 <- sem_res$w[c(3,5,10,12)] / sum(sem_res$w[c(3,5,10,12)])

coef_sem2<-data.frame(mod=c(rep("mod3", 4), rep("mod5",3), rep("mod10", 4), rep("mod12",3)),
                      w=c(rep(w2[1], 4), rep(w2[2],3), rep(w2[3], 4), rep(w2[4],3)),
                          rbind(sem3$coef, sem5$coef,sem10$coef, sem12$coef))

coef_m1=ddply(coef_sem2[which(coef_sem2$Response=="abun"),],.(Predictor),function(x) mean(x$Estimate))
coef_se1=ddply(coef_sem2[which(coef_sem2$Response=="abun"),],.(Predictor),function(x) sd(x$Estimate)/sqrt(length(x)))

# Mean and SE
coef_m<-coef_m1$V1
names(coef_m)<-c("FunDiv", "T. richness", "ZFT")

coef_se<-coef_se1$V1
names(coef_se)<-c("FD", "TRic", "ZFT")

coef_ci1<-coef_m - qnorm(0.975)*as.numeric(coef_se)
coef_ci2<-coef_m + qnorm(0.975)*as.numeric(coef_se)

#####

r2_m<-100*c(rsquared(mod1)[2,5],mean(sem_res$r2.2[c(3,5,10,12)]))
r2_se<-c(0, sd(100*sem_res$r2.2[c(3,5,10,12)])/sqrt(length(sem_res$r2.2[c(3,5,10,12)])))

r2_ci1<-r2_m - qnorm(0.975)*as.numeric(r2_se)
r2_ci2<-r2_m + qnorm(0.975)*as.numeric(r2_se)

# Setting colours for each hydrologic metric
col.hyd<-c("#56B4E9", "#7AD151FF","brown1", "gold3","dark green","green")

pdf(file=paste(plot_folder,"sem_coef.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=4.5,height=4)

dotplot(coef_m, xlab="Standardized coefficient", xlim=c(-0.70, 0.70),
        
        par.settings = list(axis.text = list(cex = 1.2, font=1), 
                            par.xlab.text = list(cex = 1.2, font=1), 
                            par.ylab.text = list(cex = 2, font=1)),
        panel=function(x,y){
          panel.segments(coef_ci1, as.numeric(y), coef_ci2, as.numeric(y), lty=1, col=col.hyd[c(5,6,1)])
          panel.xyplot(x, y, pch=15, cex=2,col=col.hyd[c(5,6,1)])
          panel.abline(v=0, col="black", lty=2)
          
        })

dev.off()

pdf(file=paste(plot_folder,"sem_r2.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=8.6,height=4)

par(mfrow=c(1,2), cex.lab=1.4, cex.axis=1.3, mar=c(4,5,4,4))

plot.new()
mtext("b)", line = 2, adj = -0.2, cex = 1.4, font = 2)

foo <- barplot(r2_m, ylim=c(0,55), main="", col=col.hyd,ylab="Explained variance (%)", names=c("ZFT","Diversity"))
arrows(x0=foo,y0=r2_ci2,y1=r2_ci1,angle=90,code=3,length=0.1)

dev.off()

# Ploting SEM results
plot(mod3)

# Saving results
write.table(sem_res,paste (plot_folder,"sem_res.txt", sep=""),sep="\t",row.names = F)
write.table(coef_sem,paste (plot_folder,"sem_coef.txt", sep=""),sep="\t",row.names = F)

