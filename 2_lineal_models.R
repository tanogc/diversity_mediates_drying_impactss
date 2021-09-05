################################################################################################################
#                                                                                                               
# R script to run models evaluating the effect and importance of drying aspects on inverebrate communities                                                                    
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
library(variancePartition)
library(plyr)
library(vegan)
library(usdm)
library(viridis)

source("0_FD_functions.R")

#loading data
metric<-read.table("final_set.txt",h=T,sep="\t", dec=".")
env<-read.table("env.txt",h=T,sep="\t", dec=".")

# Arranging dataframes
hydro<-env[,2:5]
env<-env[,-c(2:5)]

# Arranging matrices
row.names(metric)<-metric$site

# Mean correlations between environmental features and biodiversity metrics
env.cor<-round(cor(env[,-c(1,10)], data.frame(log(metric$abun),metric$ric,metric$meanRD),method="pearson"),2)
env.cor
write.table(env.cor, paste (plot_folder, "env_cor.txt", sep=""),sep="\t")

# Relationships among environmental features
as.dist(round(cor(env[,-c(1,10)]),2))
round(cor(env[,-c(1,10)], hydro),2)

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

# Creating final set of predictors
env<-dat[, c("alt","do", "ZFL","ZFP","ZFT", "RE")]

# Collinearity

# Variance Inflation Factor (VIF)
vifstep(dat[, c("alt","do","ZFT")])
vifstep(dat[, c("alt","do","ZFP")])
vifstep(dat[, c("alt","do","ZFL")])
vifstep(dat[, c("alt","do","RE")])
vifstep(dat[, c("alt","do","ZFT", "ZFP")])
vifstep(dat[, c("alt","do","ZFT", "RE")])
vifstep(dat[, c("alt","do","ZFP", "ZFL")])
vifstep(dat[, c("alt","do","ZFP", "RE")])
vifstep(dat[, c("alt","do","ZFL", "RE")])

# Pairwise Pearson correlations
round(as.dist(cor(env, method=c("pearson"))),2)
round(as.dist(cor(env [-c(2)], method=c("pearson"))),2)

sel.var<-names(dat)[c(1:2,4:16)] # selecting biodiversity metrics

############ Loop to run models and perform variance partitioning

r2_order<-res_av<-res_d<-res_d_full<-res<-list() # to store data

a=0

for (i in sel.var) {
  
  a=a+1
  
  mod1 <- lm(dat[,i] ~ ZFL + alt + do, data=dat)
  
  mod2 <- lm(dat[,i] ~ ZFP + alt + do, data=dat)
  
  mod3 <- lm(dat[,i] ~ RE + alt + do, data=dat)
  
  mod4 <- lm(dat[,i] ~ ZFT + alt + do, data=dat)
  
  mod5 <- lm(dat[,i] ~ RE + ZFT + alt + do, data=dat)
  
  mod6 <- lm(dat[,i] ~ ZFP + ZFT + alt + do, data=dat)
  
  mod7 <- lm(dat[,i] ~ ZFL + RE + alt + do, data=dat)
  
  mod8 <- lm(dat[,i] ~ ZFP + RE + alt + do, data=dat)
  
  mod9 <- lm(dat[,i] ~ ZFL + ZFP + alt + do, data=dat)
  

  mod.list<-list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8, mod9)
  
  res_d_full[[a]]<-mod_d<-model.sel(mod.list, extra=c(r2=function(x) r.squaredGLMM(x)[1]))
  
  mod_d[which(mod_d$delta<=2), ]->res_d[[a]]
  
  mod_set <- get.models (mod_d, subset=delta<=2) # subset 
 
  round(hier.r2.res(mod_set, sel.pred=c("alt","do","RE","ZFL","ZFP","ZFT"))[rownames(res_d[[a]]),],3)->res[[a]]
  
  mod_d[as.character(1:nrow(mod_d)),]->r2_order[[a]]
  
  #as.numeric(lapply(mod.list, moran_num))->moran_I[[a]]
  
  # Saving residual plots
  
  length(mod_set)->n.mod
  
  for (k in 1:n.mod) 
  {
    mod_set[[k]]->mod
    pdf(file=paste(assum_folder,"mod_int_",i,"_",k,"_assum_int",".pdf",sep=""),onefile=T,width=6.24,height=4.05)
    par(mfrow=c(1,2),cex=1, cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
    m.resid<-resid(mod)
    hist(m.resid,main="",xlab="residuals")
    mtext(paste("model_int",k),line=0.5,at=4,cex=2)
    plot(fitted(mod),m.resid,xlab="fitted values",ylab="residuals",main="")  
    dev.off()
  }
  
}

# Reshaping and reorganising results
res.df<-data.frame(var=rep(sel.var, unlist(lapply(res,nrow))), do.call(rbind.data.frame, res))
res_d.df<-data.frame(var=rep(sel.var, unlist(lapply(res_d,nrow))), do.call(rbind.data.frame, res_d))
res_d_full.df<-data.frame(var=rep(sel.var, unlist(lapply(res_d_full,nrow))), do.call(rbind.data.frame, res_d_full))
res_d.ord<-data.frame(var=rep(sel.var, unlist(lapply(r2_order,nrow))), do.call(rbind.data.frame, r2_order))

# Saving (exporting) model results
write.table(data.frame(var=res.df[,1],round(res.df[,-1], 3), res_d.df), paste(plot_folder,"varparticion_sep20.txt", sep=""),sep="\t", row.names = F)
write.table(res_d.df,paste (plot_folder, "multimodel_sep20.txt", sep=""),sep="\t")
write.table(res_d_full.df,paste (plot_folder,"multimodel_full_sep20.txt", sep=""),sep="\t")

# Weighted determination coefficients (r2)
res.df[,-1] * res_d.df$weight->w_r2

r2.full<-ddply(data.frame(var=res.df[,1], w_r2), .(var), function(x) round(colSums(x[,-1]),3))
r2.dat<-r2.full[,c(7,6,5,4)]

# Saving weighted-averaged r2 values
write.table(r2.full, paste(plot_folder,"avg_r2.txt", sep=""),sep="\t", row.names = F)

# Mean r2
abun_r2<-100*colMeans(r2.dat[c(1:7),])
div_r2<-100*colMeans(r2.dat[c(8:15),])

# SE r2
abun_se<-100*unlist(apply(r2.dat[c(1:7),], 2, function(x) if(is.na(sd(x,na.rm=T))==T) 0 else sd(x,na.rm=T)/sqrt(length(x))))
div_se<-100*unlist(apply(r2.dat[c(8:15),], 2, function(x) if(is.na(sd(x,na.rm=T))==T) 0 else sd(x,na.rm=T)/sqrt(length(x))))

# 95% CI for std. coefficients

abun_ci1<-abun_r2 - qnorm(0.975)*as.numeric(abun_se)
abun_ci2<-abun_r2 + qnorm(0.975)*as.numeric(abun_se)
div_ci1<-div_r2 - qnorm(0.975)*as.numeric(div_se)
div_ci2<-div_r2 + qnorm(0.975)*as.numeric(div_se)

# Setting colours for each hydrologic metric
col.hyd<-c("#56B4E9", "#7AD151FF","brown1", "gold3")

# Plotting barplots

pdf(file=paste(plot_folder,"r2_barplots2.pdf",sep=""),onefile=T,width=8.6,height=3.6)

par(mfrow=c(1,2), cex.lab=1.4, cex.axis=1.3, mar=c(4,5,4,4))

foo <- barplot(abun_r2, ylim=c(0,33), main="Abundance", col=col.hyd,ylab="Explained variance (%)",xlab="")
arrows(x0=foo,y0=abun_ci2,y1=abun_ci1,angle=90,code=3,length=0.1)
mtext("a)", line = 2, adj = -0.2, cex = 1.4, font = 2)

foo <- barplot(div_r2, ylim=c(0,33), main="Diversity", col=col.hyd,ylab="",xlab="")
arrows(x0=foo,y0=div_ci2,y1=div_ci1,angle=90,code=3,length=0.1)
mtext("b)", line = 2, adj = -0.2, cex = 1.4, font = 2)

dev.off()

# Mean standardized coefficients

# All coefficients
mod_coef_all<-data.frame(var=res_d.df$var, res_d.df[,c(3,4,8,6,5,7)]* res_d.df$weight)
avg_mod_coef_all<-ddply(mod_coef_all,.(var), function(x) colSums(x[,-1],na.rm=T))
# Saving weighted-averaged r2 values
write.table(avg_mod_coef_all, paste(plot_folder,"avg_coef.txt", sep=""),sep="\t", row.names = F)

# Hydrological coefficients
mod_coef<-data.frame(var=res_d.df$var, res_d.df[,c(8,6,5,7)]* res_d.df$weight)
abiotic_coef<-data.frame(var=res_d.df$var, res_d.df[,c(3,4)]* res_d.df$weight)

# Replacing NAs by zeros
for (i in 1:nrow(mod_coef)) mod_coef[i, which(is.na(mod_coef[i,])==T)]<-0

avg_mod_coef<-ddply(mod_coef,.(var), function(x) colSums(x[,-1],na.rm=T))
abi_mod_coef<-ddply(abiotic_coef,.(var), function(x) colSums(x[,-1],na.rm=T))

# DO and altitude coefficients and SE
colMeans(abi_mod_coef[,-1])
apply(abi_mod_coef[,-1],2,function(x) sd(x)/sqrt(length(x)))

abun_coef <- colMeans(avg_mod_coef[1:7, 5:2])
div_coef <- colMeans(avg_mod_coef[c(8:15), 5:2])

# SE for std. coefficients

abun_se<-unlist(apply(res_d.df[c(4:5,7:18),c(8,6,5,7)], 2, function(x) if(is.na(sd(x,na.rm=T))==T) 0 else sd(x,na.rm=T)/sqrt(length(x))))[4:1]
div_se<-unlist(apply(res_d.df[-c(4:5,7:18),c(8,6,5,7)] , 2, function(x) if(is.na(sd(x,na.rm=T))==T) 0 else sd(x,na.rm=T)/sqrt(length(x))))[4:1]

# 95% CI for std. coefficients

abun_ci1<-abun_coef - qnorm(0.975)*as.numeric(abun_se)
abun_ci2<-abun_coef + qnorm(0.975)*as.numeric(abun_se)
div_ci1<-div_coef - qnorm(0.975)*as.numeric(div_se)
div_ci2<-div_coef + qnorm(0.975)*as.numeric(div_se)

# Multitrophic abundance

pdf(file=paste(plot_folder,"abun_coef2.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=5.5)

dotplot(abun_coef, xlab="Standardized coefficient", xlim=c(-0.55, 0.35),
        
        par.settings = list(axis.text = list(cex = 2, font=1), 
                            par.xlab.text = list(cex = 2, font=1), 
                            par.ylab.text = list(cex = 2, font=1)),
        panel=function(x,y){
          panel.segments(abun_ci1, as.numeric(y), abun_ci2, as.numeric(y), lty=1, col=col.hyd[4:1])
          panel.xyplot(x, y, pch=15, cex=2,col=col.hyd[4:1])
          panel.abline(v=0, col="black", lty=2)
          
        })

dev.off()

# Diversity-type

pdf(file=paste(plot_folder,"div_coef.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=5.5)

dotplot(div_coef, xlab="Standardized coefficient", xlim=c(-0.55, 0.35),
        
        par.settings = list(axis.text = list(cex = 2, font=1), 
                            par.xlab.text = list(cex = 2, font=1), 
                            par.ylab.text = list(cex = 2, font=1)),
        panel=function(x,y){
          panel.segments(div_ci1, as.numeric(y), div_ci2, as.numeric(y), lty=1, col=col.hyd[4:1])
          panel.xyplot(x, y, pch=15, cex=2,col=col.hyd[4:1])
          panel.abline(v=0, col="black", lty=2)
          
        })

dev.off()


######################################################
# Functional groups


# Mean r2
sfp_r2<-sum(100*colMeans(r2.dat[c(2,6,7,9,13,14),]))
ggo_r2<-sum(100*colMeans(r2.dat[-c(1,2,6,7,8,9,13,14,15),]))


# SE r2
sfp_se<-sd(100*colMeans(r2.dat[c(2,6,7,9,13,14),]))/sqrt(length(r2.dat[c(2,6,7,9,13,14),]))
ggo_se<-sd(100*colMeans(r2.dat[-c(1,2,6,7,8,9,13,14,15),]))/sqrt(length(r2.dat[-c(1,2,6,7,8,9,13,14,15),]))

# 95% CI for std. coefficients

sfp_ci1<-sfp_r2 - qnorm(0.975)*as.numeric(sfp_se)
sfp_ci2<-sfp_r2 + qnorm(0.975)*as.numeric(sfp_se)
ggo_ci1<-ggo_r2 - qnorm(0.975)*as.numeric(ggo_se)
ggo_ci2<-ggo_r2 + qnorm(0.975)*as.numeric(ggo_se)

trophic_r2<-c(sfp_r2, ggo_r2)
names(trophic_r2)<-c("SHR-FIL-PRE", "GAT-GRA-OMNI")
trophic_ci1<-c(sfp_ci1, ggo_ci1)
trophic_ci2<-c(sfp_ci2, ggo_ci2)

# Setting colours for each hydrologic metric
col.hyd<-c("#56B4E9", "#7AD151FF","brown1", "gold3")

# Plotting barplots

pdf(file=paste(plot_folder,"trophic_r2_barplots2.pdf",sep=""),onefile=T,width=5,height=3.6)

par(mfrow=c(1,1), cex.lab=1.4, cex.axis=1.2, mar=c(4,5,4,4))

foo <- barplot(trophic_r2, ylim=c(0,33), main="", col=viridis(3)[2:3],ylab="Explained variance (%)",xlab="")
arrows(x0=foo,y0=trophic_ci2,y1=trophic_ci1,angle=90,code=3,length=0.1)

dev.off()

# Mean coefficients for trophic groups

sfp_coef <- colMeans(avg_mod_coef[c(2,6,7,9,13,14), 5:2])
ggo_coef <- colMeans(avg_mod_coef[c(3:5,10:13), 5:2])

# SE for std. coefficients for trophic groups

sfp_se<-unlist(apply(res_d.df[-c(1:8,10:11,15:17,19:20,23:28,31:34),c(8,6,5,7)], 2, function(x) if(is.na(sd(x,na.rm=T))==T) 0 else sd(x,na.rm=T)/sqrt(length(x))))[4:1]
ggo_se<-unlist(apply(res_d.df[c(7:8,10:11,15:17,19:20,23:28,31:34),c(8,6,5,7)] , 2, function(x) if(is.na(sd(x,na.rm=T))==T) 0 else sd(x,na.rm=T)/sqrt(length(x))))[4:1]

# 95% CI for std. coefficients

sfp_ci1<-sfp_coef - qnorm(0.975)*as.numeric(sfp_se)
sfp_ci2<-sfp_coef + qnorm(0.975)*as.numeric(sfp_se)
ggo_ci1<-ggo_coef - qnorm(0.975)*as.numeric(ggo_se)
ggo_ci2<-ggo_coef + qnorm(0.975)*as.numeric(ggo_se)

trophic_coef <- c(ggo_coef[3],ggo_coef[4], sfp_coef[3], sfp_coef[4])
names(trophic_coef) <- c("ZFP 2", "ZFT 2", "ZFP 1", "ZFT 1")
  
trophic_ci1 <- c(ggo_ci1[3], ggo_ci1[4],sfp_ci1[3], sfp_ci1[4])
trophic_ci2 <- c(ggo_ci2[3], ggo_ci2[4],sfp_ci2[3], sfp_ci2[4])

# Diversity-type

pdf(file=paste(plot_folder,"trophic_coef.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=6,height=5.5)

dotplot(trophic_coef, xlab="Standardized coefficient", xlim=c(-0.5, 0.4),
        
        par.settings = list(axis.text = list(cex = 2, font=1), 
                            par.xlab.text = list(cex = 2, font=1), 
                            par.ylab.text = list(cex = 2, font=1)),
        panel=function(x,y){
          panel.segments(trophic_ci1, as.numeric(y), trophic_ci2, as.numeric(y), lty=1, col=col.hyd[c(2,1,2,1)])
          panel.xyplot(x, y, pch=15, cex=2,col=col.hyd[c(2,1,2,1)])
          panel.abline(v=0, col="black", lty=2)
          
        })

dev.off()
