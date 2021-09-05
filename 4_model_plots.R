################################################################################################################
#                                                                                                               
# R script to plot invertebrate metric responses to drying aspects                                                              
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

gr.name <- c("FIL","GAT","GRA","PRE","SHR")

sel.var <- c("abun", "ric","meanRD",
             "FR.FIL","FR.GAT","FR.GRA","FR.PRE",
             "RD.FIL","RD.PRE","RD.SHR")


sel.pred<-c("alt","do","NF", "ZFP","RE","ZFT")

gr.col<-c(blues9[5], "#4600FFFF","#00FF66FF","#FF0000FF","#CC00FFFF")

res_d<-mod_set_list<-list()

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
  
  
  mod.list<-list(mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8,mod9)
  
  # Evaluating models
  
  mod_d<-model.sel(mod.list, extra=c(r2=function(x) r.squaredGLMM(x)[1]))
  
  mod_d[which(mod_d$delta<=2), ]->res_d[[a]]
  
  mod_set_list[[a]] <- get.models (mod_d, subset=delta<=2) # subset
  
}

days.seq<-seq(0, 300, length.out = 1000)
per.seq<-seq(0, 8, length.out = 1000)

# Response plots

pdf(file=paste(plot_folder,"metric_res2.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=9.5,height=6)

par(mfrow=c(2,3), cex.lab=1.8, cex.axis=1.6, mar=c(4,5,4,4))

###################
# Abundance
###################

plot(dat2$abun~dat2$ZFT, xlim=c(0,300), ylim=c(0,max(dat2$abun)), main="", ylab=expression(Total ~ density ~ (ind /m ^{2})), xlab="")

mtext("a)", line = 1.9, adj = -0.3, cex = 1.6, font = 2)

# Plotting fitted values

mod_set  <- mod_set_list[[1]]

mod.weights <- Weights(as.numeric(lapply(mod_set, AICc)))

mod.pred <- data.frame(lapply(mod_set, function(x) { 
  
  predict(x,data.frame(ZFT=days.seq,
                       ZFP=mean(dat$ZFP),
                       ZFL=mean(dat$ZFL),
                       RE=mean(dat$RE),
                       do=mean(dat$do),
                       alt=mean(dat$alt)))
}))

# mean weighted prediction
mean_pred <- exp(apply(mod.pred, 1, function(x) sum(x * mod.weights)))
if (var(mean_pred)>0) lines(days.seq, mean_pred, col="blue", lwd=3)

###################
# Taxonomic richness
###################

plot(dat2$ric~dat2$ZFT, xlim=c(0,300), main="", ylab="Taxon richness", xlab="Zero Flow days (ZFT)")

#mtext("b)", line = 1.0, adj = -0.2, cex = 1.6, font = 2)

# Plotting fitted values

mod_set  <- mod_set_list[[2]]

mod.weights <- Weights(as.numeric(lapply(mod_set, AICc)))

mod.pred <- data.frame(lapply(mod_set, function(x) { 
  
  predict(x,data.frame(ZFT=days.seq,
                       ZFP=mean(dat$ZFP),
                       ZFL=mean(dat$ZFL),
                       RE=mean(dat$RE),
                       do=mean(dat$do),
                       alt=mean(dat$alt)))
}))

# mean weighted prediction
mean_pred <- apply(mod.pred, 1, function(x) sum(x * mod.weights))^2
if (var(mean_pred)>0) lines(days.seq, mean_pred, col="blue", lwd=3)


###################
# Functional diversity
###################

plot(dat2$meanRD~dat2$ZFT, xlim=c(0,300), main="", ylab="Functional diversity", xlab="")

#mtext("c)", line = 1.0, adj = -0.2, cex = 1.6, font = 2)

# Plotting fitted values

mod_set  <- mod_set_list[[3]]

mod.weights <- Weights(as.numeric(lapply(mod_set, AICc)))

mod.pred <- data.frame(lapply(mod_set, function(x) { 
  
  predict(x,data.frame(ZFT=days.seq,
                       ZFP=mean(dat$ZFP),
                       ZFL=mean(dat$ZFL),
                       RE=mean(dat$RE),
                       do=mean(dat$do),
                       alt=mean(dat$alt)))
}))

# mean weighted prediction
mean_pred <- apply(mod.pred, 1, function(x) sum(x * mod.weights))
if (var(mean_pred)>0) lines(per.seq, mean_pred, col="blue", lwd=3)

#################
# Abundance
#################

plot(dat2$abun~dat2$ZFP, xlim=c(0,8), ylim=c(0,max(dat2$abun)), main="", ylab=expression(Total ~ density ~ (ind /m ^{2})), xlab="")

mtext("b)", line = 1.9, adj = -0.3, cex = 1.6, font = 2)

# Plotting fitted values

mod_set  <- mod_set_list[[1]]

mod.weights <- Weights(as.numeric(lapply(mod_set, AICc)))

mod.pred <- data.frame(lapply(mod_set, function(x) { 
  
  predict(x,data.frame(ZFT=mean(dat$ZFT),
                       ZFP=sqrt(per.seq),
                       ZFL=mean(dat$ZFL),
                       RE=mean(dat$RE),
                       do=mean(dat$do),
                       alt=mean(dat$alt)))
}))

# mean weighted prediction
mean_pred <- exp(apply(mod.pred, 1, function(x) sum(x * mod.weights)))
if (var(mean_pred)>0) lines(days.seq, mean_pred, col="blue", lwd=3)


####################
# Taxon richness
####################

plot(dat2$ric~dat2$ZFP, xlim=c(0,8), main="", ylab="Taxon richness", xlab="Zero Flow Periods (ZFP)")

#mtext("f)", line = 1.0, adj = -0.2, cex = 1.6, font = 2)

# Plotting fitted values

mod_set  <- mod_set_list[[2]]

mod.weights <- Weights(as.numeric(lapply(mod_set, AICc)))

mod.pred <- data.frame(lapply(mod_set, function(x) { 
  
  predict(x,data.frame(ZFT=mean(dat$ZFT),
                       ZFP=sqrt(per.seq),
                       ZFL=mean(dat$ZFL),
                       RE=mean(dat$RE),
                       do=mean(dat$do),
                       alt=mean(dat$alt)))
}))

# mean weighted prediction
mean_pred <- apply(mod.pred, 1, function(x) sum(x * mod.weights))^2
if (var(mean_pred)>0) lines(per.seq, mean_pred, col="blue", lwd=3)



###################
# Functional diversity
###################

plot(dat$meanRD~dat2$ZFP, xlim=c(0,8), main="", ylab="Functional diversity", xlab="")

# Plotting fitted values

mod_set  <- mod_set_list[[3]]

mod.weights <- Weights(as.numeric(lapply(mod_set, AICc)))

mod.pred <- data.frame(lapply(mod_set, function(x) { 
  
  predict(x,data.frame(ZFT=mean(dat$ZFT),
                       ZFP=sqrt(per.seq),
                       ZFL=mean(dat$ZFL),
                       RE=mean(dat$RE),
                       do=mean(dat$do),
                       alt=mean(dat$alt)))
}))

# mean weighted prediction
mean_pred <- apply(mod.pred, 1, function(x) sum(x * mod.weights))
if (var(mean_pred)>0) lines(per.seq, mean_pred, col="blue", lwd=3)

dev.off()

##########################
# Trophic groups
##########################


pdf(file=paste(plot_folder,"metric_res3.pdf",sep=""),,useDingbats=FALSE,onefile=T,width=9,height=2.8)

par(mfrow=c(1,3), cex.lab=1.5, cex.axis=1.4,mar=c(4,5,4,4))

###################
# FR groups
###################

plot(NULL, xlim=c(0,300), ylim=c(0,9000), main="", ylab=expression(Density ~ (ind /m ^{2})), xlab="Zero Flow days (ZFT)")

# Plotting fitted values

a=0

for (j in 4:7) {
  
  a=a+1
  
  mod_set  <- mod_set_list[[j]]
  
  mod.weights <- Weights(as.numeric(lapply(mod_set, AICc)))
  
  mod.pred <- data.frame(lapply(mod_set, function(x) { 
    
    predict(x,data.frame(ZFT=days.seq,
                         ZFP=mean(dat$ZFP),
                         ZFL=mean(dat$ZFL),
                         RE=mean(dat$RE),
                         do=mean(dat$do),
                         alt=mean(dat$alt)))
  }))
  
  # mean weighted prediction
  mean_pred <- exp(apply(mod.pred, 1, function(x) sum(x * mod.weights)))
  if (var(mean_pred)>0) lines(days.seq, mean_pred, col=gr.col[a], lwd=3)
  
}

legend("topright", gr.name[-5], col=gr.col, bty="n", cex=1.1, lwd=3,
       y.intersp=0.75)


###################
# FR groups
###################

plot(NULL, xlim=c(0,8), ylim=c(0,0.5), main="", ylab="Functional diversity", xlab="Zero Flow days (ZFT)")

#mtext("h)", line = 1.0, adj = -0.2, cex = 1.6, font = 2)

# Plotting fitted values

a=1

for (j in 8:9) {
  
  mod_set  <- mod_set_list[[j]]
  
  mod.weights <- Weights(as.numeric(lapply(mod_set, AICc)))
  
  mod.pred <- data.frame(lapply(mod_set, function(x) { 
    
    predict(x,data.frame(ZFT=days.seq,
                         ZFP=mean(dat$ZFP),
                         ZFL=mean(dat$ZFL),
                         RE=mean(dat$RE),
                         do=mean(dat$do),
                         alt=mean(dat$alt)))
  }))
  
  # mean weighted prediction
  mean_pred <- (apply(mod.pred, 1, function(x) sum(x * mod.weights)))
  if (var(mean_pred)>0 & j==8) lines(per.seq, mean_pred, col=gr.col[1], lwd=3)
  if (var(mean_pred)>0 & j==9) lines(per.seq, mean_pred, col=gr.col[4], lwd=3)
  
}

legend("topright", gr.name[c(1,4)], col=gr.col[c(1,4)], bty="n", cex=1.1, lwd=3,
       y.intersp=0.75)


###################
# FR groups
###################

plot(NULL, xlim=c(0,8), ylim=c(0,0.5), main="", ylab="Functional diversity", xlab="Zero Flow Periods (ZFP)")

#mtext("h)", line = 1.0, adj = -0.2, cex = 1.6, font = 2)

# Plotting fitted values

for (j in c(8,10)) {
  
  mod_set  <- mod_set_list[[j]]
  
  mod.weights <- Weights(as.numeric(lapply(mod_set, AICc)))
  
  mod.pred <- data.frame(lapply(mod_set, function(x) { 
    
    predict(x,data.frame(ZFT=mean(dat$ZFT),
                         ZFP=sqrt(per.seq),
                         ZFL=mean(dat$ZFL),
                         RE=mean(dat$RE),
                         do=mean(dat$do),
                         alt=mean(dat$alt)))
  }))
  
  # mean weighted prediction
  mean_pred <- (apply(mod.pred, 1, function(x) sum(x * mod.weights)))
  if (var(mean_pred)>0 & j==8) lines(per.seq, mean_pred, col=gr.col[1], lwd=3)
  if (var(mean_pred)>0 & j==10) lines(per.seq, mean_pred, col=gr.col[5], lwd=3)
  
}

legend("topright", gr.name[c(1,5)], col=gr.col[c(1,5)], bty="n", cex=1.1, lwd=3,
       y.intersp=0.75)


dev.off()



