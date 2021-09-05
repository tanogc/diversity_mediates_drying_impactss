################################################################################################################
#                                                                                                               
# R script to calculate taxonomic and trait-based metrics for invertebrate communities                                                                    
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

# loading libraries
library(sqldf)
library(FD)
library(vegan)
library(ade4)
library(plyr)

# Supplementary functions
source("0_FD_functions.R")
source("0_quality_funct_space_fromdist.R")

# Loading taxonomic composition matrices
kick<-read.table("inv_kick.txt",h=T,sep="\t")
surber<-read.table("inv_surber.txt",h=T,sep="\t")

# Arranging matrices
row.names(surber)<-surber$site
row.names(kick)<-kick$site

surber<-surber[,-1]
kick<-kick[,-1]

# Dividing by surber surface (m2) to get ind / m2
surber<-surber/0.04

# Loading trait matrices
tr<-read.table("traits.txt",h=T,sep="\t")
tr_codes<-read.table("traits_description.txt",h=T,sep="\t")

# Arranging trait data
trophic<-tr[,64:ncol(tr)]
rownames(trophic)<-tr$genus
tr<-tr[,-c(64:ncol(tr))]

###################
# Response traits
###################

# arranging traits
rownames(tr)<-colnames(surber)
tr[,-1]->fuzzy_tr

fuzzy_tr<-fuzzy_tr[,-c(56:ncol(fuzzy_tr))] #
fuzzy_tr<-fuzzy_tr[,which(colSums(fuzzy_tr)!=0)] #u8 is removed

traits.blo<-c(7,2,3,4,8,4,5,5,7,9)
fuzzy_tr <- prep.fuzzy(fuzzy_tr, traits.blo)

# Combining the traits
tr.ktab<-ktab.list.df(list(fuzzy_tr))
tr.dist <- dist.ktab(tr.ktab, c("F")) # fuzzy-coding adapted Gower distance (pavoine et al 2009)

# Response traits: Estimating the optimum number of dimensions.
qual_fs_tr<-quality_funct_space_fromdist(tr.dist, nbdim=10)
qual_fs_tr$meanSD

# Building functional space (PCoA)
dudi.pco(tr.dist,scannf = F,nf=10)->tr.pco
length(which(ttr.pco$eig<0)) # checking for negative eigenvalues

# Total and cumulative variance explained by each axes
tr.pco$eig[1:8]/sum(tr.pco$eig)
sum(tr.pco$eig[1:5]/sum(tr.pco$eig)) # 5D
sum(tr.pco$eig[1:8]/sum(tr.pco$eig)) # 8D

# Correlations between original Gower distance and Euclidean distance
# in the functional space
par(mfrow=c(2,2))
for (i in 3:6) plot(dist(tr.pco$li[,1:i]), tr.dist, main=paste(i,"D", sep=""))

# Correlations between PCO axesandoriginal traits
round(cor(data.frame(fuzzy_tr),tr.pco$li,use="pairwise.complete.obs"),2)->cor.res
cor.res

# Taxon coordinates in the functional space
round(tr.pco$li[,c(1:3)],2)

write.table(cor.res,"cor_resp_traits.txt",sep="\t") # saving correlation table

###################
# trophic traits
###################

# preparing trophic trait data
traits.blo<-c(10) # blocks
trophic <- prep.fuzzy(trophic, traits.blo) 

# Combining the traits
tro.ktab<-ktab.list.df(list(trophic))
tro.dist <- dist.ktab(tro.ktab, "F") # fuzzy-coding adapted Gower distance (pavoine et al 2009)

# Response traits: Estimating the optimum number of dimensions.
qual_fs_tro<-quality_funct_space_fromdist(tro.dist, nbdim=10)
qual_fs_tro$meanSD

# Setting group colours
fgr.col<-c("#4600FFFF", "gold3", "dark green", "brown1", "darkorange3", "#CC00FFFF")

# Example of visualising taxa within a group
rownames(trophic) [which(trophic$GRA>=0.50)]
rownames(trophic) [which( (trophic$AFF+trophic$PFF )>=0.50)]

# Setting a threshold to classify taxa into functional groups
th<-0.50

# Creating an empty vector
f.grs<-rep(NA, nrow(trophic))

f.grs[which(trophic$GRA>=th)]<-"GRA"
f.grs[which(trophic$SHR>=th)]<-"SHR"
f.grs[which(trophic$AFF>=th | trophic$PFF>=th)]<-"FIL"
f.grs[which(trophic$PRED>=th)]<-"PRE"
f.grs[which(trophic$GAT>=th)]<-"GAT"
f.grs[which(is.na(f.grs))]<-"OMNI"

factor(f.grs)->f.grs

mean.fgr<-ddply(data.frame(f.grs, trophic), .(f.grs), function(x) round(colMeans(x[,-1]),2))
mean.fgr

# Trophic space
dudi.pco(tro.dist,scannf = F,nf=10)->t.pco

length(which(t.pco$eig<0)) # checking for negative eigenvalues

# Total and cummulative variance explained by each axes
t.pco$eig[1:3]/sum(t.pco$eig)
sum(t.pco$eig[1:3]/sum(t.pco$eig))

par(mfrow=c(2,2))
for (i in 3:6) plot(dist(tr.pco$li[,1:i]), r.dist, main=paste(i,"D", sep=""))

# Correlations 
round(cor(trophic,t.pco$li,use="pairwise.complete.obs"),2)->cor.res.t
cor.res.t

# Taxon coordinates in the functional space
round(t.pco$li[,c(1:4)],2)

write.table(cor.res.t,"pco_cor_trophic.txt",sep="\t") # saving correlation table

#### Taxonomic-based metrics
specnumber(kick)->ric
rowSums(surber)->abun

#### Functional trait-based metrics

# Functional redundancy (community and groups)
calc.FR(surber, f.grs, gr.names=unique(f.grs))->FR

colnames(FR$ab.fgrs)<-paste("FR", colnames(FR$ab.fgrs), sep=".")

# Functional diversity for within each trophic group
RD.gr<-rd.gr(log(surber+1), r=fuzzy_tr, traits.blo=c(7,2,3,4,8,4,5,5,7,9), 
             gr.names=unique(f.grs), f.grs, m=8)

colnames(RD.gr)<-paste("RD", colnames(RD.gr), sep=".")

# Saving generated metrics
res_int<-data.frame(ric=ric,
                    abun=abun,
                    FR=FR$FR.ab,
                    meanRD=rowMeans(scale(RD.gr)),
                    FR$ab.fgrs,
                    RD.gr)

# saving dataset
write.table(res_int,"final_set.txt",sep="\t")
