#################################################################################################################
#                                                                                                               
# Functions to estimate Functional Diversity metrics                                                         
#                                                                                                               
# Code written by Cayetano Gutiérrez-Cánovas                                               
# email for queries: cayeguti@um.es                                                                             
#                                                                                                               
#################################################################################################################

# prep.fuzzy.df() Functions to prepare matrices, construct the functional space and to estimate functional diversity components
# 
# Inputs:
# traits: fuzzy coding traits
# col.blocks: blocks indicating the number of trait categories per trait
#
# Outputs
# traits: fuzzy coding traits as percentages

prep.fuzzy.df<-function (traits, col.blocks) 
{
  if (!is.data.frame(traits)) 
    stop("Data.frame expected")
  if (sum(col.blocks) != ncol(traits)) {
    stop("Non convenient data in col.blocks")
  }
  if (is.null(names(col.blocks))) {
    names(col.blocks) <- paste("FV", as.character(1:length(col.blocks)), sep = "")
  }
  f1 <- function(x) {
    a <- sum(x)
    if (is.na(a)) 
      return(rep(0, length(x)))
    if (a == 0) 
      return(rep(0, length(x)))
    return(x/a)
  }
  k2 <- 0
  col.w <- rep(1, ncol(traits))
  
  for (k in 1:(length(col.blocks))) {
    k1 <- k2 + 1
    if (col.blocks[k]==1) k2<-k1 else k2 <- k2 + col.blocks[k]
    X <- as.matrix(traits[, k1:k2])
    if (col.blocks[k]==1) X[which(X[,1]>0),]<-1 else X <- t(apply(X, 1, f1))
    X.marge <- apply(X, 1, sum)
    X.marge <- X.marge
    X.marge <- X.marge/sum(X.marge)
    X.mean <- apply(X * X.marge, 2, sum)
    nr <- sum(X.marge == 0)
    cat(nr, "missing data found in block", k, "\n")
    traits[, k1:k2] <- X
    col.w[k1:k2] <- X.mean
  }
  attr(traits, "col.blocks") <- col.blocks
  attr(traits, "col.freq") <- col.w
  col.num <- factor(rep((1:length(col.blocks)), col.blocks))
  attr(traits, "col.num") <- col.num
  return(traits)
}

# chull_3d() estimates the convex hull of a Functional Space
#
# Inputs:
# fpc: functional space
# m: number of axes to select
# prec: convex hull precision ("Qt" or "QJ")
#
# Output:
# a vector with the Functional Richness of each community

chull_3d<-function(fpc,m,prec=c("Qt","QJ")){
  
  convhulln(fpc[,1:m], c("FA",prec))$vol->fric.3d.max
  
  return(fric.3d.max)
  }
  
# fric_3d() estimates the Functional Richness of a set of communties
# This function computes the hypervolume to estimate how each community fills
# the functional space
#
# Inputs:
# taxa: community data
# fpc: functional space
# m: number of axes to select
# prec: convex hull precision ("Qt" or "QJ")
# fric.3d.max: volume of the regional convex hull (set fric.3d.max=1 to have the absolute FRic value)
#
# Output:
# a vector with the Functional Richness of each community

fric_3d<-function(taxa,fpc,m,prec=c("Qt","QJ"), fric.3d.max=NULL){
  fric.3d<-rep(NA,nrow(taxa))
  names(fric.3d)<-rownames(taxa)
  
  # Convex hull of regional pool
  if(is.null(fric.3d.max)==T) convhulln(fpc[,1:m], c("FA",prec))$vol->fric.3d.max
  
  specnumber(taxa)->ric
  
  for (com in 1:nrow(taxa)){
    fpc[which(unlist(rep(taxa[com,]))>0),1:m]->tr.com
    if (ric[com]>=m+1) convhulln(tr.com, c("FA",prec))$vol/fric.3d.max->fric.3d[com] else NA->fric.3d[com]
  }
  return(fric.3d)
}

# fric_1d() estimates the Functional Richness of a set of communties
# This function computes the hypervolume to estimate how each community fills
# the functional space
#
# Inputs:
# taxa: community data
# fpc: functional space
# m: number of axes to select
#
# Output:
# a vector with the Functional Richness of each community

fric_1d<-function(taxa,fpc,m){
  
  fric.1d<-rep(NA,nrow(taxa))
  names(fric.1d)<-rownames(taxa)
  fric.1d.max<-rep(NA,m)
  
  for (i in 1:m) sum(abs(range(fpc[,i])))->fric.1d.max[i]
  
  specnumber(taxa)->ric
  for (com in 1:nrow(taxa)){
    fpc[which(unlist(rep(taxa[com,]))>0),1:m]->tr.com
    
    if (ric[com]>=1) mean(sum(abs(range(tr.com)))/fric.1d.max)->fric.1d[com] else NA->fric.1d[com]
  }
  return(fric.1d)
}

# fdisp_k() estimates the Functional Dispersion of a set of communties
# This function computes the weighted mean distance to community centroid in
# the functional space
# 
# Function modified from Lalibert� & Legendre (2010) Ecology
#
# Inputs:
# d: trait dissimilarity matrix
# a: community data
# m: number of axes to select
# tol: tolerance threshold to test whether the distance matrix is Euclidean
#
# Output:
# FDis: a vector with the Functional Dispersion of each community
# eig: eigenvectors of each functional axis
# vectors: functional axes

fdisp_k<-function (d, a, m, tol = 1e-07) 
{
  if (!inherits(d, "dist")) 
    stop("'d' must be a 'dist' object.")
  n <- attr(d, "Size")
  if (is.null(attr(d, "Labels"))) 
    stop("'d' must have labels.", "\n")
  else sn.d <- attr(d, "Labels")
  if (missing(a)) {
    ab.names <- list("Community1", sn.d)
    a <- matrix(1, 1, n, dimnames = ab.names)
  }
  com <- nrow(a)
  if (ncol(a) != n) 
    stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
  if (is.null(colnames(a))) 
    stop("'a' must have column names", "\n")
  else sn.a <- colnames(a)
  if (any(sn.d != sn.a)) 
    warning("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
            "\n")
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  if (any(abun.sum == 0)) 
    warning("At least one community has zero-sum abundances (no species).", 
            "\n")
  abun.sum2 <- apply(a, 2, sum)
  if (any(abun.sum2 == 0)) 
    warning("At least one species does not occur in any community (zero total abundance across all communities).", 
            "\n")
  if (any(is.na(d))) 
    stop("NA's in the distance matrix.", "\n")
  A <- matrix(0, ncol = n, nrow = n)
  A[row(A) > col(A)] <- -0.5 * d^2
  A <- A + t(A)
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  w0 <- eig[n]/eig[1]
  if (w0 > -tol) 
    r <- sum(eig > (eig[1] * tol))
  else r <- length(eig)
  vectors <- vectors[, 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
                                                   r)
  dimnames(vectors) <- list(colnames(a), NULL)
  pos <- eig > 0
  if (m>length(pos)) m <- length(pos)
  if (m>0) pos<- c(pos[1:m],rep(F,length(pos)-m))
  
  avg.dist.cent <- rep(NA, nrow(a))
  names(avg.dist.cent) <- row.names(a)
  for (i in 1:com) {
    pres <- which(a[i, ] > 0)
    nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
    if (nb.sp >= 2) {
      w <- a[i, pres]
      centroid <- apply(vec, 2, weighted.mean, w = w)
      dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
      dist.pos <- rowSums(dist.pos^2)
      if (any(!pos)) {
        dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
        dist.neg <- rowSums(dist.neg^2)
      }
      else dist.neg <- 0
      zij <- sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] <- weighted.mean(zij, w)
    }
    else avg.dist.cent[i] <- 0
  }
  return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
}

# fdisp_k_sub() estimates the Functional Dispersion of a set of communties
#
# This function computes the weighted mean distance to community centroid in
# the functional space for a subset of taxa (taxonomic matric has less taxa than
# trait dissimilarity matrix)
# 
# Function modified from LalibertÃ© & Legendre (2010) Ecology
#
# Inputs:
# d: trait dissimilarity matrix
# a: community data
# m: number of axes to select
# tol: tolerance threshold to test whether the distance matrix is Euclidean
# tax_sub: taxonomic subset for which FDis will be computed
#
# Output:
# FDis: a vector with the Functional Dispersion of each community
# eig: eigenvectors of each functional axis
# vectors: functional axes

fdisp_k_sub<-function (d, a, tax_sub=NULL, m, tol = 1e-07) 
{
  if (!inherits(d, "dist")) 
    stop("'d' must be a 'dist' object.")
  n <- attr(d, "Size")
  if (is.null(attr(d, "Labels"))) 
    stop("'d' must have labels.", "\n")
  else sn.d <- attr(d, "Labels")
  if (missing(a)) {
    ab.names <- list("Community1", sn.d)
    a <- matrix(1, 1, n, dimnames = ab.names)
  }
  com <- nrow(a)
  #if (ncol(a) != (n-length(tax_sub))) 
  #  stop("Number of columns in 'a' must be equal to the number of objects in 'd'.")
  if (is.null(colnames(a))) 
    stop("'a' must have column names", "\n")
  else sn.a <- colnames(a)
  #if (any(sn.d != sn.a)) 
  #  stop("Species labels in 'd' and 'a' need to be identical and ordered alphabetically (or simply in the same order).", 
  #       "\n")
  a[which(is.na(a))] <- 0
  abun.sum <- apply(a, 1, sum)
  if (any(abun.sum == 0)) 
    warning("At least one community has zero-sum abundances (no species).", 
            "\n")
  abun.sum2 <- apply(a, 2, sum)
  if (any(abun.sum2 == 0)) 
    warning("At least one species does not occur in any community (zero total abundance across all communities).", 
            "\n")
  if (any(is.na(d))) 
    stop("NA's in the distance matrix.", "\n")
  A <- matrix(0, ncol = n, nrow = n)
  A[row(A) > col(A)] <- -0.5 * d^2
  A <- A + t(A)
  G <- bicenter.wt(A)
  e <- eigen(G, symmetric = TRUE)
  vectors <- e$vectors
  eig <- e$values
  w0 <- eig[n]/eig[1]
  if (w0 > -tol) 
    r <- sum(eig > (eig[1] * tol))
  else r <- length(eig)
  rownames(vectors) <- attr(d, "Labels")
  vectors <- vectors[(intersect(rownames(vectors), tax_sub)), 1:r, drop = FALSE] %*% diag(sqrt(abs(eig <- eig[1:r])), 
                                                   r)
  pos <- eig > 0
  if (m>0) pos<-c(pos[1:m],rep(F,length(pos)-m))
  
  avg.dist.cent <- rep(NA, nrow(a))
  names(avg.dist.cent) <- row.names(a)
  for (i in 1:com) {
    pres <- which(a[i, ] > 0)
    nb.sp <- nrow((unique(vec <- vectors[pres, , drop = F])))
    if (nb.sp >= 2) {
      w <- a[i, pres]
      centroid <- apply(vec, 2, weighted.mean, w = w)
      dist.pos <- sweep(vec[, pos, drop = F], 2, centroid[pos])
      dist.pos <- rowSums(dist.pos^2)
      if (any(!pos)) {
        dist.neg <- sweep(vec[, !pos, drop = F], 2, centroid[!pos])
        dist.neg <- rowSums(dist.neg^2)
      }
      else dist.neg <- 0
      zij <- sqrt(abs(dist.pos - dist.neg))
      avg.dist.cent[i] <- weighted.mean(zij, w)
    }
    else avg.dist.cent[i] <- 0
  }
  return(list(FDis = avg.dist.cent, eig = eig, vectors = vectors))
}

# feve_k() estimates the Functional Evenness of a set of communties
# This function computes regularity in the distribution of taxa or abundances
# across the Functional Space, using the Spanning Tree Method
# 
# Function modified from LalibertÃ© & Legendre (2010) Ecology
#
# Inputs:
# fpc: functional space axes
# taxa: community data
# m: number of axes to select
#
# Output:
# FEve: a vector with the Functional Evenness of each community

feve_k<-function(fpc,taxa,m){
  
  rdt=1
  
  # Creating variables
  nrow(taxa)->c
  FEve <- rep(NA, c) ; names(FEve) <- row.names(taxa)
  nb.sp<-specnumber(taxa)
  
  # generating the taxonomic matrix arranged according to the replicated trait values
  tax.pool<-ncol(taxa)
  taxa.rep<-data.frame(matrix(NA,nrow(taxa),tax.pool*rdt))
  spp.list<-c(1:(tax.pool*rdt))
  
  for (spp in 1:tax.pool) {paste(rep("spp",rdt),spp,sep="")->spp.list[((spp-1)*rdt+1):(spp*rdt)]}
  
  colnames(taxa.rep)<-spp.list
  
  for (spp in 1:tax.pool){taxa.rep[,((spp-1)*rdt+1):(spp*rdt)]<-taxa[,spp]/rdt}                             
  
  # Estimating Functional Evenness for each community
  
  for (i in 1:c) {
    sppres <- which(taxa.rep[i, ] > 0)
    # number of species in the community
    S <- length(sppres)
    ab <- as.matrix(taxa.rep[i, sppres])
    # scaling of abundances
    abundrel <- ab / sum(ab)
    
    # selecting the c
    tr <- data.frame(fpc[sppres,1:m ])
    
    if (nb.sp[i] > 2) {
      tr.dist <- dist(tr)
      linkmst <- mst(tr.dist)
      mstvect <- as.dist(linkmst)
      abund2 <- matrix(0, nrow = S, ncol = S)
      for (q in 1:S) for (r in 1:S) abund2[q, r] <- abundrel[q] +  abundrel[r]
      abund2vect <- as.dist(abund2)
      EW <- rep(0, S - 1)
      flag <- 1
      for (mv in 1:((S - 1) * S/2)) {
        if (mstvect[mv] != 0) {
          EW[flag] <- tr.dist[mv]/(abund2vect[mv])
          flag <- flag + 1
        }
      }
      minPEW <- rep(0, S - 1)
      OdSmO <- 1/(S - 1)
      for (l in 1:(S - 1)) minPEW[l] <- min((EW[l]/sum(EW)), 
                                            OdSmO)
      FEve[i] <- ((sum(minPEW)) - OdSmO)/(1 - OdSmO)
    } else FEve[i] <- NA
  }
  return(FEve)
}

# Function to plot 2D functional spaces
#
# coord: coordinates of two FS axes
# col_ch: colour of the FS filling
# col_points: colour of the FS vertices (points)

plot_chull2D<-function(coord,col_ch="#458B0050",border_ch="#458B00"){
  vert0<-convhulln( coord ,"Fx TO 'vert.txt'")
  vert1<-scan("vert.txt",quiet=T)
  vert_ij<-(vert1+1)[-1]
  points(coord[vert_ij,],pch=15, cex=1.5,col=border_ch)
  polygon(coord[vert_ij,], col=col_ch, border=border_ch)
}

# calc.FR() estimates the Functional Redundancy (FR) of a set of communties
# This function estimates the taxonomic richness for each taxonomic group,
# estating FR as the ratio between species richness and the number of functional groups
# in a given community.
#
# Inputs:
# taxa: community data
# groups: a grouping vector with the Funtional Groups for each taxon
#
# Outputs:
#
# $nbsp: taxonomic richness
# $ric.fgrs: the taxonomic richness for each Functional Group
# $sd.fgr: the among-Functional-group SD for taxonomic richness for each community
# $FGR: the number of Functional Gruoups for each community
# $FR: the Functional Richness for each community

calc.FR<-function(taxa,groups,gr.names=unique(groups)){
  list()->res
  
  if (ncol(taxa)!=length(groups)) stop("Trait and taxonomic matrix have different number of taxa")
  specnumber(taxa)->res$nbsp
  
  res$ab.fgrs<-res$ric.fgrs<-data.frame(matrix(NA,nrow(taxa),length(unique(groups))))
  colnames(res$ab.fgrs)<-colnames(res$ric.fgrs)<-gr.names
  
  j<-0
  
  for (i in gr.names) {
    j<-j+1
    if(is.vector(taxa[,which(groups==i)])==T) decostand(taxa[,which(groups==i)],"pa")->res$ric.fgrs[,j] else specnumber(taxa[,which(groups==i)])->res$ric.fgrs[,j]
  }
  
  j<-0
  
  for (i in gr.names) {
    j<-j+1
    if(is.vector(taxa[,which(groups==i)])==T) taxa[,which(groups==i)]->res$ab.fgrs[,j] else rowSums(taxa[,which(groups==i)])->res$ab.fgrs[,j]
  }
  
  apply(res$ric.fgrs,1,sd)->res$sd.fgr
  specnumber(res$ric.fgrs)->res$FGR # functional group richness estimation
  res$FR<-res$nbsp/res$FGR
  res$FR.ab<-rowSums(taxa)/res$FGR
  return(res)
}

# Functional diversity within functional group

# Inputs:
# taxa: community data
# r: response traits
# traits.blo: traits blocks
# gr.names: functional group names
# e.gr: a grouping vector with the Functional Groups for each taxon
# m: number of functional space dimensions to be considered
#
# Outputs:
#
# RD.g: functional diversity within each functional group for each site

rd.gr<-function(taxa, r, traits.blo, gr.names=NULL, e.gr, m){
  
  length(unique(e.gr))->cut.g
  
  if(is.null(gr.names==T)) names(e.gr)->gr.names
  
  # Response diversity per group
  RD.g<-data.frame(matrix(NA,nrow(taxa),cut.g)) # Keep the total number of groups
  colnames(RD.g)<-gr.names 
  rownames(taxa)->rownames(RD.g)
  
  for (i in unique(e.gr)){ 
    trg<-r[which(e.gr==i),]
    fuzzy_tr <- prep.fuzzy(trg, traits.blo)
    
    # Combining the traits
    tr.ktab<-ktab.list.df(list(fuzzy_tr))
    tr.dist <- dist.ktab(tr.ktab, c("F")) # fuzzy-coding adapted Gower distance (pavoine et al 2009)
    
    as.matrix(taxa[,which(e.gr==i)])->inv
    fdisp_k(tr.dist,inv,m)$FDis->RD.g[,i] 
  }
  
  return(RD.g)
}

# Function to obtain average variance partitioning from a set of models

hier.r2.res<-function(mod.set, sel.pred) {
  
  hier.r2.mod<-data.frame(matrix(0, length(mod.set), length(sel.pred)))
  names(hier.r2.mod)<-sel.pred
  rownames(hier.r2.mod)<-names(mod.set)
  
  for (j in 1:length(mod.set)){
    calcVarPart(mod.set[[j]])->rebeca
    
    rebeca[-length(rebeca)]->rebeca
    rebeca[order(names(rebeca))]->hier.r2.mod[j, which(names(hier.r2.mod) %in% names(rebeca))]
  }
  
  return(hier.r2.mod)
}  
