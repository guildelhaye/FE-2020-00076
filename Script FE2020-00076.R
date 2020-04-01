## R script related to the article Interspecific trait integration increases 
## with environmental harshness: a case study along a metal toxicity gradient 
## published in Functional Ecology by Guillaume Delhaye, David Bauman, Maxime 
## Séleck, Edouard Ilunga wa Ilunga, Grégory Mahy & Pierre Meerts.

library(ade4)
library(Hmisc)
library(FD)

######################################################################
## data preparation
traits.br <- read.table("Traits total.csv", header=T, row.names =1, sep =";") #traits data
cu <- as.data.frame(read.table("Cu.soil.csv", header=T, row.names =1, sep =";")) # soil data
cu.soil<-log10(cu$Cu)
abon <- as.data.frame(read.table("Abun.csv", header=T, row.names =1, sep =";")) # community data

traits.ef <- as.data.frame(na.omit(traits.br[,-c(1,6,8,10,11,12,15,16,17,18)]))# remove unnecessary traits

# create a data frame of transformed traits 
traits.f<-data.frame(matrix(nrow = nrow(traits.ef),ncol =ncol(traits.ef)))
traits.f[,1]=sqrt(traits.ef$HV)
traits.f[,2]=log(traits.ef$Prof)
traits.f[,3]=log(traits.ef$LA)
traits.f[,4]=sqrt(traits.ef$LDMC)
traits.f[,5]=log(traits.ef$SLA)
traits.f[,6]=log(traits.ef$SM)
traits.f[,7]=log(traits.ef$Cupl)
traits.f[,8]=log(traits.ef$Copl+0.1)
colnames(traits.f)<- c("VH", "RD", "LA", "LDMC", "SLA", "SM", "Cu", "Co")
rownames(traits.f)<-rownames(traits.ef)
traits.ef<-as.data.frame(scale(traits.f))
comm<-abon[,rownames(traits.ef)]
n_sites=nrow(cu)

#################################################################
# PCA on all traits
pc<-dudi.pca(traits.ef, scale = T, center = T, scannf = F, nf=4)
summary(pc)
(pc$eig/sum(pc$eig))*100

#################################################################
## Bivariate correlations between traits
rcorr(as.matrix(traits.ef), type = "spearman")

#################################################################
######## Multivariate trait coordination - FRic - FDis ##########
#################################################################
####
sp.rich<-NULL
FRic.obs<-NULL
FDis.obs<-NULL
FTIsd.obs<-NULL
FTIrange.obs<-NULL
FTIsd.obs.rank<-NULL
FTIrange.obs.rank<-NULL

ses.FRic<-NULL
ses.FDis<-NULL
ses.FTIsd<-NULL
ses.FTIrange<-NULL
ses.FTIsd.rank<-NULL
ses.FTIrange.rank<-NULL
nrep=1000

for (a in 1:n_sites){ 
  
  com.obs=comm[a,comm[a,]>0] # Keeping species with abundance >0
  traits.obs=na.omit(traits.ef[colnames(com.obs),]) #Keep species present 
  traits.obs.rank<-apply(traits.obs, 2, rank)
  com.obs=as.matrix(com.obs[,rownames(traits.obs)])
  sp.rich[a]<-nrow(traits.obs)
  
  ## Comptute observed FRic and FDis
  res.obs<-dbFD(traits.obs, com.obs, message=T)
  FRic.obs[a]<-res.obs$FRic
  FDis.obs[a]<-res.obs$FDis
  
  ## trait integration
  obs.pca <- dudi.pca(traits.obs, scannf=FALSE, center=TRUE, scale=TRUE)# rPCA on trait values
  FTIsd.obs[a]<-sd(obs.pca$eig/sum(obs.pca$eig))
  FTIrange.obs[a]<-max(obs.pca$eig/sum(obs.pca$eig))-min(obs.pca$eig/sum(obs.pca$eig))
  
  ## trait integration on rank data
  obs.pca.rank <- dudi.pca(traits.obs.rank, scannf=FALSE, center=TRUE, scale=TRUE)# PCA on rank of trait values
  FTIsd.obs.rank[a]<-sd(obs.pca.rank$eig/sum(obs.pca.rank$eig))
  FTIrange.obs.rank[a]<-max(obs.pca.rank$eig/sum(obs.pca.rank$eig))-min(obs.pca.rank$eig/sum(obs.pca.rank$eig))
  
  ############ NULL communities
  FRic.null<-NULL
  FTIsd.null<-NULL
  FTIrange.null<-NULL
  FTIsd.null.rank<-NULL
  FTIrange.null.rank<-NULL
  
  for(x in 1:nrep){
    
    # Local null community sampled from entier dataset
    com.null.FRic<-sort(sample(comm[a,], sp.rich[a], replace=F ))
    traits.null.FRic<-traits.ef[colnames(com.null.FRic),]
    traits.null.rank <-apply(traits.null.FRic, 2, rank)
    com.null.FRic<- com.null.FRic+1
    #Compute FRic null
    FRic.null[x]<-dbFD(traits.null.FRic, com.null.FRic, message=F)$FRic
    #Calc FTI on the same trait matrix as FRic
    null.pca<-dudi.pca(traits.null.FRic, scannf=FALSE, center=TRUE, scale=TRUE)
    FTIsd.null[x]<-sd(null.pca$eig/sum(null.pca$eig))
    FTIrange.null[x]<-max(null.pca$eig/sum(null.pca$eig))-min(null.pca$eig/sum(null.pca$eig))
    #Calc FTI on rank data
    null.pca.rank<-dudi.pca(traits.null.rank, scannf=FALSE, center=TRUE, scale=TRUE)
    FTIsd.null.rank[x]<-sd(null.pca.rank$eig/sum(null.pca.rank$eig))
    FTIrange.null.rank[x]<-max(null.pca.rank$eig/sum(null.pca.rank$eig))-min(null.pca.rank$eig/sum(null.pca.rank$eig))
  }
  
  FDis.null<-NULL
  #Local null community sampled from local dataset, shuffling abundances
  for (y in 1:nrep){
    com.null.FDis=sample(com.obs[1,])
    names(com.null.FDis)<-colnames(com.obs)
    traits.null.FDis<-traits.obs
    # Compute FDis null
    FDis.null[y]<-dbFD(traits.null.FDis, com.null.FDis, messages = F)$FDis
  }
  
  ses.FRic[a]<-(FRic.obs[a]-mean(FRic.null))/sd(FRic.null)
  ses.FDis[a]<-(FDis.obs[a]-mean(FDis.null))/sd(FDis.null)
  ses.FTIsd[a]<-(FTIsd.obs[a]-mean(FTIsd.null))/sd(FTIsd.null)
  ses.FTIrange[a]<-(FTIrange.obs[a]-mean(FTIrange.null))/sd(FTIrange.null)
  ses.FTIsd.rank[a]<-(FTIsd.obs.rank[a]-mean(FTIsd.null.rank))/sd(FTIsd.null.rank)
  ses.FTIrange.rank[a]<-(FTIrange.obs.rank[a]-mean(FTIrange.null.rank))/sd(FTIrange.null.rank)
}
cor.FRic<-cor.test(cu.soil, ses.FRic, method = "spearman")
cor.FDis<-cor.test(cu.soil, ses.FDis, method = "spearman")
cor.FTIsd<-cor.test(cu.soil, ses.FTIsd, method = "spearman")
cor.FTIrange<-cor.test(cu.soil, ses.FTIrange, method = "spearman")
cor.FTIsd.rank<-cor.test(cu.soil, ses.FTIsd.rank, method = "spearman")
cor.FTIrange.rank<-cor.test(cu.soil, ses.FTIrange.rank, method = "spearman")

(cor.FRic.sp<-cor.test(sp.rich,ses.FRic, method = "spearman"))
(cor.FDis.sp<-cor.test(sp.rich,ses.FDis, method = "spearman"))
(cor.cov.rich<-cor.test(sp.rich,ses.FTIrange.rank, method = "spearman"))
(cor.cov.rich.sd<-cor.test(sp.rich,ses.FTIsd.rank, method = "spearman"))

###############################################################
#### Figure for indices comparison
par(mfrow=c(4,2))

#ITIrange - Cu
plot(cu.soil,ses.FTIrange.rank,pch=20, cex = 1.5, bty="l", xlim = c(2,4), 
     ylim= c(-2.5,2.5), tck=0.02,
     ylab = "sesITI range", xlab="Soil Cu content (log)",
     main= c("Rho=",round(cor.FTIrange.rank$estimate, 2)))
lines(cu.soil,ses.FTIrange.rank,type="h")
abline(h=0, col = "black", lwd = 2, lty =1)

#ITIrange - sp.rich
res<-as.data.frame(cbind(ses.FTIrange.rank,sp.rich))
res<-res[order(res$ses.FTIrange.rank),]
plot(res$ses.FTIrange.rank,res$sp.rich, pch=20, cex=1.5, bty = "l", 
     ylim=c(15,35), xlim=c(-2,2), tck=0.02,
     ylab = "Species richness", xlab="ITI range", 
     main= c("Rho=",round(cor.cov.rich$estimate, 2)))
ktr<-loess(res$sp.rich~res$ses.FTIrange.rank, span = 1.3)
lines(res$ses.FTIrange.rank, ktr$fitted, lwd=2)

#FRic - Cu
plot(cu.soil,ses.FRic,pch=20, cex = 1.5, bty="l",xlim = c(2,4),ylim= c(-2.5,2.5), 
     tck=0.02, ylab = "sesFRic", xlab="Soil Cu content (log)",
     main= c("Rho=",round(cor.FRic$estimate, 2)))
lines(cu.soil,ses.FRic,type="h")
abline(h=0, col = "black", lwd = 2, lty =1)

#FRic - sp.rich
res<-as.data.frame(cbind(ses.FRic,sp.rich))
res<-res[order(res$ses.FRic),]
plot(res$ses.FRic,res$sp.rich, pch=20, cex=1.5, bty = "l", ylim=c(15,35),tck=0.02,
     ylab = "Species richness", xlab="sesFRic", 
     main= c("Rho=",round(cor.FRic.sp$estimate, 2)))

#FDis - Cu
plot(cu.soil,ses.FDis, pch=20, cex = 1.5, bty="l", xlim = c(2,4), ylim= c(-2.5,2.5), 
     tck=0.02, ylab = "sesFDis", xlab="Soil Cu content (log)",
     main= c("Rho=",round(cor.FDis$estimate, 2)))
lines(cu.soil,ses.FDis,type="h")
abline(h=0, col = "black", lwd = 2, lty =1)

#FDis - sp.rich
res<-as.data.frame(cbind(ses.FDis,sp.rich))
res<-res[order(res$ses.FDis),]
plot(res$ses.FDis,res$sp.rich, pch=20, cex=1.5, bty = "l", ylim=c(15,35), 
     tck=0.02, ylab = "Species richness", xlab="sesFDis", 
     main= c("Rho=",round(cor.FDis.sp$estimate, 2)))
ktr<-loess(res$sp.rich~res$ses.FDis, span = 1.3)
lines(res$ses.FDis, ktr$fitted, lwd=2)

#ITIsd - Cu
plot(cu.soil,ses.FTIsd.rank,pch=20, cex = 1.5, bty="l", xlim = c(2,4), ylim= c(-2,2),
     ylab = "sesITI sd", xlab="Soil Cu content (log)",
     main= c("Rho=",round(cor.FTIsd.rank$estimate, 2)))
lines(cu.soil,ses.FTIsd.rank,type="h")
abline(h=0, col = "black", lwd = 2, lty =1)

#ITIsd - sp.rich
res<-as.data.frame(cbind(ses.FTIsd.rank,sp.rich))
res<-res[order(res$ses.FTIsd.rank),]
plot(res$ses.FTIsd.rank,res$sp.rich, pch=20, cex=1.5, bty = "l", ylim=c(15,35),
     ylab = "Species richness", xlab="sesITI sd", 
     main= c("Rho=",round(cor.cov.rich.sd$estimate, 2)))
ktr<-loess(res$sp.rich~res$ses.FTIsd.rank, span = 1.3)
lines(res$ses.FTIsd.rank, ktr$fitted, lwd=2)

#################################################################
####################################################################
###### BIVARIATE ANALYSES #####################################
###############################################################

## Bivariate analysis of trait covariation based on spearman correlations
## => correlations on the rank of the trait values to investigate non-linear relationships. 
par(mfrow=c(7,7), mar=c(1,1,1,1))

result<-data.frame(matrix(nrow = ncol(traits.ef),ncol =ncol(traits.ef)))
colnames(result)<-colnames(traits.ef)
rownames(result)<- colnames(traits.ef)
null.aver<-vector()
obs.cor <- NULL
x=NULL
y=NULL
a=NULL
com.obs=NULL
min.sp <- min(sp.rich)
nreps = 1000

for (x in 2:ncol(traits.ef)){ ### Trait 1
  
  for(y in 1:(ncol(traits.ef)-1)){ ### Trait 2
    
    obs.cor<-NULL
    obs.cor.p<-NULL
    for (a in 1:n_sites){  ### On all sites
      
      com.obs=comm[a,comm[a,]>0] # Keeping species with abundance >0
      traits.obs=na.omit(traits.ef[colnames(com.obs),]) #Keep traits of species present in communities
      com.obs=com.obs[,rownames(traits.obs)] # Keep species in communities for which we have traits
      proba<-as.vector(com.obs/sum(com.obs))
      
      ## Standardized each plot by the smallest taxonomic richness
      cor.loc<-NULL 
      
      for(n.sample in 1:nreps){
        traits.sub<-traits.obs[sample(nrow(traits.obs), min.sp, prob=proba),]
        tr1 <- traits.sub[,x]
        tr2 <- traits.sub[,y]
        cor.loc[n.sample] <- rcorr(tr1,tr2, type = "spearman")$r[1,2]
        }
      
      obs.cor[a]<-mean(cor.loc)
    }
    
    if(mean(obs.cor)!=1){
      
      result[x,y] <-round(cor.test(cu.soil,obs.cor, method="spearman")$estimate,2)
      (p<- cor.test(cu.soil,obs.cor, method="spearman", exact = F)$p.value)
      
      plot(cu.soil,obs.cor, pch = "", tck=0.02, bg=1, xlab = "", ylab="", ylim=c(-0.9,0.9), xlim=c(2,4), bty="l")
      abline(h=0, lty=2)
      
      for(pr in 1:n_sites){
        points(cu.soil[pr], obs.cor[pr], col="grey", pch = 19)
      }
      if(p<0.001){
        hg <- lm(obs.cor~cu.soil)
        lines(cu.soil,hg$fitted,col="black",lwd=1.5)
        abline(h=0, lty=2)
      }
    }else{ 
      plot(cu.soil,obs.cor, pch = "", tck = 0.02, bg=1, xlab = "", ylab="", ylim=c(-0.9,0.9), xlim=c(4.5,9), bty="l")
    }
  }
}
