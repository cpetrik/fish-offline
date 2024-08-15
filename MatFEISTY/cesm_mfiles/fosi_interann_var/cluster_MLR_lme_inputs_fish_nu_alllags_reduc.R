# Cluster results of Mult linear regression of forcing on fish prod (nu)
# Includes all lags up to 2yrs
# Reduced to most sig of each driver and interaction
# From (8/13/2023) output

rm(list=ls())

library(rpart)
library(tree)
library(pvclust)
library(mclust)
library(dendextend)

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/corrs/"

### ------------------------------------------------------------
# load data
#datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/FOSI/"
#datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"
datar <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/"

fcoef <- read.csv(paste0(datar,"LMEs_mlr_nu_drivers_ALLdiv2SD_alllags_reduc_F.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef <- read.csv(paste0(datar,"LMEs_mlr_nu_drivers_ALLdiv2SD_alllags_reduc_P.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef <- read.csv(paste0(datar,"LMEs_mlr_nu_drivers_ALLdiv2SD_alllags_reduc_D.csv"),sep=",",header = T,stringsAsFactors = F)
acoef <- read.csv(paste0(datar,"LMEs_mlr_nu_drivers_ALLdiv2SD_alllags_reduc_A.csv"),sep=",",header = T,stringsAsFactors = F)

names(fcoef)[1] <- "LME"
names(pcoef)[1] <- "LME"
names(dcoef)[1] <- "LME"
names(acoef)[1] <- "LME"

### Use heatmap & clustering tree to find patterns

# Cluster LMEs
Acoef <- acoef[,2:9]
rownames(Acoef) <- acoef$LME
Acoef[is.na(Acoef)] <- 0
ad_fish <- dist(Acoef) # method="man" # is a bit better
ahc_fish <- hclust(ad_fish, method = "ward.D")
adfish <- as.dendrogram(ahc_fish)
# Color the branches based on the clusters:
adfish <- color_branches(adfish, k=6) #,k=8 
adend_6 <- cutree(adfish, k=6)
acoef$Cluster <- adend_6

## Forage
Fcoef <- fcoef[,2:9]
rownames(Fcoef) <- fcoef$LME
Fcoef[is.na(Fcoef)] <- 0
fd_fish <- dist(Fcoef) # method="man" # is a bit better
fhc_fish <- hclust(fd_fish, method = "ward.D")
fdfish <- as.dendrogram(fhc_fish)
# Color the branches based on the clusters:
fdfish <- color_branches(fdfish, k=5) #,k=8 
fdend_5 <- cutree(fdfish, k=5)
fcoef$Cluster <- fdend_5


## Lg Pel
Pcoef <- pcoef[,2:9]
rownames(Pcoef) <- pcoef$LME
Pcoef[is.na(Pcoef)] <- 0
pd_fish <- dist(Pcoef) # method="man" # is a bit better
phc_fish <- hclust(pd_fish, method = "ward.D")
pdfish <- as.dendrogram(phc_fish)
# Color the branches based on the clusters:
pdfish <- color_branches(pdfish, k=5) #,k=8 
pdend_5 <- cutree(pdfish, k=5)
pcoef$Cluster <- pdend_5


## Dem
Dcoef <- dcoef[,2:9]
rownames(Dcoef) <- dcoef$LME
Dcoef[is.na(Dcoef)] <- 0
dd_fish <- dist(Dcoef) # method="man" # is a bit better
dhc_fish <- hclust(dd_fish, method = "ward.D")
ddfish <- as.dendrogram(dhc_fish)
# Color the branches based on the clusters:
ddfish <- color_branches(ddfish, k=5) #,k=8 
ddend_5 <- cutree(ddfish, k=5)
dcoef$Cluster <- ddend_5




write.table(fcoef,paste0(datar,"LME_Fnu_mlr_coeffs_ALLdiv2SD_alllags_reduc_cluster.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datar,"LME_Pnu_mlr_coeffs_ALLdiv2SD_alllags_reduc_cluster.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datar,"LME_Dnu_mlr_coeffs_ALLdiv2SD_alllags_reduc_cluster.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datar,"LME_Anu_mlr_coeffs_ALLdiv2SD_alllags_reduc_cluster.csv"),sep=",",row.names=F)



