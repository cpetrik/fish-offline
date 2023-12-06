# Heatmaps of Mult linear regression of forcing on FishMIP catches
# Includes 2 lags: 0-1 or 1-2
# From (9/19/2023) output
# No zmeso
# Reduce to LMEs with high R2 value
# 2-4 clusters

rm(list=ls())

library(ggplot2)
library(cowplot) #plot_grid
library(gridExtra)
library(ggpubr) #cowplot and gridExtra
library(gplots)
library(reshape2) #melt

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
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/FOSI/"
datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"
#datar <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/"
dataf <- "/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/"

fcoef <- read.csv(paste0(datar,"LME_FishMIP_catch_F_mlr_coeffs_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef <- read.csv(paste0(datar,"LME_FishMIP_catch_P_mlr_coeffs_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef <- read.csv(paste0(datar,"LME_FishMIP_catch_D_mlr_coeffs_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",header = T,stringsAsFactors = F)
acoef <- read.csv(paste0(datar,"LME_FishMIP_catch_A_mlr_coeffs_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",header = T,stringsAsFactors = F)

lid <- c(1:22,24:32,34:61,63:66)
fcoef$LME <-lid
pcoef$LME <-lid
dcoef$LME <-lid
acoef$LME <-lid

# R2>0.5
fcoef25 <- subset(fcoef, R.2 >= 0.25)
pcoef25 <- subset(pcoef, R.2 >= 0.25)
dcoef25 <- subset(dcoef, R.2 >= 0.25)
acoef25 <- subset(acoef, R.2 >= 0.25)

### Use heatmap & clustering tree to find patterns

# Cluster LMEs
Acoef <- acoef25[,2:13]
rownames(Acoef) <- acoef25$LME
# anm <- c("Det1","Det2","TB1","TB2","TP1","TP2","ZmLoss1","ZmLoss2")
# Acoef <- Acoef[,anm]
Acoef[is.na(Acoef)] <- 0
ad_fish <- dist(Acoef) # method="man" # is a bit better
ahc_fish <- hclust(ad_fish, method = "ward.D")
adfish <- as.dendrogram(ahc_fish)
# Color the branches based on the clusters:
adfish <- color_branches(adfish, k=4) 
# reduce the size of the labels:
adfish <- set(adfish, "labels_cex", 0.8)

some_col_func <- colorspace::diverge_hcl(25)
col_breaks = c(-15,-5,-2,-1.1, # for low
               seq(-1,1,length=21), # for blue
               1.1,2,5,15) # for high

png(paste0(figp,"hclust_wardD_heatmap_LME_catch_A_mlr_coeffs_R2_ALLdiv2SD_alllags_4clus.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Acoef), 
                  main = "All fish",
                  dendrogram = "row",
                  Rowv = adfish,
                  #RowSideColors = ccol$col,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(6,6),      
                  srtCol = 45, #45 deg angle labels
                  offsetCol=-0.6,
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  keysize = 1.5,
                  density.info = "none",
                  col = some_col_func #, breaks=col_breaks
)
dev.off()


## Forage
Fcoef <- fcoef25[,2:13]
rownames(Fcoef) <- fcoef25$LME
Fcoef[is.na(Fcoef)] <- 0
fd_fish <- dist(Fcoef) # method="man" # is a bit better
fhc_fish <- hclust(fd_fish, method = "ward.D")
fdfish <- as.dendrogram(fhc_fish)
# Color the branches based on the clusters:
fdfish <- color_branches(fdfish,k=4) 
# reduce the size of the labels:
fdfish <- set(fdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_catch_F_mlr_coeffs_R2_ALLdiv2SD_alllags_4clus.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Fcoef), 
                  main = "Forage fish",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = fdfish,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(6,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()


## Lg Pel
Pcoef <- pcoef25[,2:13]
rownames(Pcoef) <- pcoef25$LME
Pcoef[is.na(Pcoef)] <- 0
pd_fish <- dist(Pcoef) # method="man" # is a bit better
phc_fish <- hclust(pd_fish, method = "ward.D")
pdfish <- as.dendrogram(phc_fish)
# Color the branches based on the clusters:
pdfish <- color_branches(pdfish,k=4) #
# reduce the size of the labels:
pdfish <- set(pdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_catch_P_mlr_coeffs_R2_ALLdiv2SD_alllags_4clus.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Pcoef), 
                  main = "Large Pelagics",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = pdfish,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(6,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()


## Dem
Dcoef <- dcoef25[,2:13]
rownames(Dcoef) <- dcoef25$LME
Dcoef[is.na(Dcoef)] <- 0
dd_fish <- dist(Dcoef) # method="man" # is a bit better
dhc_fish <- hclust(dd_fish, method = "ward.D")
ddfish <- as.dendrogram(dhc_fish)
# Color the branches based on the clusters:
ddfish <- color_branches(ddfish,k=4) #
# reduce the size of the labels:
ddfish <- set(ddfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_catch_D_mlr_coeffs_R2_ALLdiv2SD_alllags_4clus.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Dcoef), 
                  main = "Demersals",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = ddfish,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(6,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()

## All
adend_k <- cutree(adfish, k=4)
acoef25$Cluster <- adend_k

## Forage
fdend_k <- cutree(fdfish, k=4)
fcoef25$Cluster <- fdend_k

## Lg Pel
pdend_k <- cutree(pdfish, k=4)
pcoef25$Cluster <- pdend_k

## Dem
ddend_k <- cutree(ddfish, k=4)
dcoef25$Cluster <- ddend_k


acoef25 <- acoef25[,c("LME","Cluster")]
acoef <- merge(acoef,acoef25,by="LME",all=T)

fcoef25 <- fcoef25[,c("LME","Cluster")]
fcoef <- merge(fcoef,fcoef25,by="LME",all=T)

pcoef25 <- pcoef25[,c("LME","Cluster")]
pcoef <- merge(pcoef,pcoef25,by="LME",all=T)

dcoef25 <- dcoef25[,c("LME","Cluster")]
dcoef <- merge(dcoef,dcoef25,by="LME",all=T)


write.table(fcoef,paste0(datar,"LME_catch_F_mlr_coeffs_ALLdiv2SD_alllags_R2_025_4cluster.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datar,"LME_catch_P_mlr_coeffs_ALLdiv2SD_alllags_R2_025_4cluster.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datar,"LME_catch_D_mlr_coeffs_ALLdiv2SD_alllags_R2_025_4cluster.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datar,"LME_catch_A_mlr_coeffs_ALLdiv2SD_alllags_R2_025_4cluster.csv"),sep=",",row.names=F)
