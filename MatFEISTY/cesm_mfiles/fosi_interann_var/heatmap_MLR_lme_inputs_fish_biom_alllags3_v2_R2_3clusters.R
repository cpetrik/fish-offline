# Heatmaps of Mult linear regression of forcing on fish biomass
# Includes 3 lags: 0,1,2
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

fcoef <- read.csv(paste0(datar,"LME_F_mlr_coeffs_ALLdiv2SD_alllags3_v2_noint.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef <- read.csv(paste0(datar,"LME_P_mlr_coeffs_ALLdiv2SD_alllags3_v2_noint.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef <- read.csv(paste0(datar,"LME_D_mlr_coeffs_ALLdiv2SD_alllags3_v2_noint.csv"),sep=",",header = T,stringsAsFactors = F)
acoef <- read.csv(paste0(datar,"LME_A_mlr_coeffs_ALLdiv2SD_alllags3_v2_noint.csv"),sep=",",header = T,stringsAsFactors = F)
bcoef <- read.csv(paste0(datar,"LME_B_mlr_coeffs_ALLdiv2SD_alllags3_v2_noint.csv"),sep=",",header = T,stringsAsFactors = F)

lid <- c(1:22,24:32,34:61,63:66)
fcoef$LME <-lid
pcoef$LME <-lid
dcoef$LME <-lid
acoef$LME <-lid
bcoef$LME <-lid


# R2>0.5
fcoef50 <- subset(fcoef, R.2 >= 0.5)
pcoef50 <- subset(pcoef, R.2 >= 0.5)
dcoef50 <- subset(dcoef, R.2 >= 0.5)
acoef50 <- subset(acoef, R.2 >= 0.5)
bcoef50 <- subset(bcoef, R.2 >= 0.5)

### Use heatmap & clustering tree to find patterns

# Cluster LMEs
Acoef <- acoef50[,2:13]
rownames(Acoef) <- acoef50$LME
Acoef[is.na(Acoef)] <- 0
ad_fish <- dist(Acoef) # method="man" # is a bit better
ahc_fish <- hclust(ad_fish, method = "ward.D")
adfish <- as.dendrogram(ahc_fish)
# Color the branches based on the clusters:
adfish <- color_branches(adfish, k=3) 
# reduce the size of the labels:
adfish <- set(adfish, "labels_cex", 0.8)

some_col_func <- colorspace::diverge_hcl(25)
col_breaks = c(-15,-5,-2,-1.1, # for low
               seq(-1,1,length=21), # for blue
               1.1,2,5,15) # for high

png(paste0(figp,"hclust_wardD_heatmap_LME_A_mlr_coeffs_R2_ALLdiv2SD_alllags3_3clus.png"),    # create PNG for the heat map        
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
Fcoef <- fcoef50[,2:13]
rownames(Fcoef) <- fcoef50$LME
Fcoef[is.na(Fcoef)] <- 0
fd_fish <- dist(Fcoef) # method="man" # is a bit better
fhc_fish <- hclust(fd_fish, method = "ward.D")
fdfish <- as.dendrogram(fhc_fish)
# Color the branches based on the clusters:
fdfish <- color_branches(fdfish,k=3) 
# reduce the size of the labels:
fdfish <- set(fdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_F_mlr_coeffs_R2_ALLdiv2SD_alllags3_3clus.png"),    # create PNG for the heat map        
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
Pcoef <- pcoef50[,2:13]
rownames(Pcoef) <- pcoef50$LME
Pcoef[is.na(Pcoef)] <- 0
pd_fish <- dist(Pcoef) # method="man" # is a bit better
phc_fish <- hclust(pd_fish, method = "ward.D")
pdfish <- as.dendrogram(phc_fish)
# Color the branches based on the clusters:
pdfish <- color_branches(pdfish,k=4) #
# reduce the size of the labels:
pdfish <- set(pdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_P_mlr_coeffs_R2_ALLdiv2SD_alllags3_4clus.png"),    # create PNG for the heat map        
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
Dcoef <- dcoef50[,2:13]
rownames(Dcoef) <- dcoef50$LME
Dcoef[is.na(Dcoef)] <- 0
dd_fish <- dist(Dcoef) # method="man" # is a bit better
dhc_fish <- hclust(dd_fish, method = "ward.D")
ddfish <- as.dendrogram(dhc_fish)
# Color the branches based on the clusters:
ddfish <- color_branches(ddfish,k=3) #
# reduce the size of the labels:
ddfish <- set(ddfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_D_mlr_coeffs_R2_ALLdiv2SD_alllags3_3clus.png"),    # create PNG for the heat map        
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


## Bent
Bcoef <- bcoef50[,2:13]
rownames(Bcoef) <- bcoef50$LME
Bcoef[is.na(Bcoef)] <- 0
bd_fish <- dist(Bcoef) # method="man" # is a bit better
bhc_fish <- hclust(bd_fish, method = "ward.D")
bdfish <- as.dendrogram(bhc_fish)
# Color the branches based on the clusters:
bdfish <- color_branches(bdfish,k=3) #
# reduce the size of the labels:
bdfish <- set(bdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_B_mlr_coeffs_R2_ALLdiv2SD_alllags3_3clus.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Bcoef), 
                  main = "Benthos",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = bdfish,
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
adend_k <- cutree(adfish, k=3)
acoef50$Cluster <- adend_k

## Forage
fdend_k <- cutree(fdfish, k=3)
fcoef50$Cluster <- fdend_k

## Lg Pel
pdend_k <- cutree(pdfish, k=4)
pcoef50$Cluster <- pdend_k

## Dem
ddend_k <- cutree(ddfish, k=3)
dcoef50$Cluster <- ddend_k

## Bent
bdend_k <- cutree(bdfish, k=3)
bcoef50$Cluster <- bdend_k


### Save 0.5 (not 0.7)
acoef50 <- acoef50[,c("LME","Cluster")]
acoef <- merge(acoef,acoef50,by="LME",all=T)

fcoef50 <- fcoef50[,c("LME","Cluster")]
fcoef <- merge(fcoef,fcoef50,by="LME",all=T)

pcoef50 <- pcoef50[,c("LME","Cluster")]
pcoef <- merge(pcoef,pcoef50,by="LME",all=T)

dcoef50 <- dcoef50[,c("LME","Cluster")]
dcoef <- merge(dcoef,dcoef50,by="LME",all=T)

bcoef50 <- bcoef50[,c("LME","Cluster")]
bcoef <- merge(bcoef,bcoef50,by="LME",all=T)


write.table(fcoef,paste0(datar,"LME_F_mlr_coeffs_ALLdiv2SD_alllags3_R2_05_3cluster.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datar,"LME_P_mlr_coeffs_ALLdiv2SD_alllags3_R2_05_4cluster.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datar,"LME_D_mlr_coeffs_ALLdiv2SD_alllags3_R2_05_3cluster.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datar,"LME_A_mlr_coeffs_ALLdiv2SD_alllags3_R2_05_3cluster.csv"),sep=",",row.names=F)
write.table(bcoef,paste0(datar,"LME_B_mlr_coeffs_ALLdiv2SD_alllags3_R2_05_3cluster.csv"),sep=",",row.names=F)






