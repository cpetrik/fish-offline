# Heatmaps of Pearson corrs of single forcing on fish "recruitment" (gamma)
# Includes 3 lags: 0,1,2
# From (2/7/2024) output
# No zmeso
# Reduce to lag of covar with highest coeff ?!?!
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

fcoef <- read.csv(paste0(datar,"LME_corr_maxlag_driver_biom_F.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef <- read.csv(paste0(datar,"LME_corr_maxlag_driver_biom_P.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef <- read.csv(paste0(datar,"LME_corr_maxlag_driver_biom_D.csv"),sep=",",header = T,stringsAsFactors = F)
acoef <- read.csv(paste0(datar,"LME_corr_maxlag_driver_biom_A.csv"),sep=",",header = T,stringsAsFactors = F)


### Use heatmap & clustering tree to find patterns

### TRY R2 >= 0.5 ------------------------------------------------------------------------------------------
## Need to save p value?
# fcoef50 <- subset(fcoef, R2 >= 0.75)
# pcoef50 <- subset(pcoef, R2 >= 0.75)
# dcoef50 <- subset(dcoef, R2 >= 0.75)
# acoef50 <- subset(acoef, R2 >= 0.75)
fcoef50 <- na.omit(fcoef)
pcoef50 <- na.omit(pcoef)
dcoef50 <- na.omit(dcoef)
acoef50 <- na.omit(acoef)


Acoef <- acoef50[,1:4]
rownames(Acoef) <- acoef50$LME
Acoef[is.na(Acoef)] <- 0
ad_fish <- dist(Acoef) # method="man" # is a bit better
ahc_fish <- hclust(ad_fish, method = "ward.D2")
adfish <- as.dendrogram(ahc_fish)
# Color the branches based on the clusters:
adfish <- color_branches(adfish, k=4) 
adend_6 <- cutree(adfish, k=4)
acoef50$Cluster <- adend_6
# reduce the size of the labels:
adfish <- set(adfish, "labels_cex", 0.8)

some_col_func <- colorspace::diverge_hcl(25)

png(paste0(figp,"heatmap_LME_driver_biom_A_corr_coeffs_4clus.png"),    # create PNG for the heat map        
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
Fcoef <- fcoef50[,1:4]
rownames(Fcoef) <- fcoef50$LME
Fcoef[is.na(Fcoef)] <- 0
fd_fish <- dist(Fcoef) # method="man" # is a bit better
fhc_fish <- hclust(fd_fish, method = "ward.D2")
fdfish <- as.dendrogram(fhc_fish)
# Color the branches based on the clusters:
fdfish <- color_branches(fdfish,k=4) 
fdend_6 <- cutree(fdfish, k=4)
fcoef50$Cluster <- fdend_6
# reduce the size of the labels:
fdfish <- set(fdfish, "labels_cex", 0.8)

png(paste0(figp,"heatmap_LME_driver_biom_F_corr_coeffs_4clus.png"),    # create PNG for the heat map        
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
Pcoef <- pcoef50[,1:4]
rownames(Pcoef) <- pcoef50$LME
Pcoef[is.na(Pcoef)] <- 0
pd_fish <- dist(Pcoef) # method="man" # is a bit better
phc_fish <- hclust(pd_fish, method = "ward.D2")
pdfish <- as.dendrogram(phc_fish)
# Color the branches based on the clusters:
pdfish <- color_branches(pdfish,k=4) #
pdend_5 <- cutree(pdfish, k=4)
pcoef50$Cluster <- pdend_5
# reduce the size of the labels:
pdfish <- set(pdfish, "labels_cex", 0.8)

png(paste0(figp,"heatmap_LME_driver_biom_P_corr_coeffs_4clus.png"),    # create PNG for the heat map        
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
Dcoef <- dcoef50[,1:4]
rownames(Dcoef) <- dcoef50$LME
Dcoef[is.na(Dcoef)] <- 0
dd_fish <- dist(Dcoef) # method="man" # is a bit better
dhc_fish <- hclust(dd_fish, method = "ward.D2")
ddfish <- as.dendrogram(dhc_fish)
# Color the branches based on the clusters:
ddfish <- color_branches(ddfish,k=4) #
ddend_7 <- cutree(ddfish, k=4)
dcoef50$Cluster <- ddend_7
# reduce the size of the labels:
ddfish <- set(ddfish, "labels_cex", 0.8)

png(paste0(figp,"heatmap_LME_driver_biom_D_corr_coeffs_4clus.png"),    # create PNG for the heat map        
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


acoef50 <- acoef50[,c("LME","Cluster")]
acoef <- merge(acoef,acoef50,by="LME",all=T)

fcoef50 <- fcoef50[,c("LME","Cluster")]
fcoef <- merge(fcoef,fcoef50,by="LME",all=T)

pcoef50 <- pcoef50[,c("LME","Cluster")]
pcoef <- merge(pcoef,pcoef50,by="LME",all=T)

dcoef50 <- dcoef50[,c("LME","Cluster")]
dcoef <- merge(dcoef,dcoef50,by="LME",all=T)

### Save clustering
write.table(acoef,paste0(datar,"LME_driver_biom_A_corr_coeffs_4cluster.csv"),sep=",",row.names=F)
write.table(fcoef,paste0(datar,"LME_driver_biom_F_corr_coeffs_4cluster.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datar,"LME_driver_biom_P_corr_coeffs_4cluster.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datar,"LME_driver_biom_D_corr_coeffs_4cluster.csv"),sep=",",row.names=F)


