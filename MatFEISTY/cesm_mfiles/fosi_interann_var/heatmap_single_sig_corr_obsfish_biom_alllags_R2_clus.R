# Heatmaps of Pearson corrs of single forcing on fish biomass
# Simulation with observed fishing effort
# Includes 3 lags: 0,1,2
# From (7/23/2024) output
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
figp <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/corrs/"

### ------------------------------------------------------------
# load data
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/FOSI/"
datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"
#datar <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/"

fcoef <- read.csv(paste0(datar,"LME_sig_corr_maxlag_driver_biom_obsfish_F.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef <- read.csv(paste0(datar,"LME_sig_corr_maxlag_driver_biom_obsfish_P.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef <- read.csv(paste0(datar,"LME_sig_corr_maxlag_driver_biom_obsfish_D.csv"),sep=",",header = T,stringsAsFactors = F)
acoef <- read.csv(paste0(datar,"LME_sig_corr_maxlag_driver_biom_obsfish_A.csv"),sep=",",header = T,stringsAsFactors = F)


### Use heatmap & clustering tree to find patterns

### p-val <= 0.5 ------------------------------------------------------------------------------------------
# Remove interior seas LMEs
iid <- c(23,33,62)
## Need to save p value?
fcoef50 <- fcoef[-iid,]
pcoef50 <- pcoef[-iid,]
dcoef50 <- dcoef[-iid,]
acoef50 <- acoef[-iid,]

fcoef50[is.na(fcoef50)] <- 0
pcoef50[is.na(pcoef50)] <- 0
dcoef50[is.na(dcoef50)] <- 0
acoef50[is.na(acoef50)] <- 0


Acoef <- acoef50[,1:4]
rownames(Acoef) <- acoef50$LME
Acoef[is.na(Acoef)] <- 0
ad_fish <- dist(Acoef) # method="man" # is a bit better
ahc_fish <- hclust(ad_fish, method = "ward.D2")
adfish <- as.dendrogram(ahc_fish)
# Color the branches based on the clusters:
adfish <- color_branches(adfish, k=3) 
adend_6 <- cutree(adfish, k=3)
acoef50$Cluster <- adend_6
# reduce the size of the labels:
adfish <- set(adfish, "labels_cex", 0.8)

some_col_func <- colorspace::diverge_hcl(25)

png(paste0(figp,"heatmap_LME_driver_biom_obsfish_A_sig_corr_coeffs_3clus.png"),    # create PNG for the heat map        
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

png(paste0(figp,"heatmap_LME_driver_biom_obsfish_F_sig_corr_coeffs_4clus.png"),    # create PNG for the heat map        
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
pdfish <- color_branches(pdfish,k=5) #
pdend_5 <- cutree(pdfish, k=5)
pcoef50$Cluster <- pdend_5
# reduce the size of the labels:
pdfish <- set(pdfish, "labels_cex", 0.8)

png(paste0(figp,"heatmap_LME_driver_biom_obsfish_P_sig_corr_coeffs_5clus.png"),    # create PNG for the heat map        
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

png(paste0(figp,"heatmap_LME_driver_biom_obsfish_D_sig_corr_coeffs_4clus.png"),    # create PNG for the heat map        
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
write.table(acoef,paste0(datar,"LME_driver_biom_obsfish_A_sig_corr_coeffs_cluster.csv"),sep=",",row.names=F)
write.table(fcoef,paste0(datar,"LME_driver_biom_obsfish_F_sig_corr_coeffs_cluster.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datar,"LME_driver_biom_obsfish_P_sig_corr_coeffs_cluster.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datar,"LME_driver_biom_obsfish_D_sig_corr_coeffs_cluster.csv"),sep=",",row.names=F)


