# Heatmaps of Mult linear regression of forcing on fish biomass
# Includes 2 lags: 0-1 or 1-2
# From (8/21/2023) output
# No zmeso
# Reduce to LMEs with high R2 value

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

fcoef <- read.csv(paste0(datar,"LME_F_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef <- read.csv(paste0(datar,"LME_P_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef <- read.csv(paste0(datar,"LME_D_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
acoef <- read.csv(paste0(datar,"LME_A_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
bcoef <- read.csv(paste0(datar,"LME_B_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)

# R2>0.5
fcoef <- subset(fcoef, R.2 >= 0.5)
pcoef <- subset(pcoef, R.2 >= 0.5)
dcoef <- subset(dcoef, R.2 >= 0.5)
acoef <- subset(acoef, R.2 >= 0.5)
bcoef <- subset(bcoef, R.2 >= 0.5)

### Use heatmap & clustering tree to find patterns

# Cluster LMEs
Acoef <- acoef[,2:9]
rownames(Acoef) <- acoef$LME
anm <- c("Det1","Det2","TB1","TB2","TP1","TP2","ZmLoss1","ZmLoss2")
Acoef <- Acoef[,anm]
Acoef[is.na(Acoef)] <- 0
ad_fish <- dist(Acoef) # method="man" # is a bit better
ahc_fish <- hclust(ad_fish, method = "ward.D")
adfish <- as.dendrogram(ahc_fish)
# Color the branches based on the clusters:
adfish <- color_branches(adfish, k=6) 
# reduce the size of the labels:
adfish <- set(adfish, "labels_cex", 0.8)

some_col_func <- colorspace::diverge_hcl(25)
col_breaks = c(-15,-5,-2,-1.1, # for low
               seq(-1,1,length=21), # for blue
               1.1,2,5,15) # for high

png(paste0(figp,"hclust_wardD_heatmap_LME_A_mlr_coeffs_R2_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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
Fcoef <- fcoef[,2:22]
rownames(Fcoef) <- fcoef$LME
fnm <- c("Det0","Det1","TB0","TB1","TP0","TP1","ZmLoss0","ZmLoss1",
         "Det0.TB0","Det0.ZmLoss0","Det1.TB1","Det1.ZmLoss1","TP0.ZmLoss0",  
         "TP1.ZmLoss1")
Fcoef <- Fcoef[,fnm]
Fcoef[is.na(Fcoef)] <- 0
fd_fish <- dist(Fcoef) # method="man" # is a bit better
fhc_fish <- hclust(fd_fish, method = "ward.D")
fdfish <- as.dendrogram(fhc_fish)
# Color the branches based on the clusters:
fdfish <- color_branches(fdfish,k=7) 
# reduce the size of the labels:
fdfish <- set(fdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_F_mlr_coeffs_R2_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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
Pcoef <- pcoef[,2:22]
rownames(Pcoef) <- pcoef$LME
Pcoef <- Pcoef[,anm]
Pcoef[is.na(Pcoef)] <- 0
pd_fish <- dist(Pcoef) # method="man" # is a bit better
phc_fish <- hclust(pd_fish, method = "ward.D")
pdfish <- as.dendrogram(phc_fish)
# Color the branches based on the clusters:
pdfish <- color_branches(pdfish,k=5) #
# reduce the size of the labels:
pdfish <- set(pdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_P_mlr_coeffs_R2_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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
Dcoef <- dcoef[,2:22]
rownames(Dcoef) <- dcoef$LME
Dcoef <- Dcoef[,anm]
Dcoef[is.na(Dcoef)] <- 0
dd_fish <- dist(Dcoef) # method="man" # is a bit better
dhc_fish <- hclust(dd_fish, method = "ward.D")
ddfish <- as.dendrogram(dhc_fish)
# Color the branches based on the clusters:
ddfish <- color_branches(ddfish,k=5) #
# reduce the size of the labels:
ddfish <- set(ddfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_D_mlr_coeffs_R2_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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
Bcoef <- bcoef[,2:22]
rownames(Bcoef) <- bcoef$LME
Bcoef <- Bcoef[,fnm]
Bcoef[is.na(Bcoef)] <- 0
bd_fish <- dist(Bcoef) # method="man" # is a bit better
bhc_fish <- hclust(bd_fish, method = "ward.D")
bdfish <- as.dendrogram(bhc_fish)
# Color the branches based on the clusters:
bdfish <- color_branches(bdfish,k=5) #
# reduce the size of the labels:
bdfish <- set(bdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_B_mlr_coeffs_R2_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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


### TRY R2 >= 0.7? ------------------------------------------------------------------------------------------

fcoef <- subset(fcoef, R.2 >= 0.7)
pcoef <- subset(pcoef, R.2 >= 0.7)
dcoef <- subset(dcoef, R.2 >= 0.7)
acoef <- subset(acoef, R.2 >= 0.7)
bcoef <- subset(bcoef, R.2 >= 0.7)

### Use heatmap & clustering tree to find patterns

# Cluster LMEs
Acoef <- acoef[,2:22]
rownames(Acoef) <- acoef$LME
anm <- c("Det1","Det2","TB1","TB2","TP1","TP2","ZmLoss1","ZmLoss2",
         "Det1.TB1","Det1.ZmLoss1","Det2.TB2","Det2.ZmLoss2","TP1.ZmLoss1",  
         "TP2.ZmLoss2")
Acoef <- Acoef[,anm]
Acoef[is.na(Acoef)] <- 0
ad_fish <- dist(Acoef) # method="man" # is a bit better
ahc_fish <- hclust(ad_fish, method = "ward.D")
adfish <- as.dendrogram(ahc_fish)
# Color the branches based on the clusters:
adfish <- color_branches(adfish, k=6) 
adend_6 <- cutree(adfish, k=6)
acoef$Cluster <- adend_6
# reduce the size of the labels:
adfish <- set(adfish, "labels_cex", 0.8)

some_col_func <- colorspace::diverge_hcl(25)
col_breaks = c(-15,-5,-2,-1.1, # for low
               seq(-1,1,length=21), # for blue
               1.1,2,5,15) # for high

png(paste0(figp,"hclust_wardD_heatmap_LME_A_mlr_coeffs_R2_07_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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
Fcoef <- fcoef[,2:22]
rownames(Fcoef) <- fcoef$LME
fnm <- c("Det0","Det1","TB0","TB1","TP0","TP1","ZmLoss0","ZmLoss1",
         "Det0.TB0","Det0.ZmLoss0","Det1.TB1","Det1.ZmLoss1","TP0.ZmLoss0",  
         "TP1.ZmLoss1")
Fcoef <- Fcoef[,fnm]
Fcoef[is.na(Fcoef)] <- 0
fd_fish <- dist(Fcoef) # method="man" # is a bit better
fhc_fish <- hclust(fd_fish, method = "ward.D")
fdfish <- as.dendrogram(fhc_fish)
# Color the branches based on the clusters:
fdfish <- color_branches(fdfish,k=7) 
fdend_7 <- cutree(fdfish, k=7)
fcoef$Cluster <- fdend_7
# reduce the size of the labels:
fdfish <- set(fdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_F_mlr_coeffs_R2_07_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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
Pcoef <- pcoef[,2:22]
rownames(Pcoef) <- pcoef$LME
Pcoef <- Pcoef[,anm]
Pcoef[is.na(Pcoef)] <- 0
pd_fish <- dist(Pcoef) # method="man" # is a bit better
phc_fish <- hclust(pd_fish, method = "ward.D")
pdfish <- as.dendrogram(phc_fish)
# Color the branches based on the clusters:
pdfish <- color_branches(pdfish,k=5) #
pdend_5 <- cutree(pdfish, k=5)
pcoef$Cluster <- pdend_5
# reduce the size of the labels:
pdfish <- set(pdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_P_mlr_coeffs_R2_07_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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
Dcoef <- dcoef[,2:22]
rownames(Dcoef) <- dcoef$LME
Dcoef <- Dcoef[,anm]
Dcoef[is.na(Dcoef)] <- 0
dd_fish <- dist(Dcoef) # method="man" # is a bit better
dhc_fish <- hclust(dd_fish, method = "ward.D")
ddfish <- as.dendrogram(dhc_fish)
# Color the branches based on the clusters:
ddfish <- color_branches(ddfish,k=7) #
ddend_7 <- cutree(ddfish, k=7)
dcoef$Cluster <- ddend_7
# reduce the size of the labels:
ddfish <- set(ddfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_D_mlr_coeffs_R2_07_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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
Bcoef <- bcoef[,2:22]
rownames(Bcoef) <- bcoef$LME
Bcoef <- Bcoef[,fnm]
Bcoef[is.na(Bcoef)] <- 0
bd_fish <- dist(Bcoef) # method="man" # is a bit better
bhc_fish <- hclust(bd_fish, method = "ward.D")
bdfish <- as.dendrogram(bhc_fish)
# Color the branches based on the clusters:
bdfish <- color_branches(bdfish,k=5) #
bdend_5 <- cutree(bdfish, k=5)
bcoef$Cluster <- bdend_5
# reduce the size of the labels:
bdfish <- set(bdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_B_mlr_coeffs_R2_07_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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


write.table(fcoef,paste0(datar,"LME_F_mlr_coeffs_ALLdiv2SD_alllags_R2_07_cluster.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datar,"LME_P_mlr_coeffs_ALLdiv2SD_alllags_R2_07_cluster.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datar,"LME_D_mlr_coeffs_ALLdiv2SD_alllags_R2_07_cluster.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datar,"LME_A_mlr_coeffs_ALLdiv2SD_alllags_R2_07_cluster.csv"),sep=",",row.names=F)
write.table(bcoef,paste0(datar,"LME_B_mlr_coeffs_ALLdiv2SD_alllags_R2_07_cluster.csv"),sep=",",row.names=F)
