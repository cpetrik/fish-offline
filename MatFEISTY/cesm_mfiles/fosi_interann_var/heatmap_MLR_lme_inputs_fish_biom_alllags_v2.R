# Heatmaps of Mult linear regression of forcing on fish biomass
# Includes 2 lags: 0-1 or 1-2
# From (8/21/2023) output
# No zmeso

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

fpval <- read.csv(paste0(datar,"LME_F_mlr_pvals_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
ppval <- read.csv(paste0(datar,"LME_P_mlr_pvals_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
dpval <- read.csv(paste0(datar,"LME_D_mlr_pvals_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
apval <- read.csv(paste0(datar,"LME_A_mlr_pvals_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
bpval <- read.csv(paste0(datar,"LME_B_mlr_pvals_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)



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
# reduce the size of the labels:
adfish <- set(adfish, "labels_cex", 0.8)

some_col_func <- colorspace::diverge_hcl(25)
col_breaks = c(-15,-5,-2,-1.1, # for low
               seq(-1,1,length=21), # for blue
               1.1,2,5,15) # for high

png(paste0(figp,"hclust_wardD_heatmap_LME_A_mlr_coeffs_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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
fdfish <- color_branches(fdfish,k=8) 
# reduce the size of the labels:
fdfish <- set(fdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_F_mlr_coeffs_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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

png(paste0(figp,"hclust_wardD_heatmap_LME_P_mlr_coeffs_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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

png(paste0(figp,"hclust_wardD_heatmap_LME_D_mlr_coeffs_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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
bdfish <- color_branches(bdfish,k=6) #
# reduce the size of the labels:
bdfish <- set(bdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_B_mlr_coeffs_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
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


### ------------------------- Just US --------------------------------------
lid <- c('X54','X1','X2','X65','X10','X3','X5','X6','X7')
lall <- as.data.frame(lid)
lall$lname <- c('CHK','EBS','GAK','AI','HI','CCE','GMX','SE','NE')
names(lall)[1] <- 'LME'

Fcoef2 <- Fcoef
Fcoef2$LME <-rownames(Fcoef)
Fcoef2 <- subset(Fcoef2,Fcoef2$LME %in% lid)
Fcoef2 <- merge(Fcoef2,lall,by='LME')

Pcoef2 <- Pcoef
Pcoef2$LME <-rownames(Pcoef)
Pcoef2 <- subset(Pcoef2,Pcoef2$LME %in% lid)
Pcoef2 <- merge(Pcoef2,lall,by='LME')

Dcoef2 <- Dcoef
Dcoef2$LME <-rownames(Dcoef)
Dcoef2 <- subset(Dcoef2,Dcoef2$LME %in% lid)
Dcoef2 <- merge(Dcoef2,lall,by='LME')

Acoef2 <- Acoef
Acoef2$LME <-rownames(Acoef)
Acoef2 <- subset(Acoef2,Acoef2$LME %in% lid)
Acoef2 <- merge(Acoef2,lall,by='LME')

Bcoef2 <- Bcoef
Bcoef2$LME <-rownames(Bcoef)
Bcoef2 <- subset(Bcoef2,Bcoef2$LME %in% lid)
Bcoef2 <- merge(Bcoef2,lall,by='LME')


### Cluster
Acoef3 <- Acoef2[,2:15]
rownames(Acoef3) <- Acoef2$lname
Acoef3[is.na(Acoef3)] <- 0
d_fish <- dist(Acoef3)
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish) #, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)


png(paste0(figp,"hclust_wardD_heatmap_USLME_A_mlr_coeffs_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Acoef3), 
                  main = "All fish",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
                  #RowSideColors = ccol$col,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(6,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  #add.expr=text(x=0.2, y=-0.1, srt=45, xpd=NA, adj=0, labels="Mag"),
                  col = some_col_func #, breaks=col_breaks
)
dev.off()


## Forage
Fcoef3 <- Fcoef2[,2:15]
rownames(Fcoef3) <- Fcoef2$lname
Fcoef3[is.na(Fcoef3)] <- 0
d_fish <- dist(Fcoef3)
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish) #, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_F_mlr_coeffs_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Fcoef3), 
                  main = "Forage fish",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
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
Pcoef3 <- Pcoef2[,2:15]
rownames(Pcoef3) <- Pcoef2$lname
Pcoef3[is.na(Pcoef3)] <- 0
d_fish <- dist(Pcoef3)
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish) #, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_P_mlr_coeffs_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Pcoef3), 
                  main = "Large Pelagics",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
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
Dcoef3 <- Dcoef2[,2:15]
rownames(Dcoef3) <- Dcoef2$lname
Dcoef3[is.na(Dcoef3)] <- 0
d_fish <- dist(Dcoef3)
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish) #, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_D_mlr_coeffs_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Dcoef3), 
                  main = "Demersals",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
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
Bcoef3 <- Bcoef2[,2:15]
rownames(Bcoef3) <- Bcoef2$lname
Bcoef3[is.na(Bcoef3)] <- 0
d_fish <- dist(Bcoef3)
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish) #, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_B_mlr_coeffs_ALLdiv2SD_alllags_v2.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Bcoef3), 
                  main = "Benthos",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
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



