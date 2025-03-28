# Heatmaps of Mult linear regression of forcing on fish biomass
# Includes all lags up to 2yrs
# Reduced to most sig of each driver and interaction
# From (8/21/2023) output

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
#datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/FOSI/"
#datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"
datar <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/"

fcoef <- read.csv(paste0(datar,"LMEs_mlr_drivers_ALLdiv2SD_alllags_reduc_F.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef <- read.csv(paste0(datar,"LMEs_mlr_drivers_ALLdiv2SD_alllags_reduc_P.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef <- read.csv(paste0(datar,"LMEs_mlr_drivers_ALLdiv2SD_alllags_reduc_D.csv"),sep=",",header = T,stringsAsFactors = F)
acoef <- read.csv(paste0(datar,"LMEs_mlr_drivers_ALLdiv2SD_alllags_reduc_A.csv"),sep=",",header = T,stringsAsFactors = F)
bcoef <- read.csv(paste0(datar,"LMEs_mlr_drivers_ALLdiv2SD_alllags_reduc_B.csv"),sep=",",header = T,stringsAsFactors = F)

names(fcoef)[1] <- "LME"
names(pcoef)[1] <- "LME"
names(dcoef)[1] <- "LME"
names(acoef)[1] <- "LME"
names(bcoef)[1] <- "LME"

### Use heatmap & clustering tree to find patterns

# Cluster LMEs
Acoef <- acoef[,2:9]
rownames(Acoef) <- acoef$LME
Acoef[is.na(Acoef)] <- 0
ad_fish <- dist(Acoef) # method="man" # is a bit better
ahc_fish <- hclust(ad_fish, method = "ward.D")
adfish <- as.dendrogram(ahc_fish)
# Color the branches based on the clusters:
adfish <- color_branches(adfish, k=6) 
# reduce the size of the labels:
adfish <- set(adfish, "labels_cex", 0.8)

adend_6 <- cutree(adfish, k=6)

some_col_func <- colorspace::diverge_hcl(25)
col_breaks = c(-15,-5,-2,-1.1, # for low
               seq(-1,1,length=21), # for blue
               1.1,2,5,15) # for high

png(paste0(figp,"hclust_wardD_heatmap_LME_A_mlr_coeffs_ALLdiv2SD_alllags_reduc.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Acoef), 
                  main = "All fish",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = adfish,
                  #RowSideColors = ccol$col,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
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
Fcoef <- fcoef[,2:9]
rownames(Fcoef) <- fcoef$LME
Fcoef[is.na(Fcoef)] <- 0
fd_fish <- dist(Fcoef) # method="man" # is a bit better
fhc_fish <- hclust(fd_fish, method = "ward.D")
fdfish <- as.dendrogram(fhc_fish)
# Color the branches based on the clusters:
fdfish <- color_branches(fdfish, k=4) #,k=8 
# reduce the size of the labels:
fdfish <- set(fdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_F_mlr_coeffs_ALLdiv2SD_alllags_reduc.png"),    # create PNG for the heat map        
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
                  margins =c(3,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()


## Lg Pel
Pcoef <- pcoef[,2:9]
rownames(Pcoef) <- pcoef$LME
Pcoef[is.na(Pcoef)] <- 0
pd_fish <- dist(Pcoef) # method="man" # is a bit better
phc_fish <- hclust(pd_fish, method = "ward.D")
pdfish <- as.dendrogram(phc_fish)
# Color the branches based on the clusters:
pdfish <- color_branches(pdfish, k=5) #,k=8 
# reduce the size of the labels:
pdfish <- set(pdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_P_mlr_coeffs_ALLdiv2SD_alllags_reduc.png"),    # create PNG for the heat map        
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
                  margins =c(3,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()


## Dem
Dcoef <- dcoef[,2:9]
rownames(Dcoef) <- dcoef$LME
Dcoef[is.na(Dcoef)] <- 0
dd_fish <- dist(Dcoef) # method="man" # is a bit better
dhc_fish <- hclust(dd_fish, method = "ward.D")
ddfish <- as.dendrogram(dhc_fish)
# Color the branches based on the clusters:
ddfish <- color_branches(ddfish, k=4) #,k=8 
# reduce the size of the labels:
ddfish <- set(ddfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_D_mlr_coeffs_ALLdiv2SD_alllags_reduc.png"),    # create PNG for the heat map        
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
                  margins =c(3,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()


## Bent
Bcoef <- bcoef[,2:9]
rownames(Bcoef) <- bcoef$LME
Bcoef[is.na(Bcoef)] <- 0
bd_fish <- dist(Bcoef) # method="man" # is a bit better
bhc_fish <- hclust(bd_fish, method = "ward.D")
bdfish <- as.dendrogram(bhc_fish)
# Color the branches based on the clusters:
bdfish <- color_branches(bdfish, k=5) #,k=8 
# reduce the size of the labels:
bdfish <- set(bdfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_B_mlr_coeffs_ALLdiv2SD_alllags_reduc.png"),    # create PNG for the heat map        
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
                  margins =c(3,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()


### ------------------------- Just US --------------------------------------
lid <- c('54','1','2','65','10','3','5','6','7')
lall <- as.data.frame(lid)
lall$lname <- c('CHK','EBS','GAK','AI','HI','CCE','GMX','SE','NE')
names(lall)[1] <- 'LME'


Fcoef <- subset(fcoef,fcoef$LME %in% lid)
Fcoef <- merge(Fcoef,lall,by='LME')
Fcoef <- Fcoef[,2:10]
mFcoef <- melt(Fcoef,id="lname")
names(mFcoef) <- c("LME","driver","coef")

Pcoef <- subset(pcoef,pcoef$LME %in% lid)
Pcoef <- merge(Pcoef,lall,by='LME')
Pcoef <- Pcoef[,2:10]
mPcoef <- melt(Pcoef,id="lname")
names(mPcoef) <- c("LME","driver","coef")

Dcoef <- subset(dcoef,dcoef$LME %in% lid)
Dcoef <- merge(Dcoef,lall,by='LME')
Dcoef <- Dcoef[,2:10]
mDcoef <- melt(Dcoef,id="lname")
names(mDcoef) <- c("LME","driver","coef")

Acoef <- subset(acoef,acoef$LME %in% lid)
Acoef <- merge(Acoef,lall,by='LME')
Acoef <- Acoef[,2:10]
mAcoef <- melt(Acoef,id="lname")
names(mAcoef) <- c("LME","driver","coef")

Bcoef <- subset(bcoef,bcoef$LME %in% lid)
Bcoef <- merge(Bcoef,lall,by='LME')
Bcoef <- Bcoef[,2:10]
mBcoef <- melt(Bcoef,id="lname")
names(mBcoef) <- c("LME","driver","coef")


f1 <- ggplot(data = mFcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-4.1,4.1), 
                       breaks=c(-3,-2,-1,0,1,2,3),
                       labels=c(-3,-2,-1,0,1,2,3), name="coeff") +
  theme_minimal() + labs(x="") + labs(y="LME") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Forage fish")  

p1 <- ggplot(data = mPcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") + labs(x="") + labs(y="") +
  scale_fill_distiller(palette = "RdBu", limit = c(-4.1,4.1), 
                       breaks=c(-3,-2,-1,0,1,2,3),
                       labels=c(-3,-2,-1,0,1,2,3), name="coeff") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Large pelagics") 

d1 <- ggplot(data = mDcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") + labs(x="Driver") + labs(y="LME") +
  scale_fill_distiller(palette = "RdBu", limit = c(-4.1,4.1), 
                       breaks=c(-3,-2,-1,0,1,2,3),
                       labels=c(-3,-2,-1,0,1,2,3), name="coeff") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Demersals") 

a1 <- ggplot(data = mAcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") + labs(x="Driver") + labs(y="") +
  scale_fill_distiller(palette = "RdBu", limit = c(-4.1,4.1), 
                       breaks=c(-3,-2,-1,0,1,2,3),
                       labels=c(-3,-2,-1,0,1,2,3), name="coeff") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("All fish")  

b1 <- ggplot(data = mBcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") + labs(x="Driver") + labs(y="") +
  scale_fill_distiller(palette = "RdBu", limit = c(-4.1,4.1), 
                       breaks=c(-3,-2,-1,0,1,2,3),
                       labels=c(-3,-2,-1,0,1,2,3), name="coeff") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Benthos") 

png(paste0(figp,'Heatmaps_USLME_mlr_coeffs_ALLdiv2SD_alllags_reduc.png'), 
    width = 12*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
ggarrange(f1,p1,a1,
          d1,b1,
          nrow = 2, ncol = 3,labels = "auto",
          common.legend = TRUE, legend = "right")
dev.off()


### Cluster
Acoef2 <- Acoef[,1:8]
rownames(Acoef2) <- Acoef$lname
Acoef2[is.na(Acoef2)] <- 0
d_fish <- dist(Acoef2)
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish) #, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_A_mlr_coeffs_ALLdiv2SD_alllags_reduc.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Acoef2), 
                  main = "All fish",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
                  #RowSideColors = ccol$col,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
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
Fcoef2 <- Fcoef[,1:8]
rownames(Fcoef2) <- Fcoef$lname
Fcoef2[is.na(Fcoef2)] <- 0
d_fish <- dist(Fcoef2) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish) #, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_F_mlr_coeffs_ALLdiv2SD_alllags_reduc.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Fcoef2), 
                  main = "Forage fish",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()


## Lg Pel
Pcoef2 <- Pcoef[,1:8]
rownames(Pcoef2) <- Pcoef$lname
Pcoef2[is.na(Pcoef2)] <- 0
d_fish <- dist(Pcoef2) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish) #, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_P_mlr_coeffs_ALLdiv2SD_alllags_reduc.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Pcoef2), 
                  main = "Large Pelagics",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()


## Dem
Dcoef2 <- Dcoef[,1:8]
rownames(Dcoef2) <- Dcoef$lname
Dcoef2[is.na(Dcoef2)] <- 0
d_fish <- dist(Dcoef2) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish) #, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_D_mlr_coeffs_ALLdiv2SD_alllags_reduc.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Dcoef2), 
                  main = "Demersals",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()


## Bent
Bcoef2 <- Bcoef[,1:8]
rownames(Bcoef2) <- Bcoef$lname
Bcoef2[is.na(Bcoef2)] <- 0
d_fish <- dist(Bcoef2) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish) #, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_B_mlr_coeffs_ALLdiv2SD_alllags_reduc.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Bcoef2), 
                  main = "Benthos",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()



