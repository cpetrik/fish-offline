# Heatmaps of Mult linear regression of forcing on fish nu
# From (7/2023) output

rm(list=ls())

library(ggplot2)
library(cowplot) #plot_grid
library(gridExtra)
library(ggpubr) #cowplot and gridExtra
library(gplots)

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/corrs/"

### ------------------------------------------------------------
# load data
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/FOSI/"
datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"
fcoef <- read.csv(paste0(datar,"LME_Fnu_mlr_coeffs_ALLdiv2SD.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef <- read.csv(paste0(datar,"LME_Pnu_mlr_coeffs_ALLdiv2SD.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef <- read.csv(paste0(datar,"LME_Dnu_mlr_coeffs_ALLdiv2SD.csv"),sep=",",header = T,stringsAsFactors = F)
acoef <- read.csv(paste0(datar,"LME_Anu_mlr_coeffs_ALLdiv2SD.csv"),sep=",",header = T,stringsAsFactors = F)

fpval <- read.csv(paste0(datar,"LME_Fnu_mlr_pvals_ALLdiv2SD.csv"),sep=",",header = T,stringsAsFactors = F)
ppval <- read.csv(paste0(datar,"LME_Pnu_mlr_pvals_ALLdiv2SD.csv"),sep=",",header = T,stringsAsFactors = F)
dpval <- read.csv(paste0(datar,"LME_Dnu_mlr_pvals_ALLdiv2SD.csv"),sep=",",header = T,stringsAsFactors = F)
apval <- read.csv(paste0(datar,"LME_Anu_mlr_pvals_ALLdiv2SD.csv"),sep=",",header = T,stringsAsFactors = F)


# Cluster LMEs
library(rpart)
library(tree)
library(pvclust)
library(mclust)
library(dendextend)

Acoef <- acoef[,2:9]
rownames(Acoef) <- acoef$LME
Acoef[is.na(Acoef)] <- 0
d_fish <- dist(Acoef) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish, k=8) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

some_col_func <- colorspace::diverge_hcl(25)
col_breaks = c(-15,-5,-2,-1.1, # for low
               seq(-1,1,length=21), # for blue
               1.1,2,5,15) # for high

png(paste0(figp,"hclust_wardD_heatmap_LME_Anu_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Acoef), 
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
Fcoef <- fcoef[,2:9]
rownames(Fcoef) <- fcoef$LME
Fcoef[is.na(Fcoef)] <- 0
d_fish <- dist(Fcoef) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish, k=8) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_Fnu_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Fcoef), 
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
Pcoef <- pcoef[,2:9]
rownames(Pcoef) <- pcoef$LME
Pcoef[is.na(Pcoef)] <- 0
d_fish <- dist(Pcoef) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish, k=8) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_Pnu_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Pcoef), 
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
Dcoef <- dcoef[,2:9]
rownames(Dcoef) <- dcoef$LME
Dcoef[is.na(Dcoef)] <- 0
d_fish <- dist(Dcoef) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish, k=8) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_Dnu_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Dcoef), 
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



### ------------------------- Just US --------------------------------------
lid <- c('LME54','LME1','LME2','LME65','LME10','LME3','LME5','LME6','LME7')
lall <- as.data.frame(lid)
lall$lname <- c('CHK','EBS','GAK','AI','HI','CCE','GMX','SE','NE')
names(lall)[1] <- 'LME'

hcol <- c("Det","TB","TP","Zmeso","ZmLoss","Det.TB","TP.Zmeso","TP.ZmLoss","lname")

Fcoef <- subset(fcoef,fcoef$LME %in% lid)
Fcoef <- merge(Fcoef,lall,by='LME')
Fcoef <- Fcoef[,hcol]
mFcoef <- melt(Fcoef,id="lname")
names(mFcoef) <- c("LME","driver","coef")
Fpval <- merge(lall,fpval,by='LME',all=F)
Fpval <- Fpval[,hcol]
mFpval <- melt(Fpval,id="lname")
names(mFpval) <- c("LME","driver","pval")
mFpval$sym <- mFpval$pval
mFpval$sym[mFpval$pval <= 0.05] <- "*"
mFpval$sym[mFpval$pval > 0.05] <- NA
mFcoef$sym <- mFpval$sym

Pcoef <- subset(pcoef,pcoef$LME %in% lid)
Pcoef <- merge(Pcoef,lall,by='LME')
Pcoef <- Pcoef[,hcol]
mPcoef <- melt(Pcoef,id="lname")
names(mPcoef) <- c("LME","driver","coef")
Ppval <- merge(lall,ppval,by='LME',all=F)
Ppval <- Ppval[,hcol]
mPpval <- melt(Ppval,id="lname")
names(mPpval) <- c("LME","driver","pval")
mPpval$sym <- mPpval$pval
mPpval$sym[mPpval$pval <= 0.05] <- "*"
mPpval$sym[mPpval$pval > 0.05] <- NA
mPcoef$sym <- mPpval$sym

Dcoef <- subset(dcoef,dcoef$LME %in% lid)
Dcoef <- merge(Dcoef,lall,by='LME')
Dcoef <- Dcoef[,hcol]
mDcoef <- melt(Dcoef,id="lname")
names(mDcoef) <- c("LME","driver","coef")
Dpval <- merge(lall,dpval,by='LME',all=F)
Dpval <- Dpval[,hcol]
mDpval <- melt(Dpval,id="lname")
names(mDpval) <- c("LME","driver","pval")
mDpval$sym <- mDpval$pval
mDpval$sym[mDpval$pval <= 0.05] <- "*"
mDpval$sym[mDpval$pval > 0.05] <- NA
mDcoef$sym <- mDpval$sym

Acoef <- subset(acoef,acoef$LME %in% lid)
Acoef <- merge(Acoef,lall,by='LME')
Acoef <- Acoef[,hcol]
mAcoef <- melt(Acoef,id="lname")
names(mAcoef) <- c("LME","driver","coef")
Apval <- merge(lall,apval,by='LME',all=F)
Apval <- Apval[,hcol]
mApval <- melt(Apval,id="lname")
names(mApval) <- c("LME","driver","pval")
mApval$sym <- mApval$pval
mApval$sym[mApval$pval <= 0.05] <- "*"
mApval$sym[mApval$pval > 0.05] <- NA
mAcoef$sym <- mApval$sym


f1 <- ggplot(data = mFcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-4.8,4.8), 
                       breaks=c(-3,-2,-1,0,1,2,3),
                       labels=c(-3,-2,-1,0,1,2,3), name="coeff") +
  theme_minimal() + labs(x="") + labs(y="LME") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Forage fish") +
  geom_text(aes(driver, LME, label = sym), color = "black", size = 5) 

p1 <- ggplot(data = mPcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") + labs(x="") + labs(y="") +
  scale_fill_distiller(palette = "RdBu", limit = c(-4.8,4.8), 
                       breaks=c(-3,-2,-1,0,1,2,3),
                       labels=c(-3,-2,-1,0,1,2,3), name="coeff") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Large pelagics") +
  geom_text(aes(driver, LME, label = sym), color = "black", size = 5) 

d1 <- ggplot(data = mDcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") + labs(x="Driver") + labs(y="LME") +
  scale_fill_distiller(palette = "RdBu", limit = c(-4.8,4.8), 
                       breaks=c(-3,-2,-1,0,1,2,3),
                       labels=c(-3,-2,-1,0,1,2,3), name="coeff") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Demersals") +
  geom_text(aes(driver, LME, label = sym), color = "black", size = 5) 

a1 <- ggplot(data = mAcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") + labs(x="Driver") + labs(y="") +
  scale_fill_distiller(palette = "RdBu", limit = c(-4.8,4.8), 
                       breaks=c(-3,-2,-1,0,1,2,3),
                       labels=c(-3,-2,-1,0,1,2,3), name="coeff") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("All fish") +
  geom_text(aes(driver, LME, label = sym), color = "black", size = 5) 



png(paste0(figp,'Heatmaps_USLME_nu_mlr_coeffs_ALLdiv2SD.png'), 
    width = 6*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
ggarrange(f1,p1,
          a1,d1,
          nrow = 2, ncol = 2,labels = "auto",
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
dfish <- color_branches(dfish, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

some_col_func <- colorspace::diverge_hcl(25)
col_breaks = c(-5,-2,-1,-0.5, # for low
               seq(-0.4,0.4,length=21), # for blue
               0.5,1,2,5) # for high
# test <- colorspace::diverge_hcl(28)
# hv <- heatmap.2(as.matrix(tvec), 
#                 trace="none",          
#                 col = some_col_func)

png(paste0(figp,"hclust_wardD_heatmap_USLME_Anu_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
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
dfish <- color_branches(dfish, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_Fnu_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
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
dfish <- color_branches(dfish, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_Pnu_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
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
dfish <- color_branches(dfish, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_Dnu_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
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



