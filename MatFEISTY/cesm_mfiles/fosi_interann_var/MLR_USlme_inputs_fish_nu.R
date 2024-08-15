# Mult linear regression of forcing on fish nu
# Just US LMES
# Giving different results then when all LMEs done together
# STOP using

rm(list=ls())

library("MuMIn")


source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/corrs/"

### ------------------------------------------------------------
# load data
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/FOSI/"
TP <- read.csv(paste0(datap,"LMEs_anoms_div2SD_TP.csv"),sep=",",header = T,stringsAsFactors = F)
TB <- read.csv(paste0(datap,"LMEs_anoms_div2SD_TB.csv"),sep=",",header = T,stringsAsFactors = F)
Det <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Det.csv"),sep=",",header = T,stringsAsFactors = F)
Zmeso <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Zmeso.csv"),sep=",",header = T,stringsAsFactors = F)
ZmLoss <- read.csv(paste0(datap,"LMEs_anoms_div2SD_ZmLoss.csv"),sep=",",header = T,stringsAsFactors = F)
FF <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Fnu.csv"),sep=",",header = T,stringsAsFactors = F)
P <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Pnu.csv"),sep=",",header = T,stringsAsFactors = F)
D <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Dnu.csv"),sep=",",header = T,stringsAsFactors = F)
A <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Anu.csv"),sep=",",header = T,stringsAsFactors = F)

datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"
Flag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_Fnu.csv"),sep=",",header = T,stringsAsFactors = F)
Plag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_Pnu.csv"),sep=",",header = T,stringsAsFactors = F)
Dlag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_Dnu.csv"),sep=",",header = T,stringsAsFactors = F)
Alag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_Anu.csv"),sep=",",header = T,stringsAsFactors = F)


### MLR of drivers ---------------------------------------------------------
dname <- c('TP','TB','Zmeso','ZmLoss','Det')

#Start with US LMEs
lid <- c(54,1:2,65,10,3,5:7)
cname <- c('LME54','LME1','LME2','LME65','LME10','LME3','LME5','LME6','LME7')
lname <- c('CHK','EBS','GAK','AI','HI','CCE','GMX','SE','NE')
nlme <- length(lid)

#Loop over LMEs 
for (i in 1:nlme) {
  id <- lid[i]
  ln <- lname[i]
  cn <- cname[i]
  
  ### Drivers
  drive <- TP[,c('Year',cn)]
  names(drive) <- c('Year','TP')
  drive$TB <- TB[,cn]
  drive$Det <- Det[,cn]
  drive$Zmeso <- Zmeso[,cn]
  drive$ZmLoss <- ZmLoss[,cn]
  
  
  ### F ----------------------------------------------------------------
  maxf <- max(Flag[i,])
  yst <- 1+maxf
  yen <- 68+maxf
  ffish <- data.frame(matrix(ncol = 6, nrow = yen))
  names(ffish) <- c('Fish','TP','TB','Det','Zmeso','ZmLoss')
  
  ffish$TP[(yst-Flag$TP[i]):(yen-Flag$TP[i])] <- drive$TP 
  ffish$TB[(yst-Flag$TB[i]):(yen-Flag$TB[i])] <- drive$TB 
  ffish$Det[(yst-Flag$Det[i]):(yen-Flag$Det[i])] <- drive$Det
  ffish$Zmeso[(yst-Flag$Zmeso[i]):(yen-Flag$Zmeso[i])] <- drive$Zmeso 
  ffish$ZmLoss[(yst-Flag$ZmLoss[i]):(yen-Flag$ZmLoss[i])] <- drive$ZmLoss
  ffish$Fish[(yst):(yen)] <- FF[,cn]
  ffish <- na.omit(ffish)
  
  options(na.action = "na.fail")
  fmod <- lm(Fish ~ TP + TB + Det + Zmeso + ZmLoss + TP*Zmeso + TP*ZmLoss + TB*Det, data=ffish)
  summary(fmod)
  fcombo <- dredge(fmod)
  ## Create arrays for coefficients & p-vals
  if (i==1) {
    fcoef <- data.frame(matrix(ncol = length(fcombo), nrow = nlme))
    names(fcoef) <- names(fcombo)
    fpval <- fcoef[,1:9]
    fcoef$LME <- lname
    fpval$LME <- lname
    pcoef <- fcoef
    ppval <- fpval
    dcoef <- fcoef
    dpval <- fpval
    acoef <- fcoef
    apval <- fpval
    bcoef <- fcoef
    bpval <- fpval
  }
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (fcombo)[1]
  fcoef[i, match(names(xx), colnames(fcoef))] = xx
  y <- summary(get.models(fcombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  fpval[i, match(names(yy), colnames(fpval))] = yy
  
  ### P ----------------------------------------------------------------
  maxp <- max(Plag[i,])
  yst <- 1+maxp
  yen <- 68+maxp
  pfish <- data.frame(matrix(ncol = 6, nrow = yen))
  names(pfish) <- c('Fish','TP','TB','Det','Zmeso','ZmLoss')
  
  pfish$TP[(yst-Plag$TP[i]):(yen-Plag$TP[i])] <- drive$TP 
  pfish$TB[(yst-Plag$TB[i]):(yen-Plag$TB[i])] <- drive$TB 
  pfish$Det[(yst-Plag$Det[i]):(yen-Plag$Det[i])] <- drive$Det
  pfish$Zmeso[(yst-Plag$Zmeso[i]):(yen-Plag$Zmeso[i])] <- drive$Zmeso 
  pfish$ZmLoss[(yst-Plag$ZmLoss[i]):(yen-Plag$ZmLoss[i])] <- drive$ZmLoss
  pfish$Fish[(yst):(yen)] <- P[,cn]
  pfish <- na.omit(pfish)
  
  pmod <- lm(Fish ~ TP + TB + Det + Zmeso + ZmLoss + TP*Zmeso + TP*ZmLoss + TB*Det, data=pfish)
  pcombo <- dredge(pmod)
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (pcombo)[1]
  pcoef[i, match(names(xx), colnames(pcoef))] = xx
  y <- summary(get.models(pcombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  ppval[i, match(names(yy), colnames(ppval))] = yy
  
  
  ### D ----------------------------------------------------------------
  maxd <- max(Dlag[i,])
  yst <- 1+maxd
  yen <- 68+maxd
  dfish <- data.frame(matrix(ncol = 6, nrow = yen))
  names(dfish) <- c('Fish','TP','TB','Det','Zmeso','ZmLoss')
  
  dfish$TP[(yst-Dlag$TP[i]):(yen-Dlag$TP[i])] <- drive$TP 
  dfish$TB[(yst-Dlag$TB[i]):(yen-Dlag$TB[i])] <- drive$TB 
  dfish$Det[(yst-Dlag$Det[i]):(yen-Dlag$Det[i])] <- drive$Det
  dfish$Zmeso[(yst-Dlag$Zmeso[i]):(yen-Dlag$Zmeso[i])] <- drive$Zmeso 
  dfish$ZmLoss[(yst-Dlag$ZmLoss[i]):(yen-Dlag$ZmLoss[i])] <- drive$ZmLoss
  dfish$Fish[(yst):(yen)] <- D[,cn]
  dfish <- na.omit(dfish)
  
  dmod <- lm(Fish ~ TP + TB + Det + Zmeso + ZmLoss + TP*Zmeso + TP*ZmLoss + TB*Det, data=dfish)
  dcombo <- dredge(dmod)
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (dcombo)[1]
  dcoef[i, match(names(xx), colnames(dcoef))] = xx
  y <- summary(get.models(dcombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  dpval[i, match(names(yy), colnames(dpval))] = yy
  
  ### A ----------------------------------------------------------------
  maxa <- max(Alag[i,])
  yst <- 1+maxa
  yen <- 68+maxa
  afish <- data.frame(matrix(ncol = 6, nrow = yen))
  names(afish) <- c('Fish','TP','TB','Det','Zmeso','ZmLoss')
  
  afish$TP[(yst-Alag$TP[i]):(yen-Alag$TP[i])] <- drive$TP 
  afish$TB[(yst-Alag$TB[i]):(yen-Alag$TB[i])] <- drive$TB 
  afish$Det[(yst-Alag$Det[i]):(yen-Alag$Det[i])] <- drive$Det
  afish$Zmeso[(yst-Alag$Zmeso[i]):(yen-Alag$Zmeso[i])] <- drive$Zmeso 
  afish$ZmLoss[(yst-Alag$ZmLoss[i]):(yen-Alag$ZmLoss[i])] <- drive$ZmLoss
  afish$Fish[(yst):(yen)] <- A[,cn]
  afish <- na.omit(afish)
  
  amod <- lm(Fish ~ TP + TB + Det + Zmeso + ZmLoss + TP*Zmeso + TP*ZmLoss + TB*Det, data=afish)
  acombo <- dredge(amod)
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (acombo)[1]
  acoef[i, match(names(xx), colnames(acoef))] = xx
  y <- summary(get.models(acombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  apval[i, match(names(yy), colnames(apval))] = yy
  
  
} # LMEs

write.table(fcoef,paste0(datar,"USlme_Fnu_mlr_coeffs_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datar,"USlme_Pnu_mlr_coeffs_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datar,"USlme_Dnu_mlr_coeffs_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datar,"USlme_Anu_mlr_coeffs_ALLdiv2SD.csv"),sep=",",row.names=F)

write.table(fpval,paste0(datar,"USlme_Fnu_mlr_pvals_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(ppval,paste0(datar,"USlme_Pnu_mlr_pvals_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(dpval,paste0(datar,"USlme_Dnu_mlr_pvals_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(apval,paste0(datar,"USlme_Anu_mlr_pvals_ALLdiv2SD.csv"),sep=",",row.names=F)



### Would heatmap or something help find patterns?
library(ggplot2)
library(cowplot) #plot_grid
library(gridExtra)
library(ggpubr) #cowplot and gridExtra
library(gplots)
library(reshape2)

Fcoef <- fcoef[,2:9]
Fcoef$LME <- fcoef$LME
mFcoef <- melt(Fcoef,id="LME")
names(mFcoef) <- c("LME","driver","coef")
Fpval <- fpval[,2:9]
Fpval$LME <- fpval$LME
mFpval <- melt(Fpval,id="LME")
names(mFpval) <- c("LME","driver","pval")
mFpval$sym <- mFpval$pval
mFpval$sym[mFpval$pval <= 0.05] <- "*"
mFpval$sym[mFpval$pval > 0.05] <- NA
mFcoef$sym <- mFpval$sym

Pcoef <- pcoef[,2:9]
Pcoef$LME <- pcoef$LME
mPcoef <- melt(Pcoef,id="LME")
names(mPcoef) <- c("LME","driver","coef")
Ppval <- ppval[,2:9]
Ppval$LME <- ppval$LME
mPpval <- melt(Ppval,id="LME")
names(mPpval) <- c("LME","driver","pval")
mPpval$sym <- mPpval$pval
mPpval$sym[mPpval$pval <= 0.05] <- "*"
mPpval$sym[mPpval$pval > 0.05] <- NA
mPcoef$sym <- mPpval$sym

Dcoef <- dcoef[,2:9]
Dcoef$LME <- dcoef$LME
mDcoef <- melt(Dcoef,id="LME")
names(mDcoef) <- c("LME","driver","coef")
Dpval <- dpval[,2:9]
Dpval$LME <- dpval$LME
mDpval <- melt(Dpval,id="LME")
names(mDpval) <- c("LME","driver","pval")
mDpval$sym <- mDpval$pval
mDpval$sym[mDpval$pval <= 0.05] <- "*"
mDpval$sym[mDpval$pval > 0.05] <- NA
mDcoef$sym <- mDpval$sym

Acoef <- acoef[,2:9]
Acoef$LME <- acoef$LME
mAcoef <- melt(Acoef,id="LME")
names(mAcoef) <- c("LME","driver","coef")
Apval <- apval[,2:9]
Apval$LME <- apval$LME
mApval <- melt(Apval,id="LME")
names(mApval) <- c("LME","driver","pval")
mApval$sym <- mApval$pval
mApval$sym[mApval$pval <= 0.05] <- "*"
mApval$sym[mApval$pval > 0.05] <- NA
mAcoef$sym <- mApval$sym


f1 <- ggplot(data = mFcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-2.2,2.2), 
                       breaks=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5),
                       labels=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5), name="coeff") +
  theme_minimal() + labs(x="") + labs(y="LME") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Forage fish") +
  geom_text(aes(driver, LME, label = sym), color = "black", size = 5) 

p1 <- ggplot(data = mPcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") + labs(x="") + labs(y="") +
  scale_fill_distiller(palette = "RdBu", limit = c(-2.2,2.2), 
                       breaks=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5),
                       labels=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5), name="coeff") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Large pelagics") +
  geom_text(aes(driver, LME, label = sym), color = "black", size = 5) 

d1 <- ggplot(data = mDcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") + labs(x="Driver") + labs(y="LME") +
  scale_fill_distiller(palette = "RdBu", limit = c(-2.2,2.2), 
                       breaks=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5),
                       labels=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5), name="coeff") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Demersals") +
  geom_text(aes(driver, LME, label = sym), color = "black", size = 5) 

a1 <- ggplot(data = mAcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") + labs(x="Driver") + labs(y="") +
  scale_fill_distiller(palette = "RdBu", limit = c(-2.2,2.2), 
                       breaks=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5),
                       labels=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5), name="coeff") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("All fish") +
  geom_text(aes(driver, LME, label = sym), color = "black", size = 5) 

b1 <- ggplot(data = mBcoef, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") + labs(x="Driver") + labs(y="") +
  scale_fill_distiller(palette = "RdBu", limit = c(-2.2,2.2), 
                       breaks=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5),
                       labels=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5), name="coeff") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Benthos") +
  geom_text(aes(driver, LME, label = sym), color = "black", size = 5) 

png(paste0(figp,'Heatmaps_USLME_mlr_coeffs_nu_ALLdiv2SD.png'), 
    width = 6*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
ggarrange(f1,p1,
          d1,a1,
          nrow = 2, ncol = 2,labels = "auto",
          common.legend = TRUE, legend = "right")
dev.off()



# Cluster LMEs
library(rpart)
library(tree)
library(pvclust)
library(mclust)
library(dendextend)

tvec <- Acoef[,1:8]
rownames(tvec) <- Acoef[,9]
tvec[is.na(tvec)] <- 0
d_fish <- dist(tvec) # method="man" # is a bit better
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
gplots::heatmap.2(as.matrix(tvec), 
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
dfish <- color_branches(dfish, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_Fnu_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
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
dfish <- color_branches(dfish, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_Pnu_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
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
dfish <- color_branches(dfish, k=3) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_USLME_Dnu_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
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

