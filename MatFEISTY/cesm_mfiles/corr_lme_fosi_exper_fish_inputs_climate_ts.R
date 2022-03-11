# Find patterns in climate correlations

rm(list=ls())

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/"

library(Hmisc)
library(ggplot2)
library(gridExtra)
library(corrgram)
library(scatterplot3d)
library(rpart)
library(tree)
library(pvclust)
library(mclust)
library(dendextend)
library(gplots)
library("corrplot")
library(scico)
library(reshape2)

### -------------------------------- Raw 4rt -------------------------------
# load data
#FEISTY
ddir <- "/Volumes/MIP/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/"
sigS <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigS.csv"),sep=",",header = F,stringsAsFactors = T)
sigM <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigM.csv"),sep=",",header = F,stringsAsFactors = T)
sigL <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigL.csv"),sep=",",header = F,stringsAsFactors = T)
sigF <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigF.csv"),sep=",",header = F,stringsAsFactors = T)
sigP <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigP.csv"),sep=",",header = F,stringsAsFactors = T)
sigD <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigD.csv"),sep=",",header = F,stringsAsFactors = T)
sigA <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigA.csv"),sep=",",header = F,stringsAsFactors = T)
sigB <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigB.csv"),sep=",",header = F,stringsAsFactors = T)

names(sigS) <- c("LME","climate","lag","rS","pS")
names(sigM) <- c("LME","climate","lag","rM","pM")
names(sigL) <- c("LME","climate","lag","rL","pL")
names(sigF) <- c("LME","climate","lag","rF","pF")
names(sigP) <- c("LME","climate","lag","rP","pP")
names(sigD) <- c("LME","climate","lag","rD","pD")
names(sigA) <- c("LME","climate","lag","rA","pA")
names(sigB) <- c("LME","climate","lag","rB","pB")

sig1 <- merge(sigS,sigM,by=c("LME","climate","lag"),all=T)
sig2 <- merge(sig1,sigL,by=c("LME","climate","lag"),all=T)
sig3 <- merge(sig2,sigF,by=c("LME","climate","lag"),all=T)
sig4 <- merge(sig3,sigP,by=c("LME","climate","lag"),all=T)
sig5 <- merge(sig4,sigD,by=c("LME","climate","lag"),all=T)
sig6 <- merge(sig5,sigA,by=c("LME","climate","lag"),all=T)
sig7 <- merge(sig6,sigB,by=c("LME","climate","lag"),all=T)

#Inputs
cdir <- "/Volumes/MIP/GCM_DATA/CESM/FOSI/"
sigM <- read.csv(paste0(cdir,"LME_fosi_inputs_sigLzooC.csv"),sep=",",header = F,stringsAsFactors = T)
sigL <- read.csv(paste0(cdir,"LME_fosi_inputs_sigZooLoss.csv"),sep=",",header = F,stringsAsFactors = T)
sigP <- read.csv(paste0(cdir,"LME_fosi_inputs_sigTP.csv"),sep=",",header = F,stringsAsFactors = T)
sigD <- read.csv(paste0(cdir,"LME_fosi_inputs_sigDet.csv"),sep=",",header = F,stringsAsFactors = T)
sigB <- read.csv(paste0(cdir,"LME_fosi_inputs_sigTB.csv"),sep=",",header = F,stringsAsFactors = T)

names(sigM) <- c("LME","climate","lag","rZ","pZ")
names(sigL) <- c("LME","climate","lag","rLos","pLos")
names(sigP) <- c("LME","climate","lag","rTP","pTP")
names(sigD) <- c("LME","climate","lag","rDet","pDet")
names(sigB) <- c("LME","climate","lag","rTB","pTB")

sig1 <- merge(sig7,sigM,by=c("LME","climate","lag"),all=T)
sig2 <- merge(sig1,sigL,by=c("LME","climate","lag"),all=T)
sig4 <- merge(sig2,sigP,by=c("LME","climate","lag"),all=T)
sig5 <- merge(sig4,sigB,by=c("LME","climate","lag"),all=T)
sigALL <- merge(sig5,sigD,by=c("LME","climate","lag"),all=T)


## Subset by LME
sigCHK <- subset(sigALL, LME=="CHK")
sigEBS <- subset(sigALL, LME=="EBS")
sigGAK <- subset(sigALL, LME=="GAK")
sigCCE <- subset(sigALL, LME=="CCE")
sigHI <- subset(sigALL, LME=="HI")
sigGMX <- subset(sigALL, LME=="GMX")
sigSE <- subset(sigALL, LME=="SE")
sigNE <- subset(sigALL, LME=="NE")

write.table(sigCHK,paste0(ddir,"LME_fosi_inputs_fished_v14_All_fish03_sig_climate_corrs_CHK.csv"),sep=",",row.names=F)
write.table(sigEBS,paste0(ddir,"LME_fosi_inputs_fished_v14_All_fish03_sig_climate_corrs_EBS.csv"),sep=",",row.names=F)
write.table(sigGAK,paste0(ddir,"LME_fosi_inputs_fished_v14_All_fish03_sig_climate_corrs_GAK.csv"),sep=",",row.names=F)
write.table(sigCCE,paste0(ddir,"LME_fosi_inputs_fished_v14_All_fish03_sig_climate_corrs_CCE.csv"),sep=",",row.names=F)
write.table(sigHI,paste0(ddir,"LME_fosi_inputs_fished_v14_All_fish03_sig_climate_corrs_HI.csv"),sep=",",row.names=F)
write.table(sigGMX,paste0(ddir,"LME_fosi_inputs_fished_v14_All_fish03_sig_climate_corrs_GMX.csv"),sep=",",row.names=F)
write.table(sigSE,paste0(ddir,"LME_fosi_inputs_fished_v14_All_fish03_sig_climate_corrs_SE.csv"),sep=",",row.names=F)
write.table(sigNE,paste0(ddir,"LME_fosi_inputs_fished_v14_All_fish03_sig_climate_corrs_NE.csv"),sep=",",row.names=F)

### ------------------------------------------------------------------------ ###
### ----------------------------- Heatmaps --------------------------------- ###
### ------------------------------------------------------------------------ ###
library(cowplot) #plot_grid

### Rearrange =====================================================
## CCE
#Size
cce1NOI <- subset(sigCCE, climate=="NOI")
cce1NOI <- cce1NOI[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(cce1NOI) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
cce1NOI$lag <- as.factor(cce1NOI$lag)
rcce1NOI <- melt(cce1NOI)
names(rcce1NOI) <- c("Lag","Type","Corr")

cce1PDO <- subset(sigCCE, climate=="PDO")
cce1PDO <- cce1PDO[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(cce1PDO) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
cce1PDO$lag <- as.factor(cce1PDO$lag)
rcce1PDO <- melt(cce1PDO)
names(rcce1PDO) <- c("Lag","Type","Corr")

cce1SOI <- subset(sigCCE, climate=="SOI")
cce1SOI <- cce1SOI[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(cce1SOI) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
cce1SOI$lag <- as.factor(cce1SOI$lag)
rcce1SOI <- melt(cce1SOI)
names(rcce1SOI) <- c("Lag","Type","Corr")

#Type
cce2NOI <- subset(sigCCE, climate=="NOI")
cce2NOI <- cce2NOI[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(cce2NOI) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
cce2NOI$lag <- as.factor(cce2NOI$lag)
rcce2NOI <- melt(cce2NOI)
names(rcce2NOI) <- c("Lag","Type","Corr")

cce2PDO <- subset(sigCCE, climate=="PDO")
cce2PDO <- cce2PDO[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(cce2PDO) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
cce2PDO$lag <- as.factor(cce2PDO$lag)
rcce2PDO <- melt(cce2PDO)
names(rcce2PDO) <- c("Lag","Type","Corr")

cce2SOI <- subset(sigCCE, climate=="SOI")
cce2SOI <- cce2SOI[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(cce2SOI) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
cce2SOI$lag <- as.factor(cce2SOI$lag)
rcce2SOI <- melt(cce2SOI)
names(rcce2SOI) <- c("Lag","Type","Corr")


## EBS
#Size
ebs1NOI <- subset(sigEBS, climate=="NOI")
ebs1NOI <- ebs1NOI[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(ebs1NOI) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
ebs1NOI$lag <- as.factor(ebs1NOI$lag)
rebs1NOI <- melt(ebs1NOI)
names(rebs1NOI) <- c("Lag","Type","Corr")

ebs1PDO <- subset(sigEBS, climate=="PDO")
ebs1PDO <- ebs1PDO[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(ebs1PDO) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
ebs1PDO$lag <- as.factor(ebs1PDO$lag)
rebs1PDO <- melt(ebs1PDO)
names(rebs1PDO) <- c("Lag","Type","Corr")

ebs1SOI <- subset(sigEBS, climate=="SOI")
ebs1SOI <- ebs1SOI[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(ebs1SOI) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
ebs1SOI$lag <- as.factor(ebs1SOI$lag)
rebs1SOI <- melt(ebs1SOI)
names(rebs1SOI) <- c("Lag","Type","Corr")

#Type
ebs2NOI <- subset(sigEBS, climate=="NOI")
ebs2NOI <- ebs2NOI[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(ebs2NOI) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
ebs2NOI$lag <- as.factor(ebs2NOI$lag)
rebs2NOI <- melt(ebs2NOI)
names(rebs2NOI) <- c("Lag","Type","Corr")

ebs2PDO <- subset(sigEBS, climate=="PDO")
ebs2PDO <- ebs2PDO[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(ebs2PDO) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
ebs2PDO$lag <- as.factor(ebs2PDO$lag)
rebs2PDO <- melt(ebs2PDO)
names(rebs2PDO) <- c("Lag","Type","Corr")

ebs2SOI <- subset(sigEBS, climate=="SOI")
ebs2SOI <- ebs2SOI[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(ebs2SOI) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
ebs2SOI$lag <- as.factor(ebs2SOI$lag)
rebs2SOI <- melt(ebs2SOI)
names(rebs2SOI) <- c("Lag","Type","Corr")


## GAK
#Size
gak1MEI <- subset(sigGAK, climate=="MEI")
gak1MEI <- gak1MEI[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(gak1MEI) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
gak1MEI$lag <- as.factor(gak1MEI$lag)
rgak1MEI <- melt(gak1MEI)
names(rgak1MEI) <- c("Lag","Type","Corr")

gak1NOI <- subset(sigGAK, climate=="NOI")
gak1NOI <- gak1NOI[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(gak1NOI) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
gak1NOI$lag <- as.factor(gak1NOI$lag)
rgak1NOI <- melt(gak1NOI)
names(rgak1NOI) <- c("Lag","Type","Corr")

gak1PDO <- subset(sigGAK, climate=="PDO")
gak1PDO <- gak1PDO[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(gak1PDO) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
gak1PDO$lag <- as.factor(gak1PDO$lag)
rgak1PDO <- melt(gak1PDO)
names(rgak1PDO) <- c("Lag","Type","Corr")

#Type
gak2MEI <- subset(sigGAK, climate=="MEI")
gak2MEI <- gak2MEI[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(gak2MEI) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
gak2MEI$lag <- as.factor(gak2MEI$lag)
rgak2MEI <- melt(gak2MEI)
names(rgak2MEI) <- c("Lag","Type","Corr")

gak2NOI <- subset(sigGAK, climate=="NOI")
gak2NOI <- gak2NOI[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(gak2NOI) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
gak2NOI$lag <- as.factor(gak2NOI$lag)
rgak2NOI <- melt(gak2NOI)
names(rgak2NOI) <- c("Lag","Type","Corr")

gak2PDO <- subset(sigGAK, climate=="PDO")
gak2PDO <- gak2PDO[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(gak2PDO) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
gak2PDO$lag <- as.factor(gak2PDO$lag)
rgak2PDO <- melt(gak2PDO)
names(rgak2PDO) <- c("Lag","Type","Corr")


## HI
#Size
hi1NOI <- subset(sigHI, climate=="NOI")
hi1NOI <- hi1NOI[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(hi1NOI) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
hi1NOI$lag <- as.factor(hi1NOI$lag)
rhi1NOI <- melt(hi1NOI)
names(rhi1NOI) <- c("Lag","Type","Corr")

hi1Nino3 <- subset(sigHI, climate=="Nino3")
hi1Nino3 <- hi1Nino3[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(hi1Nino3) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
hi1Nino3$lag <- as.factor(hi1Nino3$lag)
rhi1Nino3 <- melt(hi1Nino3)
names(rhi1Nino3) <- c("Lag","Type","Corr")

hi1Nino4 <- subset(sigHI, climate=="Nino4")
hi1Nino4 <- hi1Nino4[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(hi1Nino4) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
hi1Nino4$lag <- as.factor(hi1Nino4$lag)
rhi1Nino4 <- melt(hi1Nino4)
names(rhi1Nino4) <- c("Lag","Type","Corr")

#Type
hi2NOI <- subset(sigHI, climate=="NOI")
hi2NOI <- hi2NOI[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(hi2NOI) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
hi2NOI$lag <- as.factor(hi2NOI$lag)
rhi2NOI <- melt(hi2NOI)
names(rhi2NOI) <- c("Lag","Type","Corr")

hi2Nino3 <- subset(sigHI, climate=="Nino3")
hi2Nino3 <- hi2Nino3[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(hi2Nino3) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
hi2Nino3$lag <- as.factor(hi2Nino3$lag)
rhi2Nino3 <- melt(hi2Nino3)
names(rhi2Nino3) <- c("Lag","Type","Corr")

hi2Nino4 <- subset(sigHI, climate=="Nino4")
hi2Nino4 <- hi2Nino4[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(hi2Nino4) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
hi2Nino4$lag <- as.factor(hi2Nino4$lag)
rhi2Nino4 <- melt(hi2Nino4)
names(rhi2Nino4) <- c("Lag","Type","Corr")


## CHK
#Size
chk1AMO <- subset(sigCHK, climate=="AMO")
chk1AMO <- chk1AMO[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(chk1AMO) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
chk1AMO$lag <- as.factor(chk1AMO$lag)
rchk1AMO <- melt(chk1AMO)
names(rchk1AMO) <- c("Lag","Type","Corr")

chk1PDO <- subset(sigCHK, climate=="PDO")
chk1PDO <- chk1PDO[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(chk1PDO) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
chk1PDO$lag <- as.factor(chk1PDO$lag)
rchk1PDO <- melt(chk1PDO)
names(rchk1PDO) <- c("Lag","Type","Corr")

chk1MEI <- subset(sigCHK, climate=="MEI")
chk1MEI <- chk1MEI[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(chk1MEI) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
chk1MEI$lag <- as.factor(chk1MEI$lag)
rchk1MEI <- melt(chk1MEI)
names(rchk1MEI) <- c("Lag","Type","Corr")

#Type
chk2AMO <- subset(sigCHK, climate=="AMO")
chk2AMO <- chk2AMO[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(chk2AMO) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
chk2AMO$lag <- as.factor(chk2AMO$lag)
rchk2AMO <- melt(chk2AMO)
names(rchk2AMO) <- c("Lag","Type","Corr")

chk2PDO <- subset(sigCHK, climate=="PDO")
chk2PDO <- chk2PDO[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(chk2PDO) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
chk2PDO$lag <- as.factor(chk2PDO$lag)
rchk2PDO <- melt(chk2PDO)
names(rchk2PDO) <- c("Lag","Type","Corr")

chk2MEI <- subset(sigCHK, climate=="MEI")
chk2MEI <- chk2MEI[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(chk2MEI) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
chk2MEI$lag <- as.factor(chk2MEI$lag)
rchk2MEI <- melt(chk2MEI)
names(rchk2MEI) <- c("Lag","Type","Corr")


##GMX
#Size
gmx1NOI <- subset(sigGMX, climate=="NOI")
gmx1NOI <- gmx1NOI[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(gmx1NOI) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
gmx1NOI$lag <- as.factor(gmx1NOI$lag)
rgmx1NOI <- melt(gmx1NOI)
names(rgmx1NOI) <- c("Lag","Type","Corr")

gmx1PDO <- subset(sigGMX, climate=="PDO")
gmx1PDO <- gmx1PDO[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(gmx1PDO) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
gmx1PDO$lag <- as.factor(gmx1PDO$lag)
rgmx1PDO <- melt(gmx1PDO)
names(rgmx1PDO) <- c("Lag","Type","Corr")

gmx1AMO <- subset(sigGMX, climate=="AMO")
gmx1AMO <- gmx1AMO[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(gmx1AMO) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
gmx1AMO$lag <- as.factor(gmx1AMO$lag)
rgmx1AMO <- melt(gmx1AMO)
names(rgmx1AMO) <- c("Lag","Type","Corr")

#Type
gmx2NOI <- subset(sigGMX, climate=="NOI")
gmx2NOI <- gmx2NOI[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(gmx2NOI) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
gmx2NOI$lag <- as.factor(gmx2NOI$lag)
rgmx2NOI <- melt(gmx2NOI)
names(rgmx2NOI) <- c("Lag","Type","Corr")

gmx2PDO <- subset(sigGMX, climate=="PDO")
gmx2PDO <- gmx2PDO[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(gmx2PDO) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
gmx2PDO$lag <- as.factor(gmx2PDO$lag)
rgmx2PDO <- melt(gmx2PDO)
names(rgmx2PDO) <- c("Lag","Type","Corr")

gmx2AMO <- subset(sigGMX, climate=="AMO")
gmx2AMO <- gmx2AMO[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(gmx2AMO) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
gmx2AMO$lag <- as.factor(gmx2AMO$lag)
rgmx2AMO <- melt(gmx2AMO)
names(rgmx2AMO) <- c("Lag","Type","Corr")


## NE
#Size
ne1NAO <- subset(sigNE, climate=="NAO")
ne1NAO <- ne1NAO[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(ne1NAO) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
ne1NAO$lag <- as.factor(ne1NAO$lag)
rne1NAO <- melt(ne1NAO)
names(rne1NAO) <- c("Lag","Type","Corr")

ne1AMO <- subset(sigNE, climate=="AMO")
ne1AMO <- ne1AMO[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(ne1AMO) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
ne1AMO$lag <- as.factor(ne1AMO$lag)
rne1AMO <- melt(ne1AMO)
names(rne1AMO) <- c("Lag","Type","Corr")

#Type
ne2NAO <- subset(sigNE, climate=="NAO")
ne2NAO <- ne2NAO[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(ne2NAO) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
ne2NAO$lag <- as.factor(ne2NAO$lag)
rne2NAO <- melt(ne2NAO)
names(rne2NAO) <- c("Lag","Type","Corr")

ne2AMO <- subset(sigNE, climate=="AMO")
ne2AMO <- ne2AMO[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(ne2AMO) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
ne2AMO$lag <- as.factor(ne2AMO$lag)
rne2AMO <- melt(ne2AMO)
names(rne2AMO) <- c("Lag","Type","Corr")


##SE
#Size
se1NAO <- subset(sigSE, climate=="NAO")
se1NAO <- se1NAO[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(se1NAO) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
se1NAO$lag <- as.factor(se1NAO$lag)
rse1NAO <- melt(se1NAO)
names(rse1NAO) <- c("Lag","Type","Corr")

se1PDO <- subset(sigSE, climate=="PDO")
se1PDO <- se1PDO[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(se1PDO) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
se1PDO$lag <- as.factor(se1PDO$lag)
rse1PDO <- melt(se1PDO)
names(rse1PDO) <- c("Lag","Type","Corr")

se1AMO <- subset(sigSE, climate=="AMO")
se1AMO <- se1AMO[,c("lag","rTP","rZ","rLos","rS","rM","rL","rA","rTB","rDet")]
names(se1AMO) <- c("lag","TP","LZ","LZm","S","M","L","A","TB","Det")
se1AMO$lag <- as.factor(se1AMO$lag)
rse1AMO <- melt(se1AMO)
names(rse1AMO) <- c("Lag","Type","Corr")

#Type
se2NAO <- subset(sigSE, climate=="NAO")
se2NAO <- se2NAO[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(se2NAO) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
se2NAO$lag <- as.factor(se2NAO$lag)
rse2NAO <- melt(se2NAO)
names(rse2NAO) <- c("Lag","Type","Corr")

se2PDO <- subset(sigSE, climate=="PDO")
se2PDO <- se2PDO[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(se2PDO) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
se2PDO$lag <- as.factor(se2PDO$lag)
rse2PDO <- melt(se2PDO)
names(rse2PDO) <- c("Lag","Type","Corr")

se2AMO <- subset(sigSE, climate=="AMO")
se2AMO <- se2AMO[,c("lag","rTP","rZ","rLos","rF","rP","rTB","rDet","rB","rD")]
names(se2AMO) <- c("lag","TP","LZ","LZm","F","P","TB","Det","B","D")
se2AMO$lag <- as.factor(se2AMO$lag)
rse2AMO <- melt(se2AMO)
names(rse2AMO) <- c("Lag","Type","Corr")




###====== Plots ===================================================================
## CCE
c0 <- ggplot(data = rcce1NOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

c1 <- ggplot(data = rcce1PDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

c2 <- ggplot(data = rcce1SOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="SOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

c3 <- ggplot(data = rcce2NOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

c4 <- ggplot(data = rcce2PDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

c5 <- ggplot(data = rcce2SOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="SOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 


png(paste0(figp,'Heatmaps_CCE_corr_clims_fish_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c0,c1,c2,c3,c4,c5,
           nrow = 2, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## EBS
e0 <- ggplot(data = rebs1NOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e1 <- ggplot(data = rebs1PDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e2 <- ggplot(data = rebs1SOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="SOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e3 <- ggplot(data = rebs2NOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e4 <- ggplot(data = rebs2PDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e5 <- ggplot(data = rebs2SOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="SOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 


png(paste0(figp,'Heatmaps_EBS_corr_clims_fish_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( e0,e1,e2,e3,e4,e5,
           nrow = 2, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## GAK
g0 <- ggplot(data = rgak1NOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g1 <- ggplot(data = rgak1PDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g2 <- ggplot(data = rgak1MEI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="MEI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g3 <- ggplot(data = rgak2NOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g4 <- ggplot(data = rgak2PDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g5 <- ggplot(data = rgak2MEI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="MEI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 


png(paste0(figp,'Heatmaps_GAK_corr_clims_fish_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( g0,g1,g2,g3,g4,g5,
           nrow = 2, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## HI
h0 <- ggplot(data = rhi1Nino3, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino3\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

h1 <- ggplot(data = rhi1Nino4, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino4\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

h2 <- ggplot(data = rhi1NOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

h3 <- ggplot(data = rhi2Nino3, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino3\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

h4 <- ggplot(data = rhi2Nino4, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino4\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

h5 <- ggplot(data = rhi2NOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_HI_corr_clims_fish_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( h0,h1,h2,h3,h4,h5,
           nrow = 2, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## CHK
k0 <- ggplot(data = rchk1MEI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="MEI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

k1 <- ggplot(data = rchk1PDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

k2 <- ggplot(data = rchk1AMO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AMO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3)

k3 <- ggplot(data = rchk2MEI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="MEI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

k4 <- ggplot(data = rchk2PDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

k5 <- ggplot(data = rchk2AMO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AMO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3)

png(paste0(figp,'Heatmaps_CHK_corr_clims_fish_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( k0,k1,k2,k3,k4,k5,
           nrow = 2, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## GMX
m0 <- ggplot(data = rgmx1AMO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AMO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

m1 <- ggplot(data = rgmx1PDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

m2 <- ggplot(data = rgmx1NOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

m3 <- ggplot(data = rgmx2AMO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AMO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

m4 <- ggplot(data = rgmx2PDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

m5 <- ggplot(data = rgmx2NOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_GMX_corr_clims_fish_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( m0,m1,m2,m3,m4,m5,
           nrow = 2, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## SE
s0 <- ggplot(data = rse1AMO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AMO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

s1 <- ggplot(data = rse1NAO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NAO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

s2 <- ggplot(data = rse1PDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

s3 <- ggplot(data = rse2AMO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AMO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

s4 <- ggplot(data = rse2NAO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NAO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

s5 <- ggplot(data = rse2PDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_SE_corr_clims_fish_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( s0,s1,s2,s3,s4,s5,
           nrow = 2, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## NE
n0 <- ggplot(data = rne1AMO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AMO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

n1 <- ggplot(data = rne1NAO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NAO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

n2 <- ggplot(data = rne2AMO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AMO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

n3 <- ggplot(data = rne2NAO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NAO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 


png(paste0(figp,'Heatmaps_NE_corr_clims_fish_inputs.png'), 
    width = 6.6*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( n0,n1,n2,n3,
           nrow = 2, ncol = 2,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()

