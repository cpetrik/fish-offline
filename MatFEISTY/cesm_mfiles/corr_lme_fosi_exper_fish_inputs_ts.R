# Find patterns in correlations between fish and inputs

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

### --------------------------------------------------------------
# load data
#FEISTY
ddir <- "/Volumes/MIP/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/"
sigS <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigS_inputs.csv"),sep=",",header = F,stringsAsFactors = T)
sigM <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigM_inputs.csv"),sep=",",header = F,stringsAsFactors = T)
sigL <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigL_inputs.csv"),sep=",",header = F,stringsAsFactors = T)
sigF <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigF_inputs.csv"),sep=",",header = F,stringsAsFactors = T)
sigP <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigP_inputs.csv"),sep=",",header = F,stringsAsFactors = T)
sigD <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigD_inputs.csv"),sep=",",header = F,stringsAsFactors = T)
sigA <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigA_inputs.csv"),sep=",",header = F,stringsAsFactors = T)
sigB <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigB_inputs.csv"),sep=",",header = F,stringsAsFactors = T)

names(sigS) <- c("LME","input","lag","rS","pS")
names(sigM) <- c("LME","input","lag","rM","pM")
names(sigL) <- c("LME","input","lag","rL","pL")
names(sigF) <- c("LME","input","lag","rF","pF")
names(sigP) <- c("LME","input","lag","rP","pP")
names(sigD) <- c("LME","input","lag","rD","pD")
names(sigA) <- c("LME","input","lag","rA","pA")
names(sigB) <- c("LME","input","lag","rB","pB")

sig1 <- merge(sigS,sigM,by=c("LME","input","lag"),all=T)
sig2 <- merge(sig1,sigL,by=c("LME","input","lag"),all=T)
sig3 <- merge(sig2,sigF,by=c("LME","input","lag"),all=T)
sig4 <- merge(sig3,sigP,by=c("LME","input","lag"),all=T)
sig5 <- merge(sig4,sigD,by=c("LME","input","lag"),all=T)
sig6 <- merge(sig5,sigA,by=c("LME","input","lag"),all=T)
sigALL <- merge(sig6,sigB,by=c("LME","input","lag"),all=T)


## Subset by LME
sigCHK <- subset(sigALL, LME=="CHK")
sigEBS <- subset(sigALL, LME=="EBS")
sigGAK <- subset(sigALL, LME=="GAK")
sigCCE <- subset(sigALL, LME=="CCE")
sigHI <- subset(sigALL, LME=="HI")
sigGMX <- subset(sigALL, LME=="GMX")
sigSE <- subset(sigALL, LME=="SE")
sigNE <- subset(sigALL, LME=="NE")

write.table(sigCHK,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_input_corrs_CHK.csv"),sep=",",row.names=F)
write.table(sigEBS,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_input_corrs_EBS.csv"),sep=",",row.names=F)
write.table(sigGAK,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_input_corrs_GAK.csv"),sep=",",row.names=F)
write.table(sigCCE,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_input_corrs_CCE.csv"),sep=",",row.names=F)
write.table(sigHI,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_input_corrs_HI.csv"),sep=",",row.names=F)
write.table(sigGMX,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_input_corrs_GMX.csv"),sep=",",row.names=F)
write.table(sigSE,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_input_corrs_SE.csv"),sep=",",row.names=F)
write.table(sigNE,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_input_corrs_NE.csv"),sep=",",row.names=F)

### ------------------------------------------------------------------------ ###
### ----------------------------- Heatmaps --------------------------------- ###
### ------------------------------------------------------------------------ ###
library(cowplot) #plot_grid

### Rearrange =====================================================
## CCE
cceTp <- subset(sigCCE, input=="Tp")
cceTp <- cceTp[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(cceTp) <- c("lag","B","S","M","L","F","P","D","A")
cceTp$lag <- as.factor(cceTp$lag)
rcceTp <- melt(cceTp)
names(rcceTp) <- c("Lag","Type","Corr")

cceTb <- subset(sigCCE, input=="Tb")
cceTb <- cceTb[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(cceTb) <- c("lag","B","S","M","L","F","P","D","A")
cceTb$lag <- as.factor(cceTb$lag)
rcceTb <- melt(cceTb)
names(rcceTb) <- c("Lag","Type","Corr")

cceDet <- subset(sigCCE, input=="Det")
cceDet <- cceDet[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(cceDet) <- c("lag","B","S","M","L","F","P","D","A")
cceDet$lag <- as.factor(cceDet$lag)
rcceDet <- melt(cceDet)
names(rcceDet) <- c("Lag","Type","Corr")

cceLZbiom <- subset(sigCCE, input=="LZbiom")
cceLZbiom <- cceLZbiom[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(cceLZbiom) <- c("lag","B","S","M","L","F","P","D","A")
cceLZbiom$lag <- as.factor(cceLZbiom$lag)
rcceLZbiom <- melt(cceLZbiom)
names(rcceLZbiom) <- c("Lag","Type","Corr")

cceLZloss <- subset(sigCCE, input=="LZloss")
cceLZloss <- cceLZloss[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(cceLZloss) <- c("lag","B","S","M","L","F","P","D","A")
cceLZloss$lag <- as.factor(cceLZloss$lag)
rcceLZloss <- melt(cceLZloss)
names(rcceLZloss) <- c("Lag","Type","Corr")


## EBS
ebsTp <- subset(sigEBS, input=="Tp")
ebsTp <- ebsTp[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(ebsTp) <- c("lag","B","S","M","L","F","P","D","A")
ebsTp$lag <- as.factor(ebsTp$lag)
rebsTp <- melt(ebsTp)
names(rebsTp) <- c("Lag","Type","Corr")

ebsTb <- subset(sigEBS, input=="Tb")
ebsTb <- ebsTb[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(ebsTb) <- c("lag","B","S","M","L","F","P","D","A")
ebsTb$lag <- as.factor(ebsTb$lag)
rebsTb <- melt(ebsTb)
names(rebsTb) <- c("Lag","Type","Corr")

ebsDet <- subset(sigEBS, input=="Det")
ebsDet <- ebsDet[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(ebsDet) <- c("lag","B","S","M","L","F","P","D","A")
ebsDet$lag <- as.factor(ebsDet$lag)
rebsDet <- melt(ebsDet)
names(rebsDet) <- c("Lag","Type","Corr")

ebsLZbiom <- subset(sigEBS, input=="LZbiom")
ebsLZbiom <- ebsLZbiom[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(ebsLZbiom) <- c("lag","B","S","M","L","F","P","D","A")
ebsLZbiom$lag <- as.factor(ebsLZbiom$lag)
rebsLZbiom <- melt(ebsLZbiom)
names(rebsLZbiom) <- c("Lag","Type","Corr")

ebsLZloss <- subset(sigEBS, input=="LZloss")
ebsLZloss <- ebsLZloss[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(ebsLZloss) <- c("lag","B","S","M","L","F","P","D","A")
ebsLZloss$lag <- as.factor(ebsLZloss$lag)
rebsLZloss <- melt(ebsLZloss)
names(rebsLZloss) <- c("Lag","Type","Corr")


## GAK
gakTp <- subset(sigGAK, input=="Tp")
gakTp <- gakTp[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gakTp) <- c("lag","B","S","M","L","F","P","D","A")
gakTp$lag <- as.factor(gakTp$lag)
rgakTp <- melt(gakTp)
names(rgakTp) <- c("Lag","Type","Corr")

gakTb <- subset(sigGAK, input=="Tb")
gakTb <- gakTb[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gakTb) <- c("lag","B","S","M","L","F","P","D","A")
gakTb$lag <- as.factor(gakTb$lag)
rgakTb <- melt(gakTb)
names(rgakTb) <- c("Lag","Type","Corr")

gakDet <- subset(sigGAK, input=="Det")
gakDet <- gakDet[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gakDet) <- c("lag","B","S","M","L","F","P","D","A")
gakDet$lag <- as.factor(gakDet$lag)
rgakDet <- melt(gakDet)
names(rgakDet) <- c("Lag","Type","Corr")

gakLZbiom <- subset(sigGAK, input=="LZbiom")
gakLZbiom <- gakLZbiom[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gakLZbiom) <- c("lag","B","S","M","L","F","P","D","A")
gakLZbiom$lag <- as.factor(gakLZbiom$lag)
rgakLZbiom <- melt(gakLZbiom)
names(rgakLZbiom) <- c("Lag","Type","Corr")

gakLZloss <- subset(sigGAK, input=="LZloss")
gakLZloss <- gakLZloss[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gakLZloss) <- c("lag","B","S","M","L","F","P","D","A")
gakLZloss$lag <- as.factor(gakLZloss$lag)
rgakLZloss <- melt(gakLZloss)
names(rgakLZloss) <- c("Lag","Type","Corr")


## HI
hiTp <- subset(sigHI, input=="Tp")
hiTp <- hiTp[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(hiTp) <- c("lag","B","S","M","L","F","P","D","A")
hiTp$lag <- as.factor(hiTp$lag)
rhiTp <- melt(hiTp)
names(rhiTp) <- c("Lag","Type","Corr")

hiTb <- subset(sigHI, input=="Tb")
hiTb <- hiTb[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(hiTb) <- c("lag","B","S","M","L","F","P","D","A")
hiTb$lag <- as.factor(hiTb$lag)
rhiTb <- melt(hiTb)
names(rhiTb) <- c("Lag","Type","Corr")

hiDet <- subset(sigHI, input=="Det")
hiDet <- hiDet[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(hiDet) <- c("lag","B","S","M","L","F","P","D","A")
hiDet$lag <- as.factor(hiDet$lag)
rhiDet <- melt(hiDet)
names(rhiDet) <- c("Lag","Type","Corr")

hiLZbiom <- subset(sigHI, input=="LZbiom")
hiLZbiom <- hiLZbiom[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(hiLZbiom) <- c("lag","B","S","M","L","F","P","D","A")
hiLZbiom$lag <- as.factor(hiLZbiom$lag)
rhiLZbiom <- melt(hiLZbiom)
names(rhiLZbiom) <- c("Lag","Type","Corr")

hiLZloss <- subset(sigHI, input=="LZloss")
hiLZloss <- hiLZloss[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(hiLZloss) <- c("lag","B","S","M","L","F","P","D","A")
hiLZloss$lag <- as.factor(hiLZloss$lag)
rhiLZloss <- melt(hiLZloss)
names(rhiLZloss) <- c("Lag","Type","Corr")


## CHK
chkTp <- subset(sigCHK, input=="Tp")
chkTp <- chkTp[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(chkTp) <- c("lag","B","S","M","L","F","P","D","A")
chkTp$lag <- as.factor(chkTp$lag)
rchkTp <- melt(chkTp)
names(rchkTp) <- c("Lag","Type","Corr")

chkTb <- subset(sigCHK, input=="Tb")
chkTb <- chkTb[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(chkTb) <- c("lag","B","S","M","L","F","P","D","A")
chkTb$lag <- as.factor(chkTb$lag)
rchkTb <- melt(chkTb)
names(rchkTb) <- c("Lag","Type","Corr")

chkDet <- subset(sigCHK, input=="Det")
chkDet <- chkDet[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(chkDet) <- c("lag","B","S","M","L","F","P","D","A")
chkDet$lag <- as.factor(chkDet$lag)
rchkDet <- melt(chkDet)
names(rchkDet) <- c("Lag","Type","Corr")

chkLZbiom <- subset(sigCHK, input=="LZbiom")
chkLZbiom <- chkLZbiom[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(chkLZbiom) <- c("lag","B","S","M","L","F","P","D","A")
chkLZbiom$lag <- as.factor(chkLZbiom$lag)
rchkLZbiom <- melt(chkLZbiom)
names(rchkLZbiom) <- c("Lag","Type","Corr")

chkLZloss <- subset(sigCHK, input=="LZloss")
chkLZloss <- chkLZloss[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(chkLZloss) <- c("lag","B","S","M","L","F","P","D","A")
chkLZloss$lag <- as.factor(chkLZloss$lag)
rchkLZloss <- melt(chkLZloss)
names(rchkLZloss) <- c("Lag","Type","Corr")


##GMX
gmxTp <- subset(sigGMX, input=="Tp")
gmxTp <- gmxTp[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gmxTp) <- c("lag","B","S","M","L","F","P","D","A")
gmxTp$lag <- as.factor(gmxTp$lag)
rgmxTp <- melt(gmxTp)
names(rgmxTp) <- c("Lag","Type","Corr")

gmxTb <- subset(sigGMX, input=="Tb")
gmxTb <- gmxTb[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gmxTb) <- c("lag","B","S","M","L","F","P","D","A")
gmxTb$lag <- as.factor(gmxTb$lag)
rgmxTb <- melt(gmxTb)
names(rgmxTb) <- c("Lag","Type","Corr")

gmxDet <- subset(sigGMX, input=="Det")
gmxDet <- gmxDet[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gmxDet) <- c("lag","B","S","M","L","F","P","D","A")
gmxDet$lag <- as.factor(gmxDet$lag)
rgmxDet <- melt(gmxDet)
names(rgmxDet) <- c("Lag","Type","Corr")

gmxLZbiom <- subset(sigGMX, input=="LZbiom")
gmxLZbiom <- gmxLZbiom[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gmxLZbiom) <- c("lag","B","S","M","L","F","P","D","A")
gmxLZbiom$lag <- as.factor(gmxLZbiom$lag)
rgmxLZbiom <- melt(gmxLZbiom)
names(rgmxLZbiom) <- c("Lag","Type","Corr")

gmxLZloss <- subset(sigGMX, input=="LZloss")
gmxLZloss <- gmxLZloss[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gmxLZloss) <- c("lag","B","S","M","L","F","P","D","A")
gmxLZloss$lag <- as.factor(gmxLZloss$lag)
rgmxLZloss <- melt(gmxLZloss)
names(rgmxLZloss) <- c("Lag","Type","Corr")


## NE
neTp <- subset(sigNE, input=="Tp")
neTp <- neTp[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(neTp) <- c("lag","B","S","M","L","F","P","D","A")
neTp$lag <- as.factor(neTp$lag)
rneTp <- melt(neTp)
names(rneTp) <- c("Lag","Type","Corr")

neTb <- subset(sigNE, input=="Tb")
neTb <- neTb[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(neTb) <- c("lag","B","S","M","L","F","P","D","A")
neTb$lag <- as.factor(neTb$lag)
rneTb <- melt(neTb)
names(rneTb) <- c("Lag","Type","Corr")

neDet <- subset(sigNE, input=="Det")
neDet <- neDet[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(neDet) <- c("lag","B","S","M","L","F","P","D","A")
neDet$lag <- as.factor(neDet$lag)
rneDet <- melt(neDet)
names(rneDet) <- c("Lag","Type","Corr")

neLZbiom <- subset(sigNE, input=="LZbiom")
neLZbiom <- neLZbiom[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(neLZbiom) <- c("lag","B","S","M","L","F","P","D","A")
neLZbiom$lag <- as.factor(neLZbiom$lag)
rneLZbiom <- melt(neLZbiom)
names(rneLZbiom) <- c("Lag","Type","Corr")

neLZloss <- subset(sigNE, input=="LZloss")
neLZloss <- neLZloss[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(neLZloss) <- c("lag","B","S","M","L","F","P","D","A")
neLZloss$lag <- as.factor(neLZloss$lag)
rneLZloss <- melt(neLZloss)
names(rneLZloss) <- c("Lag","Type","Corr")


##SE
seTp <- subset(sigSE, input=="Tp")
seTp <- seTp[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(seTp) <- c("lag","B","S","M","L","F","P","D","A")
seTp$lag <- as.factor(seTp$lag)
rseTp <- melt(seTp)
names(rseTp) <- c("Lag","Type","Corr")

seTb <- subset(sigSE, input=="Tb")
seTb <- seTb[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(seTb) <- c("lag","B","S","M","L","F","P","D","A")
seTb$lag <- as.factor(seTb$lag)
rseTb <- melt(seTb)
names(rseTb) <- c("Lag","Type","Corr")

seDet <- subset(sigSE, input=="Det")
seDet <- seDet[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(seDet) <- c("lag","B","S","M","L","F","P","D","A")
seDet$lag <- as.factor(seDet$lag)
rseDet <- melt(seDet)
names(rseDet) <- c("Lag","Type","Corr")

seLZbiom <- subset(sigSE, input=="LZbiom")
seLZbiom <- seLZbiom[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(seLZbiom) <- c("lag","B","S","M","L","F","P","D","A")
seLZbiom$lag <- as.factor(seLZbiom$lag)
rseLZbiom <- melt(seLZbiom)
names(rseLZbiom) <- c("Lag","Type","Corr")

seLZloss <- subset(sigSE, input=="LZloss")
seLZloss <- seLZloss[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(seLZloss) <- c("lag","B","S","M","L","F","P","D","A")
seLZloss$lag <- as.factor(seLZloss$lag)
rseLZloss <- melt(seLZloss)
names(rseLZloss) <- c("Lag","Type","Corr")


###====== Plots ===================================================================
## CCE
c5 <- ggplot(data = rcceTp, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tp\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

c4 <- ggplot(data = rcceTb, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tb\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

c1 <- ggplot(data = rcceDet, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Det\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

c2 <- ggplot(data = rcceLZbiom, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZbiom\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

c3 <- ggplot(data = rcceLZloss, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZloss\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_CCE_corr_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c1,c2,c3,c4,c5,
           nrow = 2, ncol = 3)
dev.off()


## EBS
e5 <- ggplot(data = rebsTp, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tp\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e4 <- ggplot(data = rebsTb, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tb\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e1 <- ggplot(data = rebsDet, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Det\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e2 <- ggplot(data = rebsLZbiom, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZbiom\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e3 <- ggplot(data = rebsLZloss, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZloss\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_EBS_corr_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( e1,e2,e3,e4,e5,
           nrow = 2, ncol = 3)
           #rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           #align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## GAK
g5 <- ggplot(data = rgakTp, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tp\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g4 <- ggplot(data = rgakTb, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tb\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g1 <- ggplot(data = rgakDet, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Det\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g2 <- ggplot(data = rgakLZbiom, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZbiom\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g3 <- ggplot(data = rgakLZloss, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZloss\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_GAK_corr_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( g1,g2,g3,g4,g5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## HI
h5 <- ggplot(data = rhiTp, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tp\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

h4 <- ggplot(data = rhiTb, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tb\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

h1 <- ggplot(data = rhiDet, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Det\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

h2 <- ggplot(data = rhiLZbiom, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZbiom\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

h3 <- ggplot(data = rhiLZloss, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZloss\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_HI_corr_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( h1,h2,h3,h4,h5,
           nrow = 2, ncol = 3)#,
           #align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## CHK
k5 <- ggplot(data = rchkTp, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tp\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

k4 <- ggplot(data = rchkTb, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tb\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

k1 <- ggplot(data = rchkDet, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Det\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

k2 <- ggplot(data = rchkLZbiom, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZbiom\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

k3 <- ggplot(data = rchkLZloss, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZloss\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_CHK_corr_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( k1,k2,k3,k4,k5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## GMX
m5 <- ggplot(data = rgmxTp, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tp\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

m4 <- ggplot(data = rgmxTb, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tb\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

m1 <- ggplot(data = rgmxDet, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Det\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

m2 <- ggplot(data = rgmxLZbiom, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.99,0.99),  
                       name="LZbiom\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

m3 <- ggplot(data = rgmxLZloss, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZloss\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_GMX_corr_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( m1,m2,m3,m4,m5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## SE
s5 <- ggplot(data = rseTp, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tp\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

s4 <- ggplot(data = rseTb, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tb\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

s1 <- ggplot(data = rseDet, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Det\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

s2 <- ggplot(data = rseLZbiom, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZbiom\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

s3 <- ggplot(data = rseLZloss, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZloss\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_SE_corr_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( s1,s2,s3,s4,s5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## NE
n5 <- ggplot(data = rneTp, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tp\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

n4 <- ggplot(data = rneTb, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tb\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

n1 <- ggplot(data = rneDet, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Det\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

n2 <- ggplot(data = rneLZbiom, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZbiom\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

n3 <- ggplot(data = rneLZloss, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZloss\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_NE_corr_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( n1,n2,n3,n4,n5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()






