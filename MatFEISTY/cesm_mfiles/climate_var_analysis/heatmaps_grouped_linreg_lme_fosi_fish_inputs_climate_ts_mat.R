# Find patterns in correlations between climate with fish and inputs
# Heatmaps of the strongest regressions in regions of interest

rm(list=ls())

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/corrs/"

library(ggplot2)
library(reshape2) #melt
library(cowplot) #plot_grid
library(gridExtra)
library(scatterplot3d)
# library(rpart)
# library(tree)
# library(pvclust)
# library(mclust)
# library(dendextend)
library(gplots)
library("coefplot")
library(scico)

### --------------------------------------------------------------
# load data
#FEISTY
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"

### ------------------------------------------------------------------------ ###
### ----------------------------- Heatmaps --------------------------------- ###
### ------------------------------------------------------------------------ ###

### Rearrange =====================================================
## CCE
ccePDO <- read.csv(paste0(datap,"CCE_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
ccePDO$Lag <- as.factor(ccePDO$Lag)
ccePDO$Type[ccePDO$Type=="Zlos"] <- "ZmLoss"
ccePDO$Type[ccePDO$Type=="Zoo"] <- "Zmeso"
ccePDO$Type <- factor(ccePDO$Type, 
                         levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rccePDO <- ccePDO[,c("Lag","Type","coef","p")]
rccePDO$p[rccePDO$p > 0.1] <- NA
rccePDO$sym <- rccePDO$p
rccePDO$sym[rccePDO$p <= 0.05] <- "*"
rccePDO$sym[rccePDO$p > 0.05] <- "."
rccePDO$sym[rccePDO$p > 0.1] <- NA

cceNino34 <- read.csv(paste0(datap,"CCE_regress_Nino34_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
cceNino34$Lag <- as.factor(cceNino34$Lag)
cceNino34$Type[cceNino34$Type=="Zlos"] <- "ZmLoss"
cceNino34$Type[cceNino34$Type=="Zoo"] <- "Zmeso"
cceNino34$Type <- factor(cceNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rcceNino34 <- cceNino34[,c("Lag","Type","coef","p")]
rcceNino34$p[rcceNino34$p > 0.1] <- NA
rcceNino34$sym <- rcceNino34$p
rcceNino34$sym[rcceNino34$p <= 0.05] <- "*"
rcceNino34$sym[rcceNino34$p > 0.05] <- "."
rcceNino34$sym[rcceNino34$p > 0.1] <- NA


## EBS
ebsPDO <- read.csv(paste0(datap,"EBS_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
ebsPDO$Lag <- as.factor(ebsPDO$Lag)
ebsPDO$Type[ebsPDO$Type=="Zlos"] <- "ZmLoss"
ebsPDO$Type[ebsPDO$Type=="Zoo"] <- "Zmeso"
ebsPDO$Type <- factor(ebsPDO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rebsPDO <- ebsPDO[,c("Lag","Type","coef","p")]
rebsPDO$p[rebsPDO$p > 0.1] <- NA
rebsPDO$sym <- rebsPDO$p
rebsPDO$sym[rebsPDO$p <= 0.05] <- "*"
rebsPDO$sym[rebsPDO$p > 0.05] <- "."
rebsPDO$sym[rebsPDO$p > 0.1] <- NA

ebsNino34 <- read.csv(paste0(datap,"EBS_regress_Nino34_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
ebsNino34$Lag <- as.factor(ebsNino34$Lag)
ebsNino34$Type[ebsNino34$Type=="Zlos"] <- "ZmLoss"
ebsNino34$Type[ebsNino34$Type=="Zoo"] <- "Zmeso"
ebsNino34$Type <- factor(ebsNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rebsNino34 <- ebsNino34[,c("Lag","Type","coef","p")]
rebsNino34$p[rebsNino34$p > 0.1] <- NA
rebsNino34$sym <- rebsNino34$p
rebsNino34$sym[rebsNino34$p <= 0.05] <- "*"
rebsNino34$sym[rebsNino34$p > 0.05] <- "."
rebsNino34$sym[rebsNino34$p > 0.1] <- NA

ebsAO <- read.csv(paste0(datap,"EBS_regress_AO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
ebsAO$Lag <- as.factor(ebsAO$Lag)
ebsAO$Type[ebsAO$Type=="Zlos"] <- "ZmLoss"
ebsAO$Type[ebsAO$Type=="Zoo"] <- "Zmeso"
ebsAO$Type <- factor(ebsAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rebsAO <- ebsAO[,c("Lag","Type","coef","p")]
rebsAO$p[rebsAO$p > 0.1] <- NA
rebsAO$sym <- rebsAO$p
rebsAO$sym[rebsAO$p <= 0.05] <- "*"
rebsAO$sym[rebsAO$p > 0.05] <- "."
rebsAO$sym[rebsAO$p > 0.1] <- NA

## AI
aiPDO <- read.csv(paste0(datap,"AI_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
aiPDO$Lag <- as.factor(aiPDO$Lag)
aiPDO$Type[aiPDO$Type=="Zlos"] <- "ZmLoss"
aiPDO$Type[aiPDO$Type=="Zoo"] <- "Zmeso"
aiPDO$Type <- factor(aiPDO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
raiPDO <- aiPDO[,c("Lag","Type","coef","p")]
raiPDO$p[raiPDO$p > 0.1] <- NA
raiPDO$sym <- raiPDO$p
raiPDO$sym[raiPDO$p <= 0.05] <- "*"
raiPDO$sym[raiPDO$p > 0.05] <- "."
raiPDO$sym[raiPDO$p > 0.1] <- NA

aiNino34 <- read.csv(paste0(datap,"AI_regress_Nino34_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
aiNino34$Lag <- as.factor(aiNino34$Lag)
aiNino34$Type[aiNino34$Type=="Zlos"] <- "ZmLoss"
aiNino34$Type[aiNino34$Type=="Zoo"] <- "Zmeso"
aiNino34$Type <- factor(aiNino34$Type, 
                         levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
raiNino34 <- aiNino34[,c("Lag","Type","coef","p")]
raiNino34$p[raiNino34$p > 0.1] <- NA
raiNino34$sym <- raiNino34$p
raiNino34$sym[raiNino34$p <= 0.05] <- "*"
raiNino34$sym[raiNino34$p > 0.05] <- "."
raiNino34$sym[raiNino34$p > 0.1] <- NA

aiAO <- read.csv(paste0(datap,"AI_regress_AO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
aiAO$Lag <- as.factor(aiAO$Lag)
aiAO$Type[aiAO$Type=="Zlos"] <- "ZmLoss"
aiAO$Type[aiAO$Type=="Zoo"] <- "Zmeso"
aiAO$Type <- factor(aiAO$Type, 
                     levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
raiAO <- aiAO[,c("Lag","Type","coef","p")]
raiAO$p[raiAO$p > 0.1] <- NA
raiAO$sym <- raiAO$p
raiAO$sym[raiAO$p <= 0.05] <- "*"
raiAO$sym[raiAO$p > 0.05] <- "."
raiAO$sym[raiAO$p > 0.1] <- NA



## GAK
gakPDO <- read.csv(paste0(datap,"GAK_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gakPDO$Lag <- as.factor(gakPDO$Lag)
gakPDO$Type[gakPDO$Type=="Zlos"] <- "ZmLoss"
gakPDO$Type[gakPDO$Type=="Zoo"] <- "Zmeso"
gakPDO$Type <- factor(gakPDO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rgakPDO <- gakPDO[,c("Lag","Type","coef","p")]
rgakPDO$p[rgakPDO$p > 0.1] <- NA
rgakPDO$sym <- rgakPDO$p
rgakPDO$sym[rgakPDO$p <= 0.05] <- "*"
rgakPDO$sym[rgakPDO$p > 0.05] <- "."
rgakPDO$sym[rgakPDO$p > 0.1] <- NA

gakNino34 <- read.csv(paste0(datap,"GAK_regress_Nino34_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gakNino34$Lag <- as.factor(gakNino34$Lag)
gakNino34$Type[gakNino34$Type=="Zlos"] <- "ZmLoss"
gakNino34$Type[gakNino34$Type=="Zoo"] <- "Zmeso"
gakNino34$Type <- factor(gakNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rgakNino34 <- gakNino34[,c("Lag","Type","coef","p")]
rgakNino34$p[rgakNino34$p > 0.1] <- NA
rgakNino34$sym <- rgakNino34$p
rgakNino34$sym[rgakNino34$p <= 0.05] <- "*"
rgakNino34$sym[rgakNino34$p > 0.05] <- "."
rgakNino34$sym[rgakNino34$p > 0.1] <- NA

gakAO <- read.csv(paste0(datap,"GAK_regress_AO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gakAO$Lag <- as.factor(gakAO$Lag)
gakAO$Type[gakAO$Type=="Zlos"] <- "ZmLoss"
gakAO$Type[gakAO$Type=="Zoo"] <- "Zmeso"
gakAO$Type <- factor(gakAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rgakAO <- gakAO[,c("Lag","Type","coef","p")]
rgakAO$p[rgakAO$p > 0.1] <- NA
rgakAO$sym <- rgakAO$p
rgakAO$sym[rgakAO$p <= 0.05] <- "*"
rgakAO$sym[rgakAO$p > 0.05] <- "."
rgakAO$sym[rgakAO$p > 0.1] <- NA


## HI
hiPDO <- read.csv(paste0(datap,"HI_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
hiPDO$Lag <- as.factor(hiPDO$Lag)
hiPDO$Type[hiPDO$Type=="Zlos"] <- "ZmLoss"
hiPDO$Type[hiPDO$Type=="Zoo"] <- "Zmeso"
hiPDO$Type <- factor(hiPDO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rhiPDO <- hiPDO[,c("Lag","Type","coef","p")]
rhiPDO$p[rhiPDO$p > 0.1] <- NA
rhiPDO$sym <- rhiPDO$p
rhiPDO$sym[rhiPDO$p <= 0.05] <- "*"
rhiPDO$sym[rhiPDO$p > 0.05] <- "."
rhiPDO$sym[rhiPDO$p > 0.1] <- NA

hiNino34 <- read.csv(paste0(datap,"HI_regress_Nino34_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
hiNino34$Lag <- as.factor(hiNino34$Lag)
hiNino34$Type[hiNino34$Type=="Zlos"] <- "ZmLoss"
hiNino34$Type[hiNino34$Type=="Zoo"] <- "Zmeso"
hiNino34$Type <- factor(hiNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rhiNino34 <- hiNino34[,c("Lag","Type","coef","p")]
rhiNino34$p[rhiNino34$p > 0.1] <- NA
rhiNino34$sym <- rhiNino34$p
rhiNino34$sym[rhiNino34$p <= 0.05] <- "*"
rhiNino34$sym[rhiNino34$p > 0.05] <- "."
rhiNino34$sym[rhiNino34$p > 0.1] <- NA


## CHK
chkPDO <- read.csv(paste0(datap,"CHK_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
chkPDO$Lag <- as.factor(chkPDO$Lag)
chkPDO$Type[chkPDO$Type=="Zlos"] <- "ZmLoss"
chkPDO$Type[chkPDO$Type=="Zoo"] <- "Zmeso"
chkPDO$Type <- factor(chkPDO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rchkPDO <- chkPDO[,c("Lag","Type","coef","p")]
rchkPDO$p[rchkPDO$p > 0.1] <- NA
rchkPDO$sym <- rchkPDO$p
rchkPDO$sym[rchkPDO$p <= 0.05] <- "*"
rchkPDO$sym[rchkPDO$p > 0.05] <- "."
rchkPDO$sym[rchkPDO$p > 0.1] <- NA

chkAMO <- read.csv(paste0(datap,"CHK_regress_AMO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
chkAMO$Lag <- as.factor(chkAMO$Lag)
chkAMO$Type[chkAMO$Type=="Zlos"] <- "ZmLoss"
chkAMO$Type[chkAMO$Type=="Zoo"] <- "Zmeso"
chkAMO$Type <- factor(chkAMO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rchkAMO <- chkAMO[,c("Lag","Type","coef","p")]
rchkAMO$p[rchkAMO$p > 0.1] <- NA
rchkAMO$sym <- rchkAMO$p
rchkAMO$sym[rchkAMO$p <= 0.05] <- "*"
rchkAMO$sym[rchkAMO$p > 0.05] <- "."
rchkAMO$sym[rchkAMO$p > 0.1] <- NA

chkAO <- read.csv(paste0(datap,"CHK_regress_AO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
chkAO$Lag <- as.factor(chkAO$Lag)
chkAO$Type[chkAO$Type=="Zlos"] <- "ZmLoss"
chkAO$Type[chkAO$Type=="Zoo"] <- "Zmeso"
chkAO$Type <- factor(chkAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rchkAO <- chkAO[,c("Lag","Type","coef","p")]
rchkAO$p[rchkAO$p > 0.1] <- NA
rchkAO$sym <- rchkAO$p
rchkAO$sym[rchkAO$p <= 0.05] <- "*"
rchkAO$sym[rchkAO$p > 0.05] <- "."
rchkAO$sym[rchkAO$p > 0.1] <- NA


##GMX
gmxPDO <- read.csv(paste0(datap,"GMX_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gmxPDO$Lag <- as.factor(gmxPDO$Lag)
gmxPDO$Type[gmxPDO$Type=="Zlos"] <- "ZmLoss"
gmxPDO$Type[gmxPDO$Type=="Zoo"] <- "Zmeso"
gmxPDO$Type <- factor(gmxPDO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rgmxPDO <- gmxPDO[,c("Lag","Type","coef","p")]
rgmxPDO$p[rgmxPDO$p > 0.1] <- NA
rgmxPDO$sym <- rgmxPDO$p
rgmxPDO$sym[rgmxPDO$p <= 0.05] <- "*"
rgmxPDO$sym[rgmxPDO$p > 0.05] <- "."
rgmxPDO$sym[rgmxPDO$p > 0.1] <- NA

gmxNAO <- read.csv(paste0(datap,"GMX_regress_NAO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gmxNAO$Lag <- as.factor(gmxNAO$Lag)
gmxNAO$Type[gmxNAO$Type=="Zlos"] <- "ZmLoss"
gmxNAO$Type[gmxNAO$Type=="Zoo"] <- "Zmeso"
gmxNAO$Type <- factor(gmxNAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rgmxNAO <- gmxNAO[,c("Lag","Type","coef","p")]
rgmxNAO$p[rgmxNAO$p > 0.1] <- NA
rgmxNAO$sym <- rgmxNAO$p
rgmxNAO$sym[rgmxNAO$p <= 0.05] <- "*"
rgmxNAO$sym[rgmxNAO$p > 0.05] <- "."
rgmxNAO$sym[rgmxNAO$p > 0.1] <- NA

gmxAO <- read.csv(paste0(datap,"GMX_regress_AO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gmxAO$Lag <- as.factor(gmxAO$Lag)
gmxAO$Type[gmxAO$Type=="Zlos"] <- "ZmLoss"
gmxAO$Type[gmxAO$Type=="Zoo"] <- "Zmeso"
gmxAO$Type <- factor(gmxAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rgmxAO <- gmxAO[,c("Lag","Type","coef","p")]
rgmxAO$p[rgmxAO$p > 0.1] <- NA
rgmxAO$sym <- rgmxAO$p
rgmxAO$sym[rgmxAO$p <= 0.05] <- "*"
rgmxAO$sym[rgmxAO$p > 0.05] <- "."
rgmxAO$sym[rgmxAO$p > 0.1] <- NA


## NE
neNAO <- read.csv(paste0(datap,"NE_regress_NAO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
neNAO$Lag <- as.factor(neNAO$Lag)
neNAO$Type[neNAO$Type=="Zlos"] <- "ZmLoss"
neNAO$Type[neNAO$Type=="Zoo"] <- "Zmeso"
neNAO$Type <- factor(neNAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rneNAO <- neNAO[,c("Lag","Type","coef","p")]
rneNAO$p[rneNAO$p > 0.1] <- NA
rneNAO$sym <- rneNAO$p
rneNAO$sym[rneNAO$p <= 0.05] <- "*"
rneNAO$sym[rneNAO$p > 0.05] <- "."
rneNAO$sym[rneNAO$p > 0.1] <- NA

neAMO <- read.csv(paste0(datap,"NE_regress_AMO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
neAMO$Lag <- as.factor(neAMO$Lag)
neAMO$Type[neAMO$Type=="Zlos"] <- "ZmLoss"
neAMO$Type[neAMO$Type=="Zoo"] <- "Zmeso"
neAMO$Type <- factor(neAMO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rneAMO <- neAMO[,c("Lag","Type","coef","p")]
rneAMO$p[rneAMO$p > 0.1] <- NA
rneAMO$sym <- rneAMO$p
rneAMO$sym[rneAMO$p <= 0.05] <- "*"
rneAMO$sym[rneAMO$p > 0.05] <- "."
rneAMO$sym[rneAMO$p > 0.1] <- NA


##SE
seNAO <- read.csv(paste0(datap,"SE_regress_NAO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
seNAO$Lag <- as.factor(seNAO$Lag)
seNAO$Type[seNAO$Type=="Zlos"] <- "ZmLoss"
seNAO$Type[seNAO$Type=="Zoo"] <- "Zmeso"
seNAO$Type <- factor(seNAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rseNAO <- seNAO[,c("Lag","Type","coef","p")]
rseNAO$p[rseNAO$p > 0.1] <- NA
rseNAO$sym <- rseNAO$p
rseNAO$sym[rseNAO$p <= 0.05] <- "*"
rseNAO$sym[rseNAO$p > 0.05] <- "."
rseNAO$sym[rseNAO$p > 0.1] <- NA

seAMO <- read.csv(paste0(datap,"SE_regress_AMO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
seAMO$Lag <- as.factor(seAMO$Lag)
seAMO$Type[seAMO$Type=="Zlos"] <- "ZmLoss"
seAMO$Type[seAMO$Type=="Zoo"] <- "Zmeso"
seAMO$Type <- factor(seAMO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","ZmLoss","Zmeso","Tp"))
rseAMO <- seAMO[,c("Lag","Type","coef","p")]
rseAMO$p[rseAMO$p > 0.1] <- NA
rseAMO$sym <- rseAMO$p
rseAMO$sym[rseAMO$p <= 0.05] <- "*"
rseAMO$sym[rseAMO$p > 0.05] <- "."
rseAMO$sym[rseAMO$p > 0.1] <- NA



###====== Plots ===================================================================
## CCE
c5 <- ggplot(data = rccePDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="PDO\nCoef") +
  theme_minimal() + labs(x="")+ labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("CCE") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

c4 <- ggplot(data = rcceNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+ 
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3) 


## EBS
e5 <- ggplot(data = rebsPDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="PDO\nCoef") +
  theme_minimal() + labs(x="")+ labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("EBS") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

e4 <- ggplot(data = rebsNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+ 
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

e3 <- ggplot(data = rebsAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="AO\nCoef") +
  theme_minimal() + labs(x="")+ labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


## CHK
k5 <- ggplot(data = rchkPDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="PDO\nCoef") +
  theme_minimal() + labs(x="")+ 
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

k2 <- ggplot(data = rchkAMO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="AMO\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+  labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

k3 <- ggplot(data = rchkAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="AO\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("CHK") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

## NE
n1 <- ggplot(data = rneNAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="NAO\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

n2 <- ggplot(data = rneAMO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="AMO\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+ labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("NEUS") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


png(paste0(figp,'Heatmaps_CCE_EBS_GAK_NE_coef_select_climate_inputs_fish.png'), 
    width = 7.5*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c4,c5,'',
           e4,e5,e3,
           k5,k3,k2,
           n1,n2,
           nrow = 4, ncol = 3,labels = c("a","b","","c","d","e","f","g","h","i","j"), label_size = 11)
dev.off()


## AI
i5 <- ggplot(data = raiPDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="PDO\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("AI") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

i4 <- ggplot(data = raiNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

i3 <- ggplot(data = raiAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="AO\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


## GAK
g5 <- ggplot(data = rgakPDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="PDO\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("GAK") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

g4 <- ggplot(data = rgakNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

g3 <- ggplot(data = rgakAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="AO\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


## HI
h5 <- ggplot(data = rhiPDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="PDO\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("HI") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

h4 <- ggplot(data = rhiNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


## GMX
m5 <- ggplot(data = rgmxPDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="PDO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

m1 <- ggplot(data = rgmxNAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="NAO\nCoef") +
  theme_minimal() + ggtitle("") + labs(x="Lag (y)")+ labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

m3 <- ggplot(data = rgmxAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="AO\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("GMX") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


## SE
s1 <- ggplot(data = rseNAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="NAO\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+ 
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

s2 <- ggplot(data = rseAMO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.6,0.6),breaks=c(-0.5,-0.2,0,0.2,0.5),
                       labels=c(-0.5,-0.2,0,0.2,0.5), name="AMO\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+ labs(y="")+
  theme(legend.key.height=unit(0.5,'cm'),legend.key.width=unit(0.25, 'cm'))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("SEUS") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


## Put some all together for SUpp
png(paste0(figp,'Heatmaps_AI_GAK_HI_GMX_SE_coef_select_climate_inputs_fish.png'), 
    width = 7.5*300,        # 5 x 300 pixels
    height = 10*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( g4,g5,g3,
           i4,i5,i3,
           h4,h5,'',
           m5,m3,m1,
           s1,s2,'',
           nrow = 5, ncol = 3, labels = c("a","b","c",
                                          "d","e","f",
                                          "g","h","",
                                          "i","j","k",
                                          "l","m"), label_size = 11)
dev.off()









