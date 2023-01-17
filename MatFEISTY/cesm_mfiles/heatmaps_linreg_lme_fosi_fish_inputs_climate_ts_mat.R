# Find patterns in correlations between climate with fish and inputs

rm(list=ls())

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/corrs/"

library(Hmisc)
library(ggplot2)
library(reshape2) #melt
library(cowplot) #plot_grid
library(gridExtra)
library(coefgram)
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
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/FOSI/"

### ------------------------------------------------------------------------ ###
### ----------------------------- Heatmaps --------------------------------- ###
### ------------------------------------------------------------------------ ###

### Rearrange =====================================================
## CCE
ccePDO <- read.csv(paste0(datap,"CCE_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
ccePDO$Lag <- as.factor(ccePDO$Lag)
ccePDO$Type <- factor(ccePDO$Type, 
                         levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rccePDO <- ccePDO[,c("Lag","Type","coef","p")]
rccePDO$p[rccePDO$p > 0.1] <- NA
rccePDO$sym <- rccePDO$p
rccePDO$sym[rccePDO$p <= 0.05] <- "*"
rccePDO$sym[rccePDO$p > 0.05] <- "+"
rccePDO$sym[rccePDO$p > 0.1] <- NA

cceNino34 <- read.csv(paste0(datap,"CCE_regress_Nino34_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
cceNino34$Lag <- as.factor(cceNino34$Lag)
cceNino34$Type <- factor(cceNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rcceNino34 <- cceNino34[,c("Lag","Type","coef","p")]
rcceNino34$p[rcceNino34$p > 0.1] <- NA
rcceNino34$sym <- rcceNino34$p
rcceNino34$sym[rcceNino34$p <= 0.05] <- "*"
rcceNino34$sym[rcceNino34$p > 0.05] <- "+"
rcceNino34$sym[rcceNino34$p > 0.1] <- NA

cceNAO <- read.csv(paste0(datap,"CCE_regress_NAO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
cceNAO$Lag <- as.factor(cceNAO$Lag)
cceNAO$Type <- factor(cceNAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rcceNAO <- cceNAO[,c("Lag","Type","coef","p")]
rcceNAO$p[rcceNAO$p > 0.1] <- NA
rcceNAO$sym <- rcceNAO$p
rcceNAO$sym[rcceNAO$p <= 0.05] <- "*"
rcceNAO$sym[rcceNAO$p > 0.05] <- "+"
rcceNAO$sym[rcceNAO$p > 0.1] <- NA

cceAMO <- read.csv(paste0(datap,"CCE_regress_AMO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
cceAMO$Lag <- as.factor(cceAMO$Lag)
cceAMO$Type <- factor(cceAMO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rcceAMO <- cceAMO[,c("Lag","Type","coef","p")]
rcceAMO$p[rcceAMO$p > 0.1] <- NA
rcceAMO$sym <- rcceAMO$p
rcceAMO$sym[rcceAMO$p <= 0.05] <- "*"
rcceAMO$sym[rcceAMO$p > 0.05] <- "+"
rcceAMO$sym[rcceAMO$p > 0.1] <- NA

cceAO <- read.csv(paste0(datap,"CCE_regress_AO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
cceAO$Lag <- as.factor(cceAO$Lag)
cceAO$Type <- factor(cceAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rcceAO <- cceAO[,c("Lag","Type","coef","p")]
rcceAO$p[rcceAO$p > 0.1] <- NA
rcceAO$sym <- rcceAO$p
rcceAO$sym[rcceAO$p <= 0.05] <- "*"
rcceAO$sym[rcceAO$p > 0.05] <- "+"
rcceAO$sym[rcceAO$p > 0.1] <- NA


## EBS
ebsPDO <- read.csv(paste0(datap,"EBS_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
ebsPDO$Lag <- as.factor(ebsPDO$Lag)
ebsPDO$Type <- factor(ebsPDO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rebsPDO <- ebsPDO[,c("Lag","Type","coef","p")]
rebsPDO$p[rebsPDO$p > 0.1] <- NA
rebsPDO$sym <- rebsPDO$p
rebsPDO$sym[rebsPDO$p <= 0.05] <- "*"
rebsPDO$sym[rebsPDO$p > 0.05] <- "+"
rebsPDO$sym[rebsPDO$p > 0.1] <- NA

ebsNino34 <- read.csv(paste0(datap,"EBS_regress_Nino34_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
ebsNino34$Lag <- as.factor(ebsNino34$Lag)
ebsNino34$Type <- factor(ebsNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rebsNino34 <- ebsNino34[,c("Lag","Type","coef","p")]
rebsNino34$p[rebsNino34$p > 0.1] <- NA
rebsNino34$sym <- rebsNino34$p
rebsNino34$sym[rebsNino34$p <= 0.05] <- "*"
rebsNino34$sym[rebsNino34$p > 0.05] <- "+"
rebsNino34$sym[rebsNino34$p > 0.1] <- NA

ebsNAO <- read.csv(paste0(datap,"EBS_regress_NAO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
ebsNAO$Lag <- as.factor(ebsNAO$Lag)
ebsNAO$Type <- factor(ebsNAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rebsNAO <- ebsNAO[,c("Lag","Type","coef","p")]
rebsNAO$p[rebsNAO$p > 0.1] <- NA
rebsNAO$sym <- rebsNAO$p
rebsNAO$sym[rebsNAO$p <= 0.05] <- "*"
rebsNAO$sym[rebsNAO$p > 0.05] <- "+"
rebsNAO$sym[rebsNAO$p > 0.1] <- NA

ebsAMO <- read.csv(paste0(datap,"EBS_regress_AMO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
ebsAMO$Lag <- as.factor(ebsAMO$Lag)
ebsAMO$Type <- factor(ebsAMO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rebsAMO <- ebsAMO[,c("Lag","Type","coef","p")]
rebsAMO$p[rebsAMO$p > 0.1] <- NA
rebsAMO$sym <- rebsAMO$p
rebsAMO$sym[rebsAMO$p <= 0.05] <- "*"
rebsAMO$sym[rebsAMO$p > 0.05] <- "+"
rebsAMO$sym[rebsAMO$p > 0.1] <- NA

ebsAO <- read.csv(paste0(datap,"EBS_regress_AO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
ebsAO$Lag <- as.factor(ebsAO$Lag)
ebsAO$Type <- factor(ebsAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rebsAO <- ebsAO[,c("Lag","Type","coef","p")]
rebsAO$p[rebsAO$p > 0.1] <- NA
rebsAO$sym <- rebsAO$p
rebsAO$sym[rebsAO$p <= 0.05] <- "*"
rebsAO$sym[rebsAO$p > 0.05] <- "+"
rebsAO$sym[rebsAO$p > 0.1] <- NA


## GAK
gakPDO <- read.csv(paste0(datap,"GAK_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gakPDO$Lag <- as.factor(gakPDO$Lag)
gakPDO$Type <- factor(gakPDO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rgakPDO <- gakPDO[,c("Lag","Type","coef","p")]
rgakPDO$p[rgakPDO$p > 0.1] <- NA
rgakPDO$sym <- rgakPDO$p
rgakPDO$sym[rgakPDO$p <= 0.05] <- "*"
rgakPDO$sym[rgakPDO$p > 0.05] <- "+"
rgakPDO$sym[rgakPDO$p > 0.1] <- NA

gakNino34 <- read.csv(paste0(datap,"GAK_regress_Nino34_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gakNino34$Lag <- as.factor(gakNino34$Lag)
gakNino34$Type <- factor(gakNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rgakNino34 <- gakNino34[,c("Lag","Type","coef","p")]
rgakNino34$p[rgakNino34$p > 0.1] <- NA
rgakNino34$sym <- rgakNino34$p
rgakNino34$sym[rgakNino34$p <= 0.05] <- "*"
rgakNino34$sym[rgakNino34$p > 0.05] <- "+"
rgakNino34$sym[rgakNino34$p > 0.1] <- NA

gakNAO <- read.csv(paste0(datap,"GAK_regress_NAO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gakNAO$Lag <- as.factor(gakNAO$Lag)
gakNAO$Type <- factor(gakNAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rgakNAO <- gakNAO[,c("Lag","Type","coef","p")]
rgakNAO$p[rgakNAO$p > 0.1] <- NA
rgakNAO$sym <- rgakNAO$p
rgakNAO$sym[rgakNAO$p <= 0.05] <- "*"
rgakNAO$sym[rgakNAO$p > 0.05] <- "+"
rgakNAO$sym[rgakNAO$p > 0.1] <- NA

gakAMO <- read.csv(paste0(datap,"GAK_regress_AMO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gakAMO$Lag <- as.factor(gakAMO$Lag)
gakAMO$Type <- factor(gakAMO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rgakAMO <- gakAMO[,c("Lag","Type","coef","p")]
rgakAMO$p[rgakAMO$p > 0.1] <- NA
rgakAMO$sym <- rgakAMO$p
rgakAMO$sym[rgakAMO$p <= 0.05] <- "*"
rgakAMO$sym[rgakAMO$p > 0.05] <- "+"
rgakAMO$sym[rgakAMO$p > 0.1] <- NA

gakAO <- read.csv(paste0(datap,"GAK_regress_AO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gakAO$Lag <- as.factor(gakAO$Lag)
gakAO$Type <- factor(gakAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rgakAO <- gakAO[,c("Lag","Type","coef","p")]
rgakAO$p[rgakAO$p > 0.1] <- NA
rgakAO$sym <- rgakAO$p
rgakAO$sym[rgakAO$p <= 0.05] <- "*"
rgakAO$sym[rgakAO$p > 0.05] <- "+"
rgakAO$sym[rgakAO$p > 0.1] <- NA


## HI
hiPDO <- read.csv(paste0(datap,"HI_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
hiPDO$Lag <- as.factor(hiPDO$Lag)
hiPDO$Type <- factor(hiPDO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rhiPDO <- hiPDO[,c("Lag","Type","coef","p")]
rhiPDO$p[rhiPDO$p > 0.1] <- NA
rhiPDO$sym <- rhiPDO$p
rhiPDO$sym[rhiPDO$p <= 0.05] <- "*"
rhiPDO$sym[rhiPDO$p > 0.05] <- "+"
rhiPDO$sym[rhiPDO$p > 0.1] <- NA

hiNino34 <- read.csv(paste0(datap,"HI_regress_Nino34_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
hiNino34$Lag <- as.factor(hiNino34$Lag)
hiNino34$Type <- factor(hiNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rhiNino34 <- hiNino34[,c("Lag","Type","coef","p")]
rhiNino34$p[rhiNino34$p > 0.1] <- NA
rhiNino34$sym <- rhiNino34$p
rhiNino34$sym[rhiNino34$p <= 0.05] <- "*"
rhiNino34$sym[rhiNino34$p > 0.05] <- "+"
rhiNino34$sym[rhiNino34$p > 0.1] <- NA

hiNAO <- read.csv(paste0(datap,"HI_regress_NAO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
hiNAO$Lag <- as.factor(hiNAO$Lag)
hiNAO$Type <- factor(hiNAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rhiNAO <- hiNAO[,c("Lag","Type","coef","p")]
rhiNAO$p[rhiNAO$p > 0.1] <- NA
rhiNAO$sym <- rhiNAO$p
rhiNAO$sym[rhiNAO$p <= 0.05] <- "*"
rhiNAO$sym[rhiNAO$p > 0.05] <- "+"
rhiNAO$sym[rhiNAO$p > 0.1] <- NA

hiAMO <- read.csv(paste0(datap,"HI_regress_AMO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
hiAMO$Lag <- as.factor(hiAMO$Lag)
hiAMO$Type <- factor(hiAMO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rhiAMO <- hiAMO[,c("Lag","Type","coef","p")]
rhiAMO$p[rhiAMO$p > 0.1] <- NA
rhiAMO$sym <- rhiAMO$p
rhiAMO$sym[rhiAMO$p <= 0.05] <- "*"
rhiAMO$sym[rhiAMO$p > 0.05] <- "+"
rhiAMO$sym[rhiAMO$p > 0.1] <- NA

hiAO <- read.csv(paste0(datap,"HI_regress_AO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
hiAO$Lag <- as.factor(hiAO$Lag)
hiAO$Type <- factor(hiAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rhiAO <- hiAO[,c("Lag","Type","coef","p")]
rhiAO$p[rhiAO$p > 0.1] <- NA
rhiAO$sym <- rhiAO$p
rhiAO$sym[rhiAO$p <= 0.05] <- "*"
rhiAO$sym[rhiAO$p > 0.05] <- "+"
rhiAO$sym[rhiAO$p > 0.1] <- NA


## CHK
chkPDO <- read.csv(paste0(datap,"CHK_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
chkPDO$Lag <- as.factor(chkPDO$Lag)
chkPDO$Type <- factor(chkPDO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rchkPDO <- chkPDO[,c("Lag","Type","coef","p")]
rchkPDO$p[rchkPDO$p > 0.1] <- NA
rchkPDO$sym <- rchkPDO$p
rchkPDO$sym[rchkPDO$p <= 0.05] <- "*"
rchkPDO$sym[rchkPDO$p > 0.05] <- "+"
rchkPDO$sym[rchkPDO$p > 0.1] <- NA

chkNino34 <- read.csv(paste0(datap,"CHK_regress_Nino34_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
chkNino34$Lag <- as.factor(chkNino34$Lag)
chkNino34$Type <- factor(chkNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rchkNino34 <- chkNino34[,c("Lag","Type","coef","p")]
rchkNino34$p[rchkNino34$p > 0.1] <- NA
rchkNino34$sym <- rchkNino34$p
rchkNino34$sym[rchkNino34$p <= 0.05] <- "*"
rchkNino34$sym[rchkNino34$p > 0.05] <- "+"
rchkNino34$sym[rchkNino34$p > 0.1] <- NA

chkNAO <- read.csv(paste0(datap,"CHK_regress_NAO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
chkNAO$Lag <- as.factor(chkNAO$Lag)
chkNAO$Type <- factor(chkNAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rchkNAO <- chkNAO[,c("Lag","Type","coef","p")]
rchkNAO$p[rchkNAO$p > 0.1] <- NA
rchkNAO$sym <- rchkNAO$p
rchkNAO$sym[rchkNAO$p <= 0.05] <- "*"
rchkNAO$sym[rchkNAO$p > 0.05] <- "+"
rchkNAO$sym[rchkNAO$p > 0.1] <- NA

chkAMO <- read.csv(paste0(datap,"CHK_regress_AMO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
chkAMO$Lag <- as.factor(chkAMO$Lag)
chkAMO$Type <- factor(chkAMO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rchkAMO <- chkAMO[,c("Lag","Type","coef","p")]
rchkAMO$p[rchkAMO$p > 0.1] <- NA
rchkAMO$sym <- rchkAMO$p
rchkAMO$sym[rchkAMO$p <= 0.05] <- "*"
rchkAMO$sym[rchkAMO$p > 0.05] <- "+"
rchkAMO$sym[rchkAMO$p > 0.1] <- NA

chkAO <- read.csv(paste0(datap,"CHK_regress_AO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
chkAO$Lag <- as.factor(chkAO$Lag)
chkAO$Type <- factor(chkAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rchkAO <- chkAO[,c("Lag","Type","coef","p")]
rchkAO$p[rchkAO$p > 0.1] <- NA
rchkAO$sym <- rchkAO$p
rchkAO$sym[rchkAO$p <= 0.05] <- "*"
rchkAO$sym[rchkAO$p > 0.05] <- "+"
rchkAO$sym[rchkAO$p > 0.1] <- NA


##GMX
gmxPDO <- read.csv(paste0(datap,"GMX_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gmxPDO$Lag <- as.factor(gmxPDO$Lag)
gmxPDO$Type <- factor(gmxPDO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rgmxPDO <- gmxPDO[,c("Lag","Type","coef","p")]
rgmxPDO$p[rgmxPDO$p > 0.1] <- NA
rgmxPDO$sym <- rgmxPDO$p
rgmxPDO$sym[rgmxPDO$p <= 0.05] <- "*"
rgmxPDO$sym[rgmxPDO$p > 0.05] <- "+"
rgmxPDO$sym[rgmxPDO$p > 0.1] <- NA

gmxNino34 <- read.csv(paste0(datap,"GMX_regress_Nino34_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gmxNino34$Lag <- as.factor(gmxNino34$Lag)
gmxNino34$Type <- factor(gmxNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rgmxNino34 <- gmxNino34[,c("Lag","Type","coef","p")]
rgmxNino34$p[rgmxNino34$p > 0.1] <- NA
rgmxNino34$sym <- rgmxNino34$p
rgmxNino34$sym[rgmxNino34$p <= 0.05] <- "*"
rgmxNino34$sym[rgmxNino34$p > 0.05] <- "+"
rgmxNino34$sym[rgmxNino34$p > 0.1] <- NA

gmxNAO <- read.csv(paste0(datap,"GMX_regress_NAO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gmxNAO$Lag <- as.factor(gmxNAO$Lag)
gmxNAO$Type <- factor(gmxNAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rgmxNAO <- gmxNAO[,c("Lag","Type","coef","p")]
rgmxNAO$p[rgmxNAO$p > 0.1] <- NA
rgmxNAO$sym <- rgmxNAO$p
rgmxNAO$sym[rgmxNAO$p <= 0.05] <- "*"
rgmxNAO$sym[rgmxNAO$p > 0.05] <- "+"
rgmxNAO$sym[rgmxNAO$p > 0.1] <- NA

gmxAMO <- read.csv(paste0(datap,"GMX_regress_AMO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gmxAMO$Lag <- as.factor(gmxAMO$Lag)
gmxAMO$Type <- factor(gmxAMO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rgmxAMO <- gmxAMO[,c("Lag","Type","coef","p")]
rgmxAMO$p[rgmxAMO$p > 0.1] <- NA
rgmxAMO$sym <- rgmxAMO$p
rgmxAMO$sym[rgmxAMO$p <= 0.05] <- "*"
rgmxAMO$sym[rgmxAMO$p > 0.05] <- "+"
rgmxAMO$sym[rgmxAMO$p > 0.1] <- NA

gmxAO <- read.csv(paste0(datap,"GMX_regress_AO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
gmxAO$Lag <- as.factor(gmxAO$Lag)
gmxAO$Type <- factor(gmxAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rgmxAO <- gmxAO[,c("Lag","Type","coef","p")]
rgmxAO$p[rgmxAO$p > 0.1] <- NA
rgmxAO$sym <- rgmxAO$p
rgmxAO$sym[rgmxAO$p <= 0.05] <- "*"
rgmxAO$sym[rgmxAO$p > 0.05] <- "+"
rgmxAO$sym[rgmxAO$p > 0.1] <- NA


## NE
nePDO <- read.csv(paste0(datap,"NE_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
nePDO$Lag <- as.factor(nePDO$Lag)
nePDO$Type <- factor(nePDO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rnePDO <- nePDO[,c("Lag","Type","coef","p")]
rnePDO$p[rnePDO$p > 0.1] <- NA
rnePDO$sym <- rnePDO$p
rnePDO$sym[rnePDO$p <= 0.05] <- "*"
rnePDO$sym[rnePDO$p > 0.05] <- "+"
rnePDO$sym[rnePDO$p > 0.1] <- NA

neNino34 <- read.csv(paste0(datap,"NE_regress_Nino34_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
neNino34$Lag <- as.factor(neNino34$Lag)
neNino34$Type <- factor(neNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rneNino34 <- neNino34[,c("Lag","Type","coef","p")]
rneNino34$p[rneNino34$p > 0.1] <- NA
rneNino34$sym <- rneNino34$p
rneNino34$sym[rneNino34$p <= 0.05] <- "*"
rneNino34$sym[rneNino34$p > 0.05] <- "+"
rneNino34$sym[rneNino34$p > 0.1] <- NA

neNAO <- read.csv(paste0(datap,"NE_regress_NAO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
neNAO$Lag <- as.factor(neNAO$Lag)
neNAO$Type <- factor(neNAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rneNAO <- neNAO[,c("Lag","Type","coef","p")]
rneNAO$p[rneNAO$p > 0.1] <- NA
rneNAO$sym <- rneNAO$p
rneNAO$sym[rneNAO$p <= 0.05] <- "*"
rneNAO$sym[rneNAO$p > 0.05] <- "+"
rneNAO$sym[rneNAO$p > 0.1] <- NA

neAMO <- read.csv(paste0(datap,"NE_regress_AMO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
neAMO$Lag <- as.factor(neAMO$Lag)
neAMO$Type <- factor(neAMO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rneAMO <- neAMO[,c("Lag","Type","coef","p")]
rneAMO$p[rneAMO$p > 0.1] <- NA
rneAMO$sym <- rneAMO$p
rneAMO$sym[rneAMO$p <= 0.05] <- "*"
rneAMO$sym[rneAMO$p > 0.05] <- "+"
rneAMO$sym[rneAMO$p > 0.1] <- NA

neAO <- read.csv(paste0(datap,"NE_regress_AO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
neAO$Lag <- as.factor(neAO$Lag)
neAO$Type <- factor(neAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rneAO <- neAO[,c("Lag","Type","coef","p")]
rneAO$p[rneAO$p > 0.1] <- NA
rneAO$sym <- rneAO$p
rneAO$sym[rneAO$p <= 0.05] <- "*"
rneAO$sym[rneAO$p > 0.05] <- "+"
rneAO$sym[rneAO$p > 0.1] <- NA


##SE
sePDO <- read.csv(paste0(datap,"SE_regress_PDO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
sePDO$Lag <- as.factor(sePDO$Lag)
sePDO$Type <- factor(sePDO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rsePDO <- sePDO[,c("Lag","Type","coef","p")]
rsePDO$p[rsePDO$p > 0.1] <- NA
rsePDO$sym <- rsePDO$p
rsePDO$sym[rsePDO$p <= 0.05] <- "*"
rsePDO$sym[rsePDO$p > 0.05] <- "+"
rsePDO$sym[rsePDO$p > 0.1] <- NA

seNino34 <- read.csv(paste0(datap,"SE_regress_Nino34_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
seNino34$Lag <- as.factor(seNino34$Lag)
seNino34$Type <- factor(seNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rseNino34 <- seNino34[,c("Lag","Type","coef","p")]
rseNino34$p[rseNino34$p > 0.1] <- NA
rseNino34$sym <- rseNino34$p
rseNino34$sym[rseNino34$p <= 0.05] <- "*"
rseNino34$sym[rseNino34$p > 0.05] <- "+"
rseNino34$sym[rseNino34$p > 0.1] <- NA

seNAO <- read.csv(paste0(datap,"SE_regress_NAO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
seNAO$Lag <- as.factor(seNAO$Lag)
seNAO$Type <- factor(seNAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rseNAO <- seNAO[,c("Lag","Type","coef","p")]
rseNAO$p[rseNAO$p > 0.1] <- NA
rseNAO$sym <- rseNAO$p
rseNAO$sym[rseNAO$p <= 0.05] <- "*"
rseNAO$sym[rseNAO$p > 0.05] <- "+"
rseNAO$sym[rseNAO$p > 0.1] <- NA

seAMO <- read.csv(paste0(datap,"SE_regress_AMO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
seAMO$Lag <- as.factor(seAMO$Lag)
seAMO$Type <- factor(seAMO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rseAMO <- seAMO[,c("Lag","Type","coef","p")]
rseAMO$p[rseAMO$p > 0.1] <- NA
rseAMO$sym <- rseAMO$p
rseAMO$sym[rseAMO$p <= 0.05] <- "*"
rseAMO$sym[rseAMO$p > 0.05] <- "+"
rseAMO$sym[rseAMO$p > 0.1] <- NA

seAO <- read.csv(paste0(datap,"SE_regress_AO_div2SD_melt_mat.csv"),sep=",",header = T,stringsAsFactors = F)
seAO$Lag <- as.factor(seAO$Lag)
seAO$Type <- factor(seAO$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","Zlos","Zoo","Tp"))
rseAO <- seAO[,c("Lag","Type","coef","p")]
rseAO$p[rseAO$p > 0.1] <- NA
rseAO$sym <- rseAO$p
rseAO$sym[rseAO$p <= 0.05] <- "*"
rseAO$sym[rseAO$p > 0.05] <- "+"
rseAO$sym[rseAO$p > 0.1] <- NA


###====== Plots ===================================================================
## CCE
c5 <- ggplot(data = rccePDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.80,0.80),  
                       name="PDO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3) 

c4 <- ggplot(data = rcceNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.80,0.80),  
                       name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3) 

c1 <- ggplot(data = rcceNAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.80,0.80),  
                       name="NAO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)

c2 <- ggplot(data = rcceAMO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.80,0.80),  
                       name="AMO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("CCE") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3) 
  

c3 <- ggplot(data = rcceAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.80,0.80),  
                       name="AO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_CCE_coef_climate_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c1,c2,c3,c4,c5,
           nrow = 2, ncol = 3)
dev.off()


## EBS
e5 <- ggplot(data = rebsPDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="PDO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

e4 <- ggplot(data = rebsNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

e1 <- ggplot(data = rebsNAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.955,0.955),  
                       name="NAO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

e2 <- ggplot(data = rebsAMO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="AMO\nCoef") +
  theme_minimal() + ggtitle("EBS") + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

e3 <- ggplot(data = rebsAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="AO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_EBS_coef_climate_inputs_fish.png'), 
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
g5 <- ggplot(data = rgakPDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

g4 <- ggplot(data = rgakNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

g1 <- ggplot(data = rgakNAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NAO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

g2 <- ggplot(data = rgakAMO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AMO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("GAK") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

g3 <- ggplot(data = rgakAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_GAK_coef_climate_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( g1,g2,g3,g4,g5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## HI
h5 <- ggplot(data = rhiPDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.50,0.50),  
                       name="PDO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

h4 <- ggplot(data = rhiNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.50,0.50),  
                       name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

h1 <- ggplot(data = rhiNAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.50,0.50),  
                       name="NAO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

h2 <- ggplot(data = rhiAMO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.50,0.50),  
                       name="AMO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("HI") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

h3 <- ggplot(data = rhiAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.50,0.50),  
                       name="AO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_HI_coef_climate_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( h1,h2,h3,h4,h5,
           nrow = 2, ncol = 3)#,
           #align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## CHK
k5 <- ggplot(data = rchkPDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.90,0.90),  
                       name="PDO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

k4 <- ggplot(data = rchkNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.90,0.90),  
                       name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

k1 <- ggplot(data = rchkNAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.91,0.91),  
                       name="NAO\nCoef") +
  theme_minimal() + ggtitle("") + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

k2 <- ggplot(data = rchkAMO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.90,0.90),  
                       name="AMO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("CHK") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

k3 <- ggplot(data = rchkAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.90,0.90),  
                       name="AO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_CHK_coef_climate_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( k1,k2,k3,k4,k5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## GMX
m5 <- ggplot(data = rgmxPDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.45,0.45),  
                       name="PDO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

m4 <- ggplot(data = rgmxNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.45,0.45),  
                       name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

m1 <- ggplot(data = rgmxNAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.45,0.45),  
                       name="NAO\nCoef") +
  theme_minimal() + ggtitle("") + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

m2 <- ggplot(data = rgmxAMO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.99,0.99),  
                       name="AMO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("GMX") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

m3 <- ggplot(data = rgmxAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.45,0.45),  
                       name="AO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_GMX_coef_climate_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( m1,m2,m3,m4,m5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## SE
s5 <- ggplot(data = rsePDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.55,0.55),  
                       name="PDO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

s4 <- ggplot(data = rseNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.55,0.55),  
                       name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

s1 <- ggplot(data = rseNAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.55,0.55),  
                       name="NAO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

s2 <- ggplot(data = rseAMO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.55,0.55),  
                       name="AMO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("SEUS") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

s3 <- ggplot(data = rseAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.55,0.55),  
                       name="AO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_SE_coef_climate_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( s1,s2,s3,s4,s5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## NE
n5 <- ggplot(data = rnePDO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.25,1.25),  
                       name="PDO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

n4 <- ggplot(data = rneNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.25,1.25),  
                       name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

n1 <- ggplot(data = rneNAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.25,1.25),  
                       name="NAO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

n2 <- ggplot(data = rneAMO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.25,1.25),  
                       name="AMO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("NEUS") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

n3 <- ggplot(data = rneAO, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.25,1.25),  
                       name="AO\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_NE_coef_climate_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( n1,n2,n3,n4,n5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()






