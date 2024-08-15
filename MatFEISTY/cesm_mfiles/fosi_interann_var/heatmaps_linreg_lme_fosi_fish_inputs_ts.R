# Find patterns in coefelations between fish and inputs

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
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"

### ------------------------------------------------------------------------ ###
### ----------------------------- Heatmaps --------------------------------- ###
### ------------------------------------------------------------------------ ###

### Rearrange =====================================================
## CCE
cceTp <- read.csv(paste0(datap,"CCE_regress_Tp_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
cceTp$Lag <- as.factor(cceTp$Lag)
cceTp$Type <- factor(cceTp$Type, 
                      levels = c("B","D","A","P","F","L","M","S"))
rcceTp <- cceTp[,c("Lag","Type","coef","p")]
rcceTp$p[rcceTp$p > 0.1] <- NA
rcceTp$sym <- rcceTp$p
rcceTp$sym[rcceTp$p <= 0.05] <- "*"
rcceTp$sym[rcceTp$p > 0.05] <- "+"

cceTb <- read.csv(paste0(datap,"CCE_regress_Tb_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
cceTb$Lag <- as.factor(cceTb$Lag)
cceTb$Type <- factor(cceTb$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rcceTb <- cceTb[,c("Lag","Type","coef","p")]
rcceTb$p[rcceTb$p > 0.1] <- NA
rcceTb$sym <- rcceTb$p
rcceTb$sym[rcceTb$p <= 0.05] <- "*"
rcceTb$sym[rcceTb$p > 0.05] <- "+"

cceDet <- read.csv(paste0(datap,"CCE_regress_Det_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
cceDet$Lag <- as.factor(cceDet$Lag)
cceDet$Type <- factor(cceDet$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rcceDet <- cceDet[,c("Lag","Type","coef","p")]
rcceDet$p[rcceDet$p > 0.1] <- NA
rcceDet$sym <- rcceDet$p
rcceDet$sym[rcceDet$p <= 0.05] <- "*"
rcceDet$sym[rcceDet$p > 0.05] <- "+"

cceLZbiom <- read.csv(paste0(datap,"CCE_regress_LZbiom_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
cceLZbiom$Lag <- as.factor(cceLZbiom$Lag)
cceLZbiom$Type <- factor(cceLZbiom$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rcceLZbiom <- cceLZbiom[,c("Lag","Type","coef","p")]
rcceLZbiom$p[rcceLZbiom$p > 0.1] <- NA
rcceLZbiom$sym <- rcceLZbiom$p
rcceLZbiom$sym[rcceLZbiom$p <= 0.05] <- "*"
rcceLZbiom$sym[rcceLZbiom$p > 0.05] <- "+"

cceLZloss <- read.csv(paste0(datap,"CCE_regress_LZloss_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
cceLZloss$Lag <- as.factor(cceLZloss$Lag)
cceLZloss$Type <- factor(cceLZloss$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rcceLZloss <- cceLZloss[,c("Lag","Type","coef","p")]
rcceLZloss$p[rcceLZloss$p > 0.1] <- NA
rcceLZloss$sym <- rcceLZloss$p
rcceLZloss$sym[rcceLZloss$p <= 0.05] <- "*"
rcceLZloss$sym[rcceLZloss$p > 0.05] <- "+"


## EBS
ebsTp <- read.csv(paste0(datap,"EBS_regress_Tp_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
ebsTp$Lag <- as.factor(ebsTp$Lag)
ebsTp$Type <- factor(ebsTp$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rebsTp <- ebsTp[,c("Lag","Type","coef","p")]
rebsTp$p[rebsTp$p > 0.1] <- NA
rebsTp$sym <- rebsTp$p
rebsTp$sym[rebsTp$p <= 0.05] <- "*"
rebsTp$sym[rebsTp$p > 0.05] <- "+"

ebsTb <- read.csv(paste0(datap,"EBS_regress_Tb_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
ebsTb$Lag <- as.factor(ebsTb$Lag)
ebsTb$Type <- factor(ebsTb$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rebsTb <- ebsTb[,c("Lag","Type","coef","p")]
rebsTb$p[rebsTb$p > 0.1] <- NA
rebsTb$sym <- rebsTb$p
rebsTb$sym[rebsTb$p <= 0.05] <- "*"
rebsTb$sym[rebsTb$p > 0.05] <- "+"

ebsDet <- read.csv(paste0(datap,"EBS_regress_Det_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
ebsDet$Lag <- as.factor(ebsDet$Lag)
ebsDet$Type <- factor(ebsDet$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rebsDet <- ebsDet[,c("Lag","Type","coef","p")]
rebsDet$p[rebsDet$p > 0.1] <- NA
rebsDet$sym <- rebsDet$p
rebsDet$sym[rebsDet$p <= 0.05] <- "*"
rebsDet$sym[rebsDet$p > 0.05] <- "+"

ebsLZbiom <- read.csv(paste0(datap,"EBS_regress_LZbiom_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
ebsLZbiom$Lag <- as.factor(ebsLZbiom$Lag)
ebsLZbiom$Type <- factor(ebsLZbiom$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rebsLZbiom <- ebsLZbiom[,c("Lag","Type","coef","p")]
rebsLZbiom$p[rebsLZbiom$p > 0.1] <- NA
rebsLZbiom$sym <- rebsLZbiom$p
rebsLZbiom$sym[rebsLZbiom$p <= 0.05] <- "*"
rebsLZbiom$sym[rebsLZbiom$p > 0.05] <- "+"

ebsLZloss <- read.csv(paste0(datap,"EBS_regress_LZloss_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
ebsLZloss$Lag <- as.factor(ebsLZloss$Lag)
ebsLZloss$Type <- factor(ebsLZloss$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rebsLZloss <- ebsLZloss[,c("Lag","Type","coef","p")]
rebsLZloss$p[rebsLZloss$p > 0.1] <- NA
rebsLZloss$sym <- rebsLZloss$p
rebsLZloss$sym[rebsLZloss$p <= 0.05] <- "*"
rebsLZloss$sym[rebsLZloss$p > 0.05] <- "+"


## GAK
gakTp <- read.csv(paste0(datap,"GAK_regress_Tp_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gakTp$Lag <- as.factor(gakTp$Lag)
gakTp$Type <- factor(gakTp$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rgakTp <- gakTp[,c("Lag","Type","coef","p")]
rgakTp$p[rgakTp$p > 0.1] <- NA
rgakTp$sym <- rgakTp$p
rgakTp$sym[rgakTp$p <= 0.05] <- "*"
rgakTp$sym[rgakTp$p > 0.05] <- "+"

gakTb <- read.csv(paste0(datap,"GAK_regress_Tb_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gakTb$Lag <- as.factor(gakTb$Lag)
gakTb$Type <- factor(gakTb$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rgakTb <- gakTb[,c("Lag","Type","coef","p")]
rgakTb$p[rgakTb$p > 0.1] <- NA
rgakTb$sym <- rgakTb$p
rgakTb$sym[rgakTb$p <= 0.05] <- "*"
rgakTb$sym[rgakTb$p > 0.05] <- "+"

gakDet <- read.csv(paste0(datap,"GAK_regress_Det_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gakDet$Lag <- as.factor(gakDet$Lag)
gakDet$Type <- factor(gakDet$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rgakDet <- gakDet[,c("Lag","Type","coef","p")]
rgakDet$p[rgakDet$p > 0.1] <- NA
rgakDet$sym <- rgakDet$p
rgakDet$sym[rgakDet$p <= 0.05] <- "*"
rgakDet$sym[rgakDet$p > 0.05] <- "+"

gakLZbiom <- read.csv(paste0(datap,"GAK_regress_LZbiom_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gakLZbiom$Lag <- as.factor(gakLZbiom$Lag)
gakLZbiom$Type <- factor(gakLZbiom$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rgakLZbiom <- gakLZbiom[,c("Lag","Type","coef","p")]
rgakLZbiom$p[rgakLZbiom$p > 0.1] <- NA
rgakLZbiom$sym <- rgakLZbiom$p
rgakLZbiom$sym[rgakLZbiom$p <= 0.05] <- "*"
rgakLZbiom$sym[rgakLZbiom$p > 0.05] <- "+"

gakLZloss <- read.csv(paste0(datap,"GAK_regress_LZloss_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gakLZloss$Lag <- as.factor(gakLZloss$Lag)
gakLZloss$Type <- factor(gakLZloss$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rgakLZloss <- gakLZloss[,c("Lag","Type","coef","p")]
rgakLZloss$p[rgakLZloss$p > 0.1] <- NA
rgakLZloss$sym <- rgakLZloss$p
rgakLZloss$sym[rgakLZloss$p <= 0.05] <- "*"
rgakLZloss$sym[rgakLZloss$p > 0.05] <- "+"


## HI
hiTp <- read.csv(paste0(datap,"HI_regress_Tp_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
hiTp$Lag <- as.factor(hiTp$Lag)
hiTp$Type <- factor(hiTp$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rhiTp <- hiTp[,c("Lag","Type","coef","p")]
rhiTp$p[rhiTp$p > 0.1] <- NA
rhiTp$sym <- rhiTp$p
rhiTp$sym[rhiTp$p <= 0.05] <- "*"
rhiTp$sym[rhiTp$p > 0.05] <- "+"

hiTb <- read.csv(paste0(datap,"HI_regress_Tb_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
hiTb$Lag <- as.factor(hiTb$Lag)
hiTb$Type <- factor(hiTb$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rhiTb <- hiTb[,c("Lag","Type","coef","p")]
rhiTb$p[rhiTb$p > 0.1] <- NA
rhiTb$sym <- rhiTb$p
rhiTb$sym[rhiTb$p <= 0.05] <- "*"
rhiTb$sym[rhiTb$p > 0.05] <- "+"

hiDet <- read.csv(paste0(datap,"HI_regress_Det_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
hiDet$Lag <- as.factor(hiDet$Lag)
hiDet$Type <- factor(hiDet$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rhiDet <- hiDet[,c("Lag","Type","coef","p")]
rhiDet$p[rhiDet$p > 0.1] <- NA
rhiDet$sym <- rhiDet$p
rhiDet$sym[rhiDet$p <= 0.05] <- "*"
rhiDet$sym[rhiDet$p > 0.05] <- "+"

hiLZbiom <- read.csv(paste0(datap,"HI_regress_LZbiom_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
hiLZbiom$Lag <- as.factor(hiLZbiom$Lag)
hiLZbiom$Type <- factor(hiLZbiom$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rhiLZbiom <- hiLZbiom[,c("Lag","Type","coef","p")]
rhiLZbiom$p[rhiLZbiom$p > 0.1] <- NA
rhiLZbiom$sym <- rhiLZbiom$p
rhiLZbiom$sym[rhiLZbiom$p <= 0.05] <- "*"
rhiLZbiom$sym[rhiLZbiom$p > 0.05] <- "+"

hiLZloss <- read.csv(paste0(datap,"HI_regress_LZloss_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
hiLZloss$Lag <- as.factor(hiLZloss$Lag)
hiLZloss$Type <- factor(hiLZloss$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rhiLZloss <- hiLZloss[,c("Lag","Type","coef","p")]
rhiLZloss$p[rhiLZloss$p > 0.1] <- NA
rhiLZloss$sym <- rhiLZloss$p
rhiLZloss$sym[rhiLZloss$p <= 0.05] <- "*"
rhiLZloss$sym[rhiLZloss$p > 0.05] <- "+"


## CHK
chkTp <- read.csv(paste0(datap,"CHK_regress_Tp_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
chkTp$Lag <- as.factor(chkTp$Lag)
chkTp$Type <- factor(chkTp$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rchkTp <- chkTp[,c("Lag","Type","coef","p")]
rchkTp$p[rchkTp$p > 0.1] <- NA
rchkTp$sym <- rchkTp$p
rchkTp$sym[rchkTp$p <= 0.05] <- "*"
rchkTp$sym[rchkTp$p > 0.05] <- "+"

chkTb <- read.csv(paste0(datap,"CHK_regress_Tb_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
chkTb$Lag <- as.factor(chkTb$Lag)
chkTb$Type <- factor(chkTb$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rchkTb <- chkTb[,c("Lag","Type","coef","p")]
rchkTb$p[rchkTb$p > 0.1] <- NA
rchkTb$sym <- rchkTb$p
rchkTb$sym[rchkTb$p <= 0.05] <- "*"
rchkTb$sym[rchkTb$p > 0.05] <- "+"

chkDet <- read.csv(paste0(datap,"CHK_regress_Det_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
chkDet$Lag <- as.factor(chkDet$Lag)
chkDet$Type <- factor(chkDet$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rchkDet <- chkDet[,c("Lag","Type","coef","p")]
rchkDet$p[rchkDet$p > 0.1] <- NA
rchkDet$sym <- rchkDet$p
rchkDet$sym[rchkDet$p <= 0.05] <- "*"
rchkDet$sym[rchkDet$p > 0.05] <- "+"

chkLZbiom <- read.csv(paste0(datap,"CHK_regress_LZbiom_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
chkLZbiom$Lag <- as.factor(chkLZbiom$Lag)
chkLZbiom$Type <- factor(chkLZbiom$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rchkLZbiom <- chkLZbiom[,c("Lag","Type","coef","p")]
rchkLZbiom$p[rchkLZbiom$p > 0.1] <- NA
rchkLZbiom$sym <- rchkLZbiom$p
rchkLZbiom$sym[rchkLZbiom$p <= 0.05] <- "*"
rchkLZbiom$sym[rchkLZbiom$p > 0.05] <- "+"

chkLZloss <- read.csv(paste0(datap,"CHK_regress_LZloss_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
chkLZloss$Lag <- as.factor(chkLZloss$Lag)
chkLZloss$Type <- factor(chkLZloss$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rchkLZloss <- chkLZloss[,c("Lag","Type","coef","p")]
rchkLZloss$p[rchkLZloss$p > 0.1] <- NA
rchkLZloss$sym <- rchkLZloss$p
rchkLZloss$sym[rchkLZloss$p <= 0.05] <- "*"
rchkLZloss$sym[rchkLZloss$p > 0.05] <- "+"


##GMX
gmxTp <- read.csv(paste0(datap,"MX_regress_Tp_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gmxTp$Lag <- as.factor(gmxTp$Lag)
gmxTp$Type <- factor(gmxTp$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rgmxTp <- gmxTp[,c("Lag","Type","coef","p")]
rgmxTp$p[rgmxTp$p > 0.1] <- NA
rgmxTp$sym <- rgmxTp$p
rgmxTp$sym[rgmxTp$p <= 0.05] <- "*"
rgmxTp$sym[rgmxTp$p > 0.05] <- "+"

gmxTb <- read.csv(paste0(datap,"MX_regress_Tb_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gmxTb$Lag <- as.factor(gmxTb$Lag)
gmxTb$Type <- factor(gmxTb$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rgmxTb <- gmxTb[,c("Lag","Type","coef","p")]
rgmxTb$p[rgmxTb$p > 0.1] <- NA
rgmxTb$sym <- rgmxTb$p
rgmxTb$sym[rgmxTb$p <= 0.05] <- "*"
rgmxTb$sym[rgmxTb$p > 0.05] <- "+"

gmxDet <- read.csv(paste0(datap,"MX_regress_Det_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gmxDet$Lag <- as.factor(gmxDet$Lag)
gmxDet$Type <- factor(gmxDet$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rgmxDet <- gmxDet[,c("Lag","Type","coef","p")]
rgmxDet$p[rgmxDet$p > 0.1] <- NA
rgmxDet$sym <- rgmxDet$p
rgmxDet$sym[rgmxDet$p <= 0.05] <- "*"
rgmxDet$sym[rgmxDet$p > 0.05] <- "+"

gmxLZbiom <- read.csv(paste0(datap,"MX_regress_LZbiom_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gmxLZbiom$Lag <- as.factor(gmxLZbiom$Lag)
gmxLZbiom$Type <- factor(gmxLZbiom$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rgmxLZbiom <- gmxLZbiom[,c("Lag","Type","coef","p")]
rgmxLZbiom$p[rgmxLZbiom$p > 0.1] <- NA
rgmxLZbiom$sym <- rgmxLZbiom$p
rgmxLZbiom$sym[rgmxLZbiom$p <= 0.05] <- "*"
rgmxLZbiom$sym[rgmxLZbiom$p > 0.05] <- "+"

gmxLZloss <- read.csv(paste0(datap,"MX_regress_LZloss_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gmxLZloss$Lag <- as.factor(gmxLZloss$Lag)
gmxLZloss$Type <- factor(gmxLZloss$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rgmxLZloss <- gmxLZloss[,c("Lag","Type","coef","p")]
rgmxLZloss$p[rgmxLZloss$p > 0.1] <- NA
rgmxLZloss$sym <- rgmxLZloss$p
rgmxLZloss$sym[rgmxLZloss$p <= 0.05] <- "*"
rgmxLZloss$sym[rgmxLZloss$p > 0.05] <- "+"


## NE
neTp <- read.csv(paste0(datap,"NE_regress_Tp_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
neTp$Lag <- as.factor(neTp$Lag)
neTp$Type <- factor(neTp$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rneTp <- neTp[,c("Lag","Type","coef","p")]
rneTp$p[rneTp$p > 0.1] <- NA
rneTp$sym <- rneTp$p
rneTp$sym[rneTp$p <= 0.05] <- "*"
rneTp$sym[rneTp$p > 0.05] <- "+"

neTb <- read.csv(paste0(datap,"NE_regress_Tb_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
neTb$Lag <- as.factor(neTb$Lag)
neTb$Type <- factor(neTb$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rneTb <- neTb[,c("Lag","Type","coef","p")]
rneTb$p[rneTb$p > 0.1] <- NA
rneTb$sym <- rneTb$p
rneTb$sym[rneTb$p <= 0.05] <- "*"
rneTb$sym[rneTb$p > 0.05] <- "+"

neDet <- read.csv(paste0(datap,"NE_regress_Det_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
neDet$Lag <- as.factor(neDet$Lag)
neTb$Type <- factor(neTb$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rneDet <- neDet[,c("Lag","Type","coef","p")]
rneDet$p[rneDet$p > 0.1] <- NA
rneDet$sym <- rneDet$p
rneDet$sym[rneDet$p <= 0.05] <- "*"
rneDet$sym[rneDet$p > 0.05] <- "+"

neLZbiom <- read.csv(paste0(datap,"NE_regress_LZbiom_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
neLZbiom$Lag <- as.factor(neLZbiom$Lag)
neLZbiom$Type <- factor(neLZbiom$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rneLZbiom <- neLZbiom[,c("Lag","Type","coef","p")]
rneLZbiom$p[rneLZbiom$p > 0.1] <- NA
rneLZbiom$sym <- rneLZbiom$p
rneLZbiom$sym[rneLZbiom$p <= 0.05] <- "*"
rneLZbiom$sym[rneLZbiom$p > 0.05] <- "+"

neLZloss <- read.csv(paste0(datap,"NE_regress_LZloss_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
neLZloss$Lag <- as.factor(neLZloss$Lag)
neLZloss$Type <- factor(neLZloss$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rneLZloss <- neLZloss[,c("Lag","Type","coef","p")]
rneLZloss$p[rneLZloss$p > 0.1] <- NA
rneLZloss$sym <- rneLZloss$p
rneLZloss$sym[rneLZloss$p <= 0.05] <- "*"
rneLZloss$sym[rneLZloss$p > 0.05] <- "+"


##SE
seTp <- read.csv(paste0(datap,"SE_regress_Tp_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
seTp$Lag <- as.factor(seTp$Lag)
seTp$Type <- factor(seTp$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rseTp <- seTp[,c("Lag","Type","coef","p")]
rseTp$p[rseTp$p > 0.1] <- NA
rseTp$sym <- rseTp$p
rseTp$sym[rseTp$p <= 0.05] <- "*"
rseTp$sym[rseTp$p > 0.05] <- "+"

seTb <- read.csv(paste0(datap,"SE_regress_Tb_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
seTb$Lag <- as.factor(seTb$Lag)
seTb$Type <- factor(seTb$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rseTb <- seTb[,c("Lag","Type","coef","p")]
rseTb$p[rseTb$p > 0.1] <- NA
rseTb$sym <- rseTb$p
rseTb$sym[rseTb$p <= 0.05] <- "*"
rseTb$sym[rseTb$p > 0.05] <- "+"

seDet <- read.csv(paste0(datap,"SE_regress_Det_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
seDet$Lag <- as.factor(seDet$Lag)
seDet$Type <- factor(seDet$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rseDet <- seDet[,c("Lag","Type","coef","p")]
rseDet$p[rseDet$p > 0.1] <- NA
rseDet$sym <- rseDet$p
rseDet$sym[rseDet$p <= 0.05] <- "*"
rseDet$sym[rseDet$p > 0.05] <- "+"

seLZbiom <- read.csv(paste0(datap,"SE_regress_LZbiom_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
seLZbiom$Lag <- as.factor(seLZbiom$Lag)
seLZbiom$Type <- factor(seLZbiom$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rseLZbiom <- seLZbiom[,c("Lag","Type","coef","p")]
rseLZbiom$p[rseLZbiom$p > 0.1] <- NA
rseLZbiom$sym <- rseLZbiom$p
rseLZbiom$sym[rseLZbiom$p <= 0.05] <- "*"
rseLZbiom$sym[rseLZbiom$p > 0.05] <- "+"

seLZloss <- read.csv(paste0(datap,"SE_regress_LZloss_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
seLZloss$Lag <- as.factor(seLZloss$Lag)
seLZloss$Type <- factor(seLZloss$Type, 
                     levels = c("B","D","A","P","F","L","M","S"))
rseLZloss <- seLZloss[,c("Lag","Type","coef","p")]
rseLZloss$p[rseLZloss$p > 0.1] <- NA
rseLZloss$sym <- rseLZloss$p
rseLZloss$sym[rseLZloss$p <= 0.05] <- "*"
rseLZloss$sym[rseLZloss$p > 0.05] <- "+"


###====== Plots ===================================================================
## CCE
c5 <- ggplot(data = rcceTp, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.35,0.35),  
                       name="Tp\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3) 

c4 <- ggplot(data = rcceTb, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.35,0.35),  
                       name="Tb\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3)

c1 <- ggplot(data = rcceDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.35,0.35),  
                       name="Det\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3)

c2 <- ggplot(data = rcceLZbiom, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.35,0.35),  
                       name="LZbiom\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("CCE") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3)
  

c3 <- ggplot(data = rcceLZloss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.35,0.35),  
                       name="LZloss\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3)

png(paste0(figp,'Heatmaps_CCE_coef_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c1,c2,c3,c4,c5,
           nrow = 2, ncol = 3)
dev.off()


## EBS
e5 <- ggplot(data = rebsTp, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tp\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

e4 <- ggplot(data = rebsTb, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="Tb\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

e1 <- ggplot(data = rebsDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.955,0.955),  
                       name="Det\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

e2 <- ggplot(data = rebsLZbiom, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZbiom\nCoef") +
  theme_minimal() + ggtitle("EBS") + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

e3 <- ggplot(data = rebsLZloss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.95,0.95),  
                       name="LZloss\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_EBS_coef_inputs_fish.png'), 
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
g5 <- ggplot(data = rgakTp, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.2,0.2),  
                       name="Tp\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

g4 <- ggplot(data = rgakTb, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.2,0.2),  
                       name="Tb\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

g1 <- ggplot(data = rgakDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.2,0.2),  
                       name="Det\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

g2 <- ggplot(data = rgakLZbiom, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.2,0.2),  
                       name="LZbiom\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("GAK") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

g3 <- ggplot(data = rgakLZloss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.2,0.2),  
                       name="LZloss\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_GAK_coef_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( g1,g2,g3,g4,g5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## HI
h5 <- ggplot(data = rhiTp, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.20,0.20),  
                       name="Tp\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

h4 <- ggplot(data = rhiTb, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.20,0.20),  
                       name="Tb\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

h1 <- ggplot(data = rhiDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.20,0.20),  
                       name="Det\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

h2 <- ggplot(data = rhiLZbiom, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.20,0.20),  
                       name="LZbiom\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("HI") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

h3 <- ggplot(data = rhiLZloss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.20,0.20),  
                       name="LZloss\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_HI_coef_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( h1,h2,h3,h4,h5,
           nrow = 2, ncol = 3)#,
           #align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## CHK
k5 <- ggplot(data = rchkTp, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.90,0.90),  
                       name="Tp\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

k4 <- ggplot(data = rchkTb, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.90,0.90),  
                       name="Tb\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

k1 <- ggplot(data = rchkDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.91,0.91),  
                       name="Det\nCoef") +
  theme_minimal() + ggtitle("") + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

k2 <- ggplot(data = rchkLZbiom, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.90,0.90),  
                       name="LZbiom\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("CHK") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

k3 <- ggplot(data = rchkLZloss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.90,0.90),  
                       name="LZloss\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_CHK_coef_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( k1,k2,k3,k4,k5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## GMX
m5 <- ggplot(data = rgmxTp, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.20,0.20),  
                       name="Tp\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

m4 <- ggplot(data = rgmxTb, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.20,0.20),  
                       name="Tb\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

m1 <- ggplot(data = rgmxDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.20,0.20),  
                       name="Det\nCoef") +
  theme_minimal() + ggtitle("") + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

m2 <- ggplot(data = rgmxLZbiom, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.99,0.99),  
                       name="LZbiom\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("GMX") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

m3 <- ggplot(data = rgmxLZloss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.20,0.20),  
                       name="LZloss\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_GMX_coef_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( m1,m2,m3,m4,m5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## SE
s5 <- ggplot(data = rseTp, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.35,0.35),  
                       name="Tp\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

s4 <- ggplot(data = rseTb, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.35,0.35),  
                       name="Tb\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

s1 <- ggplot(data = rseDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.35,0.35),  
                       name="Det\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

s2 <- ggplot(data = rseLZbiom, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.35,0.35),  
                       name="LZbiom\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("SEUS") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

s3 <- ggplot(data = rseLZloss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.35,0.35),  
                       name="LZloss\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_SE_coef_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( s1,s2,s3,s4,s5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## NE
n5 <- ggplot(data = rneTp, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.60,0.60),  
                       name="Tp\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

n4 <- ggplot(data = rneTb, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.60,0.60),  
                       name="Tb\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

n1 <- ggplot(data = rneDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.60,0.60),  
                       name="Det\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

n2 <- ggplot(data = rneLZbiom, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.60,0.60),  
                       name="LZbiom\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("NEUS") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

n3 <- ggplot(data = rneLZloss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.60,0.60),  
                       name="LZloss\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_NE_coef_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( n1,n2,n3,n4,n5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()






