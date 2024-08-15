# Find patterns in correlations between fish and inputs
# Nu instead of biomass
# ZmLoss & Det units annual nu

rm(list=ls())

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/corrs/"

library(ggplot2)
library(reshape2) #melt
library(cowplot) #plot_grid
library(gridExtra)
library(ggpubr) #cowplot and gridExtra
library(gplots)
#library("coefplot")
#library(scico)

### --------------------------------------------------------------
# load data
#FEISTY
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"

### ------------------------------------------------------------------------ ###
### ----------------------------- Heatmaps --------------------------------- ###
### ------------------------------------------------------------------------ ###

### Rearrange =====================================================
## CCE
cceDet <- read.csv(paste0(datap,"CCE_regress_Det_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
cceDet <- subset.data.frame(cceDet,Type %in% c("D","A","P","F"))
cceDet$Lag <- as.factor(cceDet$Lag)
cceDet$Type <- factor(cceDet$Type, 
                         levels = c("D","A","P","F"))
rcceDet <- cceDet[,c("Lag","Type","coef","p")]
rcceDet$p[rcceDet$p > 0.1] <- NA
rcceDet$sym <- rcceDet$p
rcceDet$sym[rcceDet$p <= 0.05] <- "*"
rcceDet$sym[rcceDet$p > 0.05] <- "."
rcceDet$sym[rcceDet$p > 0.1] <- NA

cceTB <- read.csv(paste0(datap,"CCE_regress_TB_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
cceTB <- subset.data.frame(cceTB,Type %in% c("D","A","P","F"))
cceTB$Lag <- as.factor(cceTB$Lag)
cceTB$Type <- factor(cceTB$Type, 
                      levels = c("D","A","P","F"))
rcceTB <- cceTB[,c("Lag","Type","coef","p")]
rcceTB$p[rcceTB$p > 0.1] <- NA
rcceTB$sym <- rcceTB$p
rcceTB$sym[rcceTB$p <= 0.05] <- "*"
rcceTB$sym[rcceTB$p > 0.05] <- "."
rcceTB$sym[rcceTB$p > 0.1] <- NA

cceTP <- read.csv(paste0(datap,"CCE_regress_TP_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
cceTP <- subset.data.frame(cceTP,Type %in% c("D","A","P","F"))
cceTP$Lag <- as.factor(cceTP$Lag)
cceTP$Type <- factor(cceTP$Type, 
                      levels = c("D","A","P","F"))
rcceTP <- cceTP[,c("Lag","Type","coef","p")]
rcceTP$p[rcceTP$p > 0.1] <- NA
rcceTP$sym <- rcceTP$p
rcceTP$sym[rcceTP$p <= 0.05] <- "*"
rcceTP$sym[rcceTP$p > 0.05] <- "."
rcceTP$sym[rcceTP$p > 0.1] <- NA

cceZmeso <- read.csv(paste0(datap,"CCE_regress_Zmeso_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
cceZmeso <- subset.data.frame(cceZmeso,Type %in% c("D","A","P","F"))
cceZmeso$Lag <- as.factor(cceZmeso$Lag)
cceZmeso$Type <- factor(cceZmeso$Type, 
                      levels = c("D","A","P","F"))
rcceZmeso <- cceZmeso[,c("Lag","Type","coef","p")]
rcceZmeso$p[rcceZmeso$p > 0.1] <- NA
rcceZmeso$sym <- rcceZmeso$p
rcceZmeso$sym[rcceZmeso$p <= 0.05] <- "*"
rcceZmeso$sym[rcceZmeso$p > 0.05] <- "."
rcceZmeso$sym[rcceZmeso$p > 0.1] <- NA

cceZmLoss <- read.csv(paste0(datap,"CCE_regress_ZmLoss_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
cceZmLoss <- subset.data.frame(cceZmLoss,Type %in% c("D","A","P","F"))
cceZmLoss$Lag <- as.factor(cceZmLoss$Lag)
cceZmLoss$Type <- factor(cceZmLoss$Type, 
                      levels = c("D","A","P","F"))
rcceZmLoss <- cceZmLoss[,c("Lag","Type","coef","p")]
rcceZmLoss$p[rcceZmLoss$p > 0.1] <- NA
rcceZmLoss$sym <- rcceZmLoss$p
rcceZmLoss$sym[rcceZmLoss$p <= 0.05] <- "*"
rcceZmLoss$sym[rcceZmLoss$p > 0.05] <- "."
rcceZmLoss$sym[rcceZmLoss$p > 0.1] <- NA


## EBS
ebsDet <- read.csv(paste0(datap,"EBS_regress_Det_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
ebsDet <- subset.data.frame(ebsDet,Type %in% c("D","A","P","F"))
ebsDet$Lag <- as.factor(ebsDet$Lag)
ebsDet$Type <- factor(ebsDet$Type, 
                      levels = c("D","A","P","F"))
rebsDet <- ebsDet[,c("Lag","Type","coef","p")]
rebsDet$p[rebsDet$p > 0.1] <- NA
rebsDet$sym <- rebsDet$p
rebsDet$sym[rebsDet$p <= 0.05] <- "*"
rebsDet$sym[rebsDet$p > 0.05] <- "."
rebsDet$sym[rebsDet$p > 0.1] <- NA

ebsTB <- read.csv(paste0(datap,"EBS_regress_TB_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
ebsTB <- subset.data.frame(ebsTB,Type %in% c("D","A","P","F"))
ebsTB$Lag <- as.factor(ebsTB$Lag)
ebsTB$Type <- factor(ebsTB$Type, 
                      levels = c("D","A","P","F"))
rebsTB <- ebsTB[,c("Lag","Type","coef","p")]
rebsTB$p[rebsTB$p > 0.1] <- NA
rebsTB$sym <- rebsTB$p
rebsTB$sym[rebsTB$p <= 0.05] <- "*"
rebsTB$sym[rebsTB$p > 0.05] <- "."
rebsTB$sym[rebsTB$p > 0.1] <- NA

ebsTP <- read.csv(paste0(datap,"EBS_regress_TP_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
ebsTP <- subset.data.frame(ebsTP,Type %in% c("D","A","P","F"))
ebsTP$Lag <- as.factor(ebsTP$Lag)
ebsTP$Type <- factor(ebsTP$Type, 
                      levels = c("D","A","P","F"))
rebsTP <- ebsTP[,c("Lag","Type","coef","p")]
rebsTP$p[rebsTP$p > 0.1] <- NA
rebsTP$sym <- rebsTP$p
rebsTP$sym[rebsTP$p <= 0.05] <- "*"
rebsTP$sym[rebsTP$p > 0.05] <- "."
rebsTP$sym[rebsTP$p > 0.1] <- NA

ebsZmeso <- read.csv(paste0(datap,"EBS_regress_Zmeso_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
ebsZmeso <- subset.data.frame(ebsZmeso,Type %in% c("D","A","P","F"))
ebsZmeso$Lag <- as.factor(ebsZmeso$Lag)
ebsZmeso$Type <- factor(ebsZmeso$Type, 
                      levels = c("D","A","P","F"))
rebsZmeso <- ebsZmeso[,c("Lag","Type","coef","p")]
rebsZmeso$p[rebsZmeso$p > 0.1] <- NA
rebsZmeso$sym <- rebsZmeso$p
rebsZmeso$sym[rebsZmeso$p <= 0.05] <- "*"
rebsZmeso$sym[rebsZmeso$p > 0.05] <- "."
rebsZmeso$sym[rebsZmeso$p > 0.1] <- NA

ebsZmLoss <- read.csv(paste0(datap,"EBS_regress_ZmLoss_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
ebsZmLoss <- subset.data.frame(ebsZmLoss,Type %in% c("D","A","P","F"))
ebsZmLoss$Lag <- as.factor(ebsZmLoss$Lag)
ebsZmLoss$Type <- factor(ebsZmLoss$Type, 
                      levels = c("D","A","P","F"))
rebsZmLoss <- ebsZmLoss[,c("Lag","Type","coef","p")]
rebsZmLoss$p[rebsZmLoss$p > 0.1] <- NA
rebsZmLoss$sym <- rebsZmLoss$p
rebsZmLoss$sym[rebsZmLoss$p <= 0.05] <- "*"
rebsZmLoss$sym[rebsZmLoss$p > 0.05] <- "."
rebsZmLoss$sym[rebsZmLoss$p > 0.1] <- NA


## AI
aiDet <- read.csv(paste0(datap,"AI_regress_Det_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
aiDet <- subset.data.frame(aiDet,Type %in% c("D","A","P","F"))
aiDet$Lag <- as.factor(aiDet$Lag)
aiDet$Type[aiDet$Type=="Zlos"] <- "ZmLoss"
aiDet$Type[aiDet$Type=="Zoo"] <- "Zmeso"
aiDet$Type <- factor(aiDet$Type, 
                     levels = c("D","A","P","F"))
raiDet <- aiDet[,c("Lag","Type","coef","p")]
raiDet$p[raiDet$p > 0.1] <- NA
raiDet$sym <- raiDet$p
raiDet$sym[raiDet$p <= 0.05] <- "*"
raiDet$sym[raiDet$p > 0.05] <- "."
raiDet$sym[raiDet$p > 0.1] <- NA

aiTB <- read.csv(paste0(datap,"AI_regress_TB_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
aiTB <- subset.data.frame(aiTB,Type %in% c("D","A","P","F"))
aiTB$Lag <- as.factor(aiTB$Lag)
aiTB$Type[aiTB$Type=="Zlos"] <- "ZmLoss"
aiTB$Type[aiTB$Type=="Zoo"] <- "Zmeso"
aiTB$Type <- factor(aiTB$Type, 
                        levels = c("D","A","P","F"))
raiTB <- aiTB[,c("Lag","Type","coef","p")]
raiTB$p[raiTB$p > 0.1] <- NA
raiTB$sym <- raiTB$p
raiTB$sym[raiTB$p <= 0.05] <- "*"
raiTB$sym[raiTB$p > 0.05] <- "."
raiTB$sym[raiTB$p > 0.1] <- NA

aiTP <- read.csv(paste0(datap,"AI_regress_TP_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
aiTP <- subset.data.frame(aiTP,Type %in% c("D","A","P","F"))
aiTP$Lag <- as.factor(aiTP$Lag)
aiTP$Type <- factor(aiTP$Type, 
                     levels = c("D","A","P","F"))
raiTP <- aiTP[,c("Lag","Type","coef","p")]
raiTP$p[raiTP$p > 0.1] <- NA
raiTP$sym <- raiTP$p
raiTP$sym[raiTP$p <= 0.05] <- "*"
raiTP$sym[raiTP$p > 0.05] <- "."
raiTP$sym[raiTP$p > 0.1] <- NA

aiZmeso <- read.csv(paste0(datap,"AI_regress_Zmeso_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
aiZmeso <- subset.data.frame(aiZmeso,Type %in% c("D","A","P","F"))
aiZmeso$Lag <- as.factor(aiZmeso$Lag)
aiZmeso$Type <- factor(aiZmeso$Type, 
                        levels = c("D","A","P","F"))
raiZmeso <- aiZmeso[,c("Lag","Type","coef","p")]
raiZmeso$p[raiZmeso$p > 0.1] <- NA
raiZmeso$sym <- raiZmeso$p
raiZmeso$sym[raiZmeso$p <= 0.05] <- "*"
raiZmeso$sym[raiZmeso$p > 0.05] <- "."
raiZmeso$sym[raiZmeso$p > 0.1] <- NA

aiZmLoss <- read.csv(paste0(datap,"AI_regress_ZmLoss_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
aiZmLoss <- subset.data.frame(aiZmLoss,Type %in% c("D","A","P","F"))
aiZmLoss$Lag <- as.factor(aiZmLoss$Lag)
aiZmLoss$Type[aiZmLoss$Type=="Zlos"] <- "ZmLoss"
aiZmLoss$Type[aiZmLoss$Type=="Zoo"] <- "Zmeso"
aiZmLoss$Type <- factor(aiZmLoss$Type, 
                    levels = c("D","A","P","F"))
raiZmLoss <- aiZmLoss[,c("Lag","Type","coef","p")]
raiZmLoss$p[raiZmLoss$p > 0.1] <- NA
raiZmLoss$sym <- raiZmLoss$p
raiZmLoss$sym[raiZmLoss$p <= 0.05] <- "*"
raiZmLoss$sym[raiZmLoss$p > 0.05] <- "."
raiZmLoss$sym[raiZmLoss$p > 0.1] <- NA



## GAK
gakDet <- read.csv(paste0(datap,"GAK_regress_Det_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
gakDet <- subset.data.frame(gakDet,Type %in% c("D","A","P","F"))
gakDet$Lag <- as.factor(gakDet$Lag)
gakDet$Type <- factor(gakDet$Type, 
                      levels = c("D","A","P","F"))
rgakDet <- gakDet[,c("Lag","Type","coef","p")]
rgakDet$p[rgakDet$p > 0.1] <- NA
rgakDet$sym <- rgakDet$p
rgakDet$sym[rgakDet$p <= 0.05] <- "*"
rgakDet$sym[rgakDet$p > 0.05] <- "."
rgakDet$sym[rgakDet$p > 0.1] <- NA

gakTB <- read.csv(paste0(datap,"GAK_regress_TB_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
gakTB <- subset.data.frame(gakTB,Type %in% c("D","A","P","F"))
gakTB$Lag <- as.factor(gakTB$Lag)
gakTB$Type <- factor(gakTB$Type, 
                      levels = c("D","A","P","F"))
rgakTB <- gakTB[,c("Lag","Type","coef","p")]
rgakTB$p[rgakTB$p > 0.1] <- NA
rgakTB$sym <- rgakTB$p
rgakTB$sym[rgakTB$p <= 0.05] <- "*"
rgakTB$sym[rgakTB$p > 0.05] <- "."
rgakTB$sym[rgakTB$p > 0.1] <- NA

gakTP <- read.csv(paste0(datap,"GAK_regress_TP_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
gakTP <- subset.data.frame(gakTP,Type %in% c("D","A","P","F"))
gakTP$Lag <- as.factor(gakTP$Lag)
gakTP$Type <- factor(gakTP$Type, 
                      levels = c("D","A","P","F"))
rgakTP <- gakTP[,c("Lag","Type","coef","p")]
rgakTP$p[rgakTP$p > 0.1] <- NA
rgakTP$sym <- rgakTP$p
rgakTP$sym[rgakTP$p <= 0.05] <- "*"
rgakTP$sym[rgakTP$p > 0.05] <- "."
rgakTP$sym[rgakTP$p > 0.1] <- NA

gakZmeso <- read.csv(paste0(datap,"GAK_regress_Zmeso_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
gakZmeso <- subset.data.frame(gakZmeso,Type %in% c("D","A","P","F"))
gakZmeso$Lag <- as.factor(gakZmeso$Lag)
gakZmeso$Type <- factor(gakZmeso$Type, 
                      levels = c("D","A","P","F"))
rgakZmeso <- gakZmeso[,c("Lag","Type","coef","p")]
rgakZmeso$p[rgakZmeso$p > 0.1] <- NA
rgakZmeso$sym <- rgakZmeso$p
rgakZmeso$sym[rgakZmeso$p <= 0.05] <- "*"
rgakZmeso$sym[rgakZmeso$p > 0.05] <- "."
rgakZmeso$sym[rgakZmeso$p > 0.1] <- NA

gakZmLoss <- read.csv(paste0(datap,"GAK_regress_ZmLoss_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
gakZmLoss <- subset.data.frame(gakZmLoss,Type %in% c("D","A","P","F"))
gakZmLoss$Lag <- as.factor(gakZmLoss$Lag)
gakZmLoss$Type <- factor(gakZmLoss$Type, 
                      levels = c("D","A","P","F"))
rgakZmLoss <- gakZmLoss[,c("Lag","Type","coef","p")]
rgakZmLoss$p[rgakZmLoss$p > 0.1] <- NA
rgakZmLoss$sym <- rgakZmLoss$p
rgakZmLoss$sym[rgakZmLoss$p <= 0.05] <- "*"
rgakZmLoss$sym[rgakZmLoss$p > 0.05] <- "."
rgakZmLoss$sym[rgakZmLoss$p > 0.1] <- NA


## HI
hiDet <- read.csv(paste0(datap,"HI_regress_Det_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
hiDet <- subset.data.frame(hiDet,Type %in% c("D","A","P","F"))
hiDet$Lag <- as.factor(hiDet$Lag)
hiDet$Type <- factor(hiDet$Type, 
                      levels = c("D","A","P","F"))
rhiDet <- hiDet[,c("Lag","Type","coef","p")]
rhiDet$p[rhiDet$p > 0.1] <- NA
rhiDet$sym <- rhiDet$p
rhiDet$sym[rhiDet$p <= 0.05] <- "*"
rhiDet$sym[rhiDet$p > 0.05] <- "."
rhiDet$sym[rhiDet$p > 0.1] <- NA

hiTB <- read.csv(paste0(datap,"HI_regress_TB_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
hiTB <- subset.data.frame(hiTB,Type %in% c("D","A","P","F"))
hiTB$Lag <- as.factor(hiTB$Lag)
hiTB$Type <- factor(hiTB$Type, 
                      levels = c("D","A","P","F"))
rhiTB <- hiTB[,c("Lag","Type","coef","p")]
rhiTB$p[rhiTB$p > 0.1] <- NA
rhiTB$sym <- rhiTB$p
rhiTB$sym[rhiTB$p <= 0.05] <- "*"
rhiTB$sym[rhiTB$p > 0.05] <- "."
rhiTB$sym[rhiTB$p > 0.1] <- NA

hiTP <- read.csv(paste0(datap,"HI_regress_TP_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
hiTP <- subset.data.frame(hiTP,Type %in% c("D","A","P","F"))
hiTP$Lag <- as.factor(hiTP$Lag)
hiTP$Type <- factor(hiTP$Type, 
                      levels = c("D","A","P","F"))
rhiTP <- hiTP[,c("Lag","Type","coef","p")]
rhiTP$p[rhiTP$p > 0.1] <- NA
rhiTP$sym <- rhiTP$p
rhiTP$sym[rhiTP$p <= 0.05] <- "*"
rhiTP$sym[rhiTP$p > 0.05] <- "."
rhiTP$sym[rhiTP$p > 0.1] <- NA

hiZmeso <- read.csv(paste0(datap,"HI_regress_Zmeso_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
hiZmeso <- subset.data.frame(hiZmeso,Type %in% c("D","A","P","F"))
hiZmeso$Lag <- as.factor(hiZmeso$Lag)
hiZmeso$Type <- factor(hiZmeso$Type, 
                      levels = c("D","A","P","F"))
rhiZmeso <- hiZmeso[,c("Lag","Type","coef","p")]
rhiZmeso$p[rhiZmeso$p > 0.1] <- NA
rhiZmeso$sym <- rhiZmeso$p
rhiZmeso$sym[rhiZmeso$p <= 0.05] <- "*"
rhiZmeso$sym[rhiZmeso$p > 0.05] <- "."
rhiZmeso$sym[rhiZmeso$p > 0.1] <- NA

hiZmLoss <- read.csv(paste0(datap,"HI_regress_ZmLoss_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
hiZmLoss <- subset.data.frame(hiZmLoss,Type %in% c("D","A","P","F"))
hiZmLoss$Lag <- as.factor(hiZmLoss$Lag)
hiZmLoss$Type <- factor(hiZmLoss$Type, 
                      levels = c("D","A","P","F"))
rhiZmLoss <- hiZmLoss[,c("Lag","Type","coef","p")]
rhiZmLoss$p[rhiZmLoss$p > 0.1] <- NA
rhiZmLoss$sym <- rhiZmLoss$p
rhiZmLoss$sym[rhiZmLoss$p <= 0.05] <- "*"
rhiZmLoss$sym[rhiZmLoss$p > 0.05] <- "."
rhiZmLoss$sym[rhiZmLoss$p > 0.1] <- NA


## CHK
chkDet <- read.csv(paste0(datap,"CHK_regress_Det_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
chkDet <- subset.data.frame(chkDet,Type %in% c("D","A","P","F"))
chkDet$Lag <- as.factor(chkDet$Lag)
chkDet$Type <- factor(chkDet$Type, 
                      levels = c("D","A","P","F"))
rchkDet <- chkDet[,c("Lag","Type","coef","p")]
rchkDet$p[rchkDet$p > 0.1] <- NA
rchkDet$sym <- rchkDet$p
rchkDet$sym[rchkDet$p <= 0.05] <- "*"
rchkDet$sym[rchkDet$p > 0.05] <- "."
rchkDet$sym[rchkDet$p > 0.1] <- NA

chkTB <- read.csv(paste0(datap,"CHK_regress_TB_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
chkTB <- subset.data.frame(chkTB,Type %in% c("D","A","P","F"))
chkTB$Lag <- as.factor(chkTB$Lag)
chkTB$Type <- factor(chkTB$Type, 
                      levels = c("D","A","P","F"))
rchkTB <- chkTB[,c("Lag","Type","coef","p")]
rchkTB$p[rchkTB$p > 0.1] <- NA
rchkTB$sym <- rchkTB$p
rchkTB$sym[rchkTB$p <= 0.05] <- "*"
rchkTB$sym[rchkTB$p > 0.05] <- "."
rchkTB$sym[rchkTB$p > 0.1] <- NA

chkTP <- read.csv(paste0(datap,"CHK_regress_TP_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
chkTP <- subset.data.frame(chkTP,Type %in% c("D","A","P","F"))
chkTP$Lag <- as.factor(chkTP$Lag)
chkTP$Type <- factor(chkTP$Type, 
                      levels = c("D","A","P","F"))
rchkTP <- chkTP[,c("Lag","Type","coef","p")]
rchkTP$p[rchkTP$p > 0.1] <- NA
rchkTP$sym <- rchkTP$p
rchkTP$sym[rchkTP$p <= 0.05] <- "*"
rchkTP$sym[rchkTP$p > 0.05] <- "."
rchkTP$sym[rchkTP$p > 0.1] <- NA

chkZmeso <- read.csv(paste0(datap,"CHK_regress_Zmeso_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
chkZmeso <- subset.data.frame(chkZmeso,Type %in% c("D","A","P","F"))
chkZmeso$Lag <- as.factor(chkZmeso$Lag)
chkZmeso$Type <- factor(chkZmeso$Type, 
                      levels = c("D","A","P","F"))
rchkZmeso <- chkZmeso[,c("Lag","Type","coef","p")]
rchkZmeso$p[rchkZmeso$p > 0.1] <- NA
rchkZmeso$sym <- rchkZmeso$p
rchkZmeso$sym[rchkZmeso$p <= 0.05] <- "*"
rchkZmeso$sym[rchkZmeso$p > 0.05] <- "."
rchkZmeso$sym[rchkZmeso$p > 0.1] <- NA

chkZmLoss <- read.csv(paste0(datap,"CHK_regress_ZmLoss_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
chkZmLoss <- subset.data.frame(chkZmLoss,Type %in% c("D","A","P","F"))
chkZmLoss$Lag <- as.factor(chkZmLoss$Lag)
chkZmLoss$Type <- factor(chkZmLoss$Type, 
                      levels = c("D","A","P","F"))
rchkZmLoss <- chkZmLoss[,c("Lag","Type","coef","p")]
rchkZmLoss$p[rchkZmLoss$p > 0.1] <- NA
rchkZmLoss$sym <- rchkZmLoss$p
rchkZmLoss$sym[rchkZmLoss$p <= 0.05] <- "*"
rchkZmLoss$sym[rchkZmLoss$p > 0.05] <- "."
rchkZmLoss$sym[rchkZmLoss$p > 0.1] <- NA


##GMX
gmxDet <- read.csv(paste0(datap,"GMX_regress_Det_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
gmxDet <- subset.data.frame(gmxDet,Type %in% c("D","A","P","F"))
gmxDet$Lag <- as.factor(gmxDet$Lag)
gmxDet$Type <- factor(gmxDet$Type, 
                      levels = c("D","A","P","F"))
rgmxDet <- gmxDet[,c("Lag","Type","coef","p")]
rgmxDet$p[rgmxDet$p > 0.1] <- NA
rgmxDet$sym <- rgmxDet$p
rgmxDet$sym[rgmxDet$p <= 0.05] <- "*"
rgmxDet$sym[rgmxDet$p > 0.05] <- "."
rgmxDet$sym[rgmxDet$p > 0.1] <- NA

gmxTB <- read.csv(paste0(datap,"GMX_regress_TB_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
gmxTB <- subset.data.frame(gmxTB,Type %in% c("D","A","P","F"))
gmxTB$Lag <- as.factor(gmxTB$Lag)
gmxTB$Type <- factor(gmxTB$Type, 
                      levels = c("D","A","P","F"))
rgmxTB <- gmxTB[,c("Lag","Type","coef","p")]
rgmxTB$p[rgmxTB$p > 0.1] <- NA
rgmxTB$sym <- rgmxTB$p
rgmxTB$sym[rgmxTB$p <= 0.05] <- "*"
rgmxTB$sym[rgmxTB$p > 0.05] <- "."
rgmxTB$sym[rgmxTB$p > 0.1] <- NA

gmxTP <- read.csv(paste0(datap,"GMX_regress_TP_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
gmxTP <- subset.data.frame(gmxTP,Type %in% c("D","A","P","F"))
gmxTP$Lag <- as.factor(gmxTP$Lag)
gmxTP$Type <- factor(gmxTP$Type, 
                      levels = c("D","A","P","F"))
rgmxTP <- gmxTP[,c("Lag","Type","coef","p")]
rgmxTP$p[rgmxTP$p > 0.1] <- NA
rgmxTP$sym <- rgmxTP$p
rgmxTP$sym[rgmxTP$p <= 0.05] <- "*"
rgmxTP$sym[rgmxTP$p > 0.05] <- "."
rgmxTP$sym[rgmxTP$p > 0.1] <- NA

gmxZmeso <- read.csv(paste0(datap,"GMX_regress_Zmeso_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
gmxZmeso <- subset.data.frame(gmxZmeso,Type %in% c("D","A","P","F"))
gmxZmeso$Lag <- as.factor(gmxZmeso$Lag)
gmxZmeso$Type <- factor(gmxZmeso$Type, 
                      levels = c("D","A","P","F"))
rgmxZmeso <- gmxZmeso[,c("Lag","Type","coef","p")]
rgmxZmeso$p[rgmxZmeso$p > 0.1] <- NA
rgmxZmeso$sym <- rgmxZmeso$p
rgmxZmeso$sym[rgmxZmeso$p <= 0.05] <- "*"
rgmxZmeso$sym[rgmxZmeso$p > 0.05] <- "."
rgmxZmeso$sym[rgmxZmeso$p > 0.1] <- NA

gmxZmLoss <- read.csv(paste0(datap,"GMX_regress_ZmLoss_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
gmxZmLoss <- subset.data.frame(gmxZmLoss,Type %in% c("D","A","P","F"))
gmxZmLoss$Lag <- as.factor(gmxZmLoss$Lag)
gmxZmLoss$Type <- factor(gmxZmLoss$Type, 
                      levels = c("D","A","P","F"))
rgmxZmLoss <- gmxZmLoss[,c("Lag","Type","coef","p")]
rgmxZmLoss$p[rgmxZmLoss$p > 0.1] <- NA
rgmxZmLoss$sym <- rgmxZmLoss$p
rgmxZmLoss$sym[rgmxZmLoss$p <= 0.05] <- "*"
rgmxZmLoss$sym[rgmxZmLoss$p > 0.05] <- "."
rgmxZmLoss$sym[rgmxZmLoss$p > 0.1] <- NA


## NE
neDet <- read.csv(paste0(datap,"NE_regress_Det_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
neDet <- subset.data.frame(neDet,Type %in% c("D","A","P","F"))
neDet$Lag <- as.factor(neDet$Lag)
neDet$Type <- factor(neDet$Type, 
                      levels = c("D","A","P","F"))
rneDet <- neDet[,c("Lag","Type","coef","p")]
rneDet$p[rneDet$p > 0.1] <- NA
rneDet$sym <- rneDet$p
rneDet$sym[rneDet$p <= 0.05] <- "*"
rneDet$sym[rneDet$p > 0.05] <- "."
rneDet$sym[rneDet$p > 0.1] <- NA

neTB <- read.csv(paste0(datap,"NE_regress_TB_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
neTB <- subset.data.frame(neTB,Type %in% c("D","A","P","F"))
neTB$Lag <- as.factor(neTB$Lag)
neTB$Type <- factor(neTB$Type, 
                      levels = c("D","A","P","F"))
rneTB <- neTB[,c("Lag","Type","coef","p")]
rneTB$p[rneTB$p > 0.1] <- NA
rneTB$sym <- rneTB$p
rneTB$sym[rneTB$p <= 0.05] <- "*"
rneTB$sym[rneTB$p > 0.05] <- "."
rneTB$sym[rneTB$p > 0.1] <- NA

neTP <- read.csv(paste0(datap,"NE_regress_TP_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
neTP <- subset.data.frame(neTP,Type %in% c("D","A","P","F"))
neTP$Lag <- as.factor(neTP$Lag)
neTP$Type <- factor(neTP$Type, 
                      levels = c("D","A","P","F"))
rneTP <- neTP[,c("Lag","Type","coef","p")]
rneTP$p[rneTP$p > 0.1] <- NA
rneTP$sym <- rneTP$p
rneTP$sym[rneTP$p <= 0.05] <- "*"
rneTP$sym[rneTP$p > 0.05] <- "."
rneTP$sym[rneTP$p > 0.1] <- NA

neZmeso <- read.csv(paste0(datap,"NE_regress_Zmeso_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
neZmeso <- subset.data.frame(neZmeso,Type %in% c("D","A","P","F"))
neZmeso$Lag <- as.factor(neZmeso$Lag)
neZmeso$Type <- factor(neZmeso$Type, 
                      levels = c("D","A","P","F"))
rneZmeso <- neZmeso[,c("Lag","Type","coef","p")]
rneZmeso$p[rneZmeso$p > 0.1] <- NA
rneZmeso$sym <- rneZmeso$p
rneZmeso$sym[rneZmeso$p <= 0.05] <- "*"
rneZmeso$sym[rneZmeso$p > 0.05] <- "."
rneZmeso$sym[rneZmeso$p > 0.1] <- NA

neZmLoss <- read.csv(paste0(datap,"NE_regress_ZmLoss_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
neZmLoss <- subset.data.frame(neZmLoss,Type %in% c("D","A","P","F"))
neZmLoss$Lag <- as.factor(neZmLoss$Lag)
neZmLoss$Type <- factor(neZmLoss$Type, 
                      levels = c("D","A","P","F"))
rneZmLoss <- neZmLoss[,c("Lag","Type","coef","p")]
rneZmLoss$p[rneZmLoss$p > 0.1] <- NA
rneZmLoss$sym <- rneZmLoss$p
rneZmLoss$sym[rneZmLoss$p <= 0.05] <- "*"
rneZmLoss$sym[rneZmLoss$p > 0.05] <- "."
rneZmLoss$sym[rneZmLoss$p > 0.1] <- NA


##SE
seDet <- read.csv(paste0(datap,"SE_regress_Det_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
seDet <- subset.data.frame(seDet,Type %in% c("D","A","P","F"))
seDet$Lag <- as.factor(seDet$Lag)
seDet$Type <- factor(seDet$Type, 
                      levels = c("D","A","P","F"))
rseDet <- seDet[,c("Lag","Type","coef","p")]
rseDet$p[rseDet$p > 0.1] <- NA
rseDet$sym <- rseDet$p
rseDet$sym[rseDet$p <= 0.05] <- "*"
rseDet$sym[rseDet$p > 0.05] <- "."
rseDet$sym[rseDet$p > 0.1] <- NA

seTB <- read.csv(paste0(datap,"SE_regress_TB_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
seTB <- subset.data.frame(seTB,Type %in% c("D","A","P","F"))
seTB$Lag <- as.factor(seTB$Lag)
seTB$Type <- factor(seTB$Type, 
                      levels = c("D","A","P","F"))
rseTB <- seTB[,c("Lag","Type","coef","p")]
rseTB$p[rseTB$p > 0.1] <- NA
rseTB$sym <- rseTB$p
rseTB$sym[rseTB$p <= 0.05] <- "*"
rseTB$sym[rseTB$p > 0.05] <- "."
rseTB$sym[rseTB$p > 0.1] <- NA

seTP <- read.csv(paste0(datap,"SE_regress_TP_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
seTP <- subset.data.frame(seTP,Type %in% c("D","A","P","F"))
seTP$Lag <- as.factor(seTP$Lag)
seTP$Type <- factor(seTP$Type, 
                      levels = c("D","A","P","F"))
rseTP <- seTP[,c("Lag","Type","coef","p")]
rseTP$p[rseTP$p > 0.1] <- NA
rseTP$sym <- rseTP$p
rseTP$sym[rseTP$p <= 0.05] <- "*"
rseTP$sym[rseTP$p > 0.05] <- "."
rseTP$sym[rseTP$p > 0.1] <- NA

seZmeso <- read.csv(paste0(datap,"SE_regress_Zmeso_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
seZmeso <- subset.data.frame(seZmeso,Type %in% c("D","A","P","F"))
seZmeso$Lag <- as.factor(seZmeso$Lag)
seZmeso$Type <- factor(seZmeso$Type, 
                      levels = c("D","A","P","F"))
rseZmeso <- seZmeso[,c("Lag","Type","coef","p")]
rseZmeso$p[rseZmeso$p > 0.1] <- NA
rseZmeso$sym <- rseZmeso$p
rseZmeso$sym[rseZmeso$p <= 0.05] <- "*"
rseZmeso$sym[rseZmeso$p > 0.05] <- "."
rseZmeso$sym[rseZmeso$p > 0.1] <- NA

seZmLoss <- read.csv(paste0(datap,"SE_regress_ZmLoss_div2SD_melt_mat_nu.csv"),sep=",",header = T,stringsAsFactors = F)
seZmLoss <- subset.data.frame(seZmLoss,Type %in% c("D","A","P","F"))
seZmLoss$Lag <- as.factor(seZmLoss$Lag)
seZmLoss$Type <- factor(seZmLoss$Type, 
                      levels = c("D","A","P","F"))
rseZmLoss <- seZmLoss[,c("Lag","Type","coef","p")]
rseZmLoss$p[rseZmLoss$p > 0.1] <- NA
rseZmLoss$sym <- rseZmLoss$p
rseZmLoss$sym[rseZmLoss$p <= 0.05] <- "*"
rseZmLoss$sym[rseZmLoss$p > 0.05] <- "."
rseZmLoss$sym[rseZmLoss$p > 0.1] <- NA


###====== Plots ===================================================================
## CCE
c5 <- ggplot(data = rcceDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Det\nCoef") +
  theme_minimal() + labs(x="")+ labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Det") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

c4 <- ggplot(data = rcceTB, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TB\nCoef") +
  theme_minimal() + labs(x="")+ labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Tb") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

c1 <- ggplot(data = rcceTP, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),  name="TP\nCoef") +
  theme_minimal() + labs(x="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("CCE Tp") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

c2 <- ggplot(data = rcceZmeso, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Zmeso\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Zmeso") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


c3 <- ggplot(data = rcceZmLoss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="ZmLoss\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("ZmLoss") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


## EBS
e5 <- ggplot(data = rebsDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Det\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

e4 <- ggplot(data = rebsTB, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TB\nCoef") +
  theme_minimal() + labs(x="")+ labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

e1 <- ggplot(data = rebsTP, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TP\nCoef") +
  theme_minimal() + labs(x="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("EBS") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

e2 <- ggplot(data = rebsZmeso, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +labs(y="")+
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Zmeso\nCoef") +
  theme_minimal() + ggtitle("") + labs(x="")+ labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

e3 <- ggplot(data = rebsZmLoss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="ZmLoss\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


## CHK
k5 <- ggplot(data = rchkDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.15,1.15),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Det\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

k4 <- ggplot(data = rchkTB, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TB\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

k1 <- ggplot(data = rchkTP, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),  name="TP\nCoef") +
  theme_minimal() + ggtitle("CHK") + labs(x="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

k2 <- ggplot(data = rchkZmeso, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Zmeso\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

k3 <- ggplot(data = rchkZmLoss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="ZmLoss\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


## NE
n5 <- ggplot(data = rneDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Det\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

n4 <- ggplot(data = rneTB, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TB\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

n1 <- ggplot(data = rneTP, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TP\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("NEUS") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

n2 <- ggplot(data = rneZmeso, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.71,1.71),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Zmeso\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

n3 <- ggplot(data = rneZmLoss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Coef") +
  theme_minimal() + labs(x="Lag (y)")+labs(y="")+
  #theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

n33 <- ggplot(data = rneZmLoss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Coef") +
  theme_minimal() + labs(x="Lag (y)")+labs(y="")+
  theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  #guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 



png(paste0(figp,'Heatmaps_CCE_EBS_GAK_NE_coef_inputs_fntype2_v2_nu.png'), 
    width = 7.5*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
ggarrange(c1,c2,c3,c4,c5,
          e1,e2,e3,e4,e5,
          k1,k2,k3,k4,k5,
          n1,n2,n33,n4,n5,
          nrow = 4, ncol = 5,labels = "auto",
          common.legend = TRUE, legend = "right")
dev.off()



## AI
i5 <- ggplot(data = raiDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Det\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  #theme(legend.key.height=unit(0.3,'cm'),legend.key.width=unit(0.1, 'cm'))+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

i4 <- ggplot(data = raiTB, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TB\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

i1 <- ggplot(data = raiTP, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TP\nCoef") +
  theme_minimal() + labs(x="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("AI") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

i2 <- ggplot(data = raiZmeso, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Zmeso\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

i3 <- ggplot(data = raiZmLoss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="ZmLoss\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


## GAK
g5 <- ggplot(data = rgakDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Det\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Det") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

g4 <- ggplot(data = rgakTB, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TB\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Tb") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

g1 <- ggplot(data = rgakTP, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TP\nCoef") +
  theme_minimal() + labs(x="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("GAK Tp") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

g2 <- ggplot(data = rgakZmeso, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Zmeso\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("Zmeso") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

g3 <- ggplot(data = rgakZmLoss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),  name="ZmLoss\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("ZmLoss") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


## HI
h5 <- ggplot(data = rhiDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Det\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

h4 <- ggplot(data = rhiTB, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Coef") +
  theme_minimal() + labs(x="")+labs(y="")+
  theme(legend.key.height=unit(0.4,'cm'),legend.key.width=unit(0.15, 'cm'))+
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

h1 <- ggplot(data = rhiTP, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TP\nCoef") +
  theme_minimal() + labs(x="")+ 
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("HI") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

h2 <- ggplot(data = rhiZmeso, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Zmeso\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

h3 <- ggplot(data = rhiZmLoss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="ZmLoss\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

## GMX
m5 <- ggplot(data = rgmxDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Det\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

m4 <- ggplot(data = rgmxTB, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TB\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

m1 <- ggplot(data = rgmxTP, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TP\nCoef") +
  theme_minimal() + ggtitle("GMX") + labs(x="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

m2 <- ggplot(data = rgmxZmeso, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Zmeso\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

m3 <- ggplot(data = rgmxZmLoss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="ZmLoss\nCoef") +
  theme_minimal() + labs(x="")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 



## SE
s5 <- ggplot(data = rseDet, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Det\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

s4 <- ggplot(data = rseTB, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TB\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

s1 <- ggplot(data = rseTP, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="TP\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("SEUS") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

s2 <- ggplot(data = rseZmeso, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="Zmeso\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

s3 <- ggplot(data = rseZmLoss, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.5,0.5),
                       breaks=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5),
                       labels=c(-0.5,-0.2,-0.1,0,0.1,0.2,0.5), name="ZmLoss\nCoef") +
  theme_minimal() + labs(x="Lag (y)")+labs(y="")+
  guides(fill = "none")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 


### Put some all together for SUpp

png(paste0(figp,'Heatmaps_AI_GAK_HI_GMX_SE_coef_inputs_fntypes_v2_nu.png'), 
    width = 7.5*300,        # 5 x 300 pixels
    height = 7.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
ggarrange(g1,g2,g3,g4,g5,
          i1,i2,i3,i4,i5,
          h1,h2,h3,h4,h5,
          m1,m2,m3,m4,m5,
          s1,s2,s3,s4,s5,
          nrow = 5, ncol = 5,labels = "auto",
          common.legend = TRUE, legend = "right")
dev.off()




