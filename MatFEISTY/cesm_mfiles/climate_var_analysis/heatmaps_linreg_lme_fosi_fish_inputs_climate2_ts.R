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
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"

### ------------------------------------------------------------------------ ###
### ----------------------------- Heatmaps --------------------------------- ###
### ------------------------------------------------------------------------ ###

### Rearrange =====================================================
## CCE
cceSOI <- read.csv(paste0(datap,"CCE_regress_SOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
cceSOI$Lag <- as.factor(cceSOI$Lag)
cceSOI$Type <- factor(cceSOI$Type, 
                         levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rcceSOI <- cceSOI[,c("Lag","Type","coef","p")]
rcceSOI$p[rcceSOI$p > 0.1] <- NA
rcceSOI$sym <- rcceSOI$p
rcceSOI$sym[rcceSOI$p <= 0.05] <- "*"
rcceSOI$sym[rcceSOI$p > 0.05] <- "+"
rcceSOI$sym[rcceSOI$p > 0.1] <- NA

cceNino34 <- read.csv(paste0(datap,"CCE_regress_Nino34_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
cceNino34$Lag <- as.factor(cceNino34$Lag)
cceNino34$Type <- factor(cceNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rcceNino34 <- cceNino34[,c("Lag","Type","coef","p")]
rcceNino34$p[rcceNino34$p > 0.1] <- NA
rcceNino34$sym <- rcceNino34$p
rcceNino34$sym[rcceNino34$p <= 0.05] <- "*"
rcceNino34$sym[rcceNino34$p > 0.05] <- "+"
rcceNino34$sym[rcceNino34$p > 0.1] <- NA

cceNino12 <- read.csv(paste0(datap,"CCE_regress_Nino12_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
cceNino12$Lag <- as.factor(cceNino12$Lag)
cceNino12$Type <- factor(cceNino12$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rcceNino12 <- cceNino12[,c("Lag","Type","coef","p")]
rcceNino12$p[rcceNino12$p > 0.1] <- NA
rcceNino12$sym <- rcceNino12$p
rcceNino12$sym[rcceNino12$p <= 0.05] <- "*"
rcceNino12$sym[rcceNino12$p > 0.05] <- "+"
rcceNino12$sym[rcceNino12$p > 0.1] <- NA

cceNOI <- read.csv(paste0(datap,"CCE_regress_NOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
cceNOI$Lag <- as.factor(cceNOI$Lag)
cceNOI$Type <- factor(cceNOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rcceNOI <- cceNOI[,c("Lag","Type","coef","p")]
rcceNOI$p[rcceNOI$p > 0.1] <- NA
rcceNOI$sym <- rcceNOI$p
rcceNOI$sym[rcceNOI$p <= 0.05] <- "*"
rcceNOI$sym[rcceNOI$p > 0.05] <- "+"
rcceNOI$sym[rcceNOI$p > 0.1] <- NA

cceMEI <- read.csv(paste0(datap,"CCE_regress_MEI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
cceMEI$Lag <- as.factor(cceMEI$Lag)
cceMEI$Type <- factor(cceMEI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rcceMEI <- cceMEI[,c("Lag","Type","coef","p")]
rcceMEI$p[rcceMEI$p > 0.1] <- NA
rcceMEI$sym <- rcceMEI$p
rcceMEI$sym[rcceMEI$p <= 0.05] <- "*"
rcceMEI$sym[rcceMEI$p > 0.05] <- "+"
rcceMEI$sym[rcceMEI$p > 0.1] <- NA


## EBS
ebsSOI <- read.csv(paste0(datap,"EBS_regress_SOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
ebsSOI$Lag <- as.factor(ebsSOI$Lag)
ebsSOI$Type <- factor(ebsSOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rebsSOI <- ebsSOI[,c("Lag","Type","coef","p")]
rebsSOI$p[rebsSOI$p > 0.1] <- NA
rebsSOI$sym <- rebsSOI$p
rebsSOI$sym[rebsSOI$p <= 0.05] <- "*"
rebsSOI$sym[rebsSOI$p > 0.05] <- "+"
rebsSOI$sym[rebsSOI$p > 0.1] <- NA

ebsNino34 <- read.csv(paste0(datap,"EBS_regress_Nino34_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
ebsNino34$Lag <- as.factor(ebsNino34$Lag)
ebsNino34$Type <- factor(ebsNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rebsNino34 <- ebsNino34[,c("Lag","Type","coef","p")]
rebsNino34$p[rebsNino34$p > 0.1] <- NA
rebsNino34$sym <- rebsNino34$p
rebsNino34$sym[rebsNino34$p <= 0.05] <- "*"
rebsNino34$sym[rebsNino34$p > 0.05] <- "+"
rebsNino34$sym[rebsNino34$p > 0.1] <- NA

ebsNino12 <- read.csv(paste0(datap,"EBS_regress_Nino12_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
ebsNino12$Lag <- as.factor(ebsNino12$Lag)
ebsNino12$Type <- factor(ebsNino12$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rebsNino12 <- ebsNino12[,c("Lag","Type","coef","p")]
rebsNino12$p[rebsNino12$p > 0.1] <- NA
rebsNino12$sym <- rebsNino12$p
rebsNino12$sym[rebsNino12$p <= 0.05] <- "*"
rebsNino12$sym[rebsNino12$p > 0.05] <- "+"
rebsNino12$sym[rebsNino12$p > 0.1] <- NA

ebsNOI <- read.csv(paste0(datap,"EBS_regress_NOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
ebsNOI$Lag <- as.factor(ebsNOI$Lag)
ebsNOI$Type <- factor(ebsNOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rebsNOI <- ebsNOI[,c("Lag","Type","coef","p")]
rebsNOI$p[rebsNOI$p > 0.1] <- NA
rebsNOI$sym <- rebsNOI$p
rebsNOI$sym[rebsNOI$p <= 0.05] <- "*"
rebsNOI$sym[rebsNOI$p > 0.05] <- "+"
rebsNOI$sym[rebsNOI$p > 0.1] <- NA

ebsMEI <- read.csv(paste0(datap,"EBS_regress_MEI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
ebsMEI$Lag <- as.factor(ebsMEI$Lag)
ebsMEI$Type <- factor(ebsMEI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rebsMEI <- ebsMEI[,c("Lag","Type","coef","p")]
rebsMEI$p[rebsMEI$p > 0.1] <- NA
rebsMEI$sym <- rebsMEI$p
rebsMEI$sym[rebsMEI$p <= 0.05] <- "*"
rebsMEI$sym[rebsMEI$p > 0.05] <- "+"
rebsMEI$sym[rebsMEI$p > 0.1] <- NA


## GAK
gakSOI <- read.csv(paste0(datap,"GAK_regress_SOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gakSOI$Lag <- as.factor(gakSOI$Lag)
gakSOI$Type <- factor(gakSOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rgakSOI <- gakSOI[,c("Lag","Type","coef","p")]
rgakSOI$p[rgakSOI$p > 0.1] <- NA
rgakSOI$sym <- rgakSOI$p
rgakSOI$sym[rgakSOI$p <= 0.05] <- "*"
rgakSOI$sym[rgakSOI$p > 0.05] <- "+"
rgakSOI$sym[rgakSOI$p > 0.1] <- NA

gakNino34 <- read.csv(paste0(datap,"GAK_regress_Nino34_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gakNino34$Lag <- as.factor(gakNino34$Lag)
gakNino34$Type <- factor(gakNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rgakNino34 <- gakNino34[,c("Lag","Type","coef","p")]
rgakNino34$p[rgakNino34$p > 0.1] <- NA
rgakNino34$sym <- rgakNino34$p
rgakNino34$sym[rgakNino34$p <= 0.05] <- "*"
rgakNino34$sym[rgakNino34$p > 0.05] <- "+"
rgakNino34$sym[rgakNino34$p > 0.1] <- NA

gakNino12 <- read.csv(paste0(datap,"GAK_regress_Nino12_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gakNino12$Lag <- as.factor(gakNino12$Lag)
gakNino12$Type <- factor(gakNino12$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rgakNino12 <- gakNino12[,c("Lag","Type","coef","p")]
rgakNino12$p[rgakNino12$p > 0.1] <- NA
rgakNino12$sym <- rgakNino12$p
rgakNino12$sym[rgakNino12$p <= 0.05] <- "*"
rgakNino12$sym[rgakNino12$p > 0.05] <- "+"
rgakNino12$sym[rgakNino12$p > 0.1] <- NA

gakNOI <- read.csv(paste0(datap,"GAK_regress_NOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gakNOI$Lag <- as.factor(gakNOI$Lag)
gakNOI$Type <- factor(gakNOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rgakNOI <- gakNOI[,c("Lag","Type","coef","p")]
rgakNOI$p[rgakNOI$p > 0.1] <- NA
rgakNOI$sym <- rgakNOI$p
rgakNOI$sym[rgakNOI$p <= 0.05] <- "*"
rgakNOI$sym[rgakNOI$p > 0.05] <- "+"
rgakNOI$sym[rgakNOI$p > 0.1] <- NA

gakMEI <- read.csv(paste0(datap,"GAK_regress_MEI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gakMEI$Lag <- as.factor(gakMEI$Lag)
gakMEI$Type <- factor(gakMEI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rgakMEI <- gakMEI[,c("Lag","Type","coef","p")]
rgakMEI$p[rgakMEI$p > 0.1] <- NA
rgakMEI$sym <- rgakMEI$p
rgakMEI$sym[rgakMEI$p <= 0.05] <- "*"
rgakMEI$sym[rgakMEI$p > 0.05] <- "+"
rgakMEI$sym[rgakMEI$p > 0.1] <- NA


## HI
hiSOI <- read.csv(paste0(datap,"HI_regress_SOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
hiSOI$Lag <- as.factor(hiSOI$Lag)
hiSOI$Type <- factor(hiSOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rhiSOI <- hiSOI[,c("Lag","Type","coef","p")]
rhiSOI$p[rhiSOI$p > 0.1] <- NA
rhiSOI$sym <- rhiSOI$p
rhiSOI$sym[rhiSOI$p <= 0.05] <- "*"
rhiSOI$sym[rhiSOI$p > 0.05] <- "+"
rhiSOI$sym[rhiSOI$p > 0.1] <- NA

hiNino34 <- read.csv(paste0(datap,"HI_regress_Nino34_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
hiNino34$Lag <- as.factor(hiNino34$Lag)
hiNino34$Type <- factor(hiNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rhiNino34 <- hiNino34[,c("Lag","Type","coef","p")]
rhiNino34$p[rhiNino34$p > 0.1] <- NA
rhiNino34$sym <- rhiNino34$p
rhiNino34$sym[rhiNino34$p <= 0.05] <- "*"
rhiNino34$sym[rhiNino34$p > 0.05] <- "+"
rhiNino34$sym[rhiNino34$p > 0.1] <- NA

hiNino12 <- read.csv(paste0(datap,"HI_regress_Nino12_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
hiNino12$Lag <- as.factor(hiNino12$Lag)
hiNino12$Type <- factor(hiNino12$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rhiNino12 <- hiNino12[,c("Lag","Type","coef","p")]
rhiNino12$p[rhiNino12$p > 0.1] <- NA
rhiNino12$sym <- rhiNino12$p
rhiNino12$sym[rhiNino12$p <= 0.05] <- "*"
rhiNino12$sym[rhiNino12$p > 0.05] <- "+"
rhiNino12$sym[rhiNino12$p > 0.1] <- NA

hiNOI <- read.csv(paste0(datap,"HI_regress_NOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
hiNOI$Lag <- as.factor(hiNOI$Lag)
hiNOI$Type <- factor(hiNOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rhiNOI <- hiNOI[,c("Lag","Type","coef","p")]
rhiNOI$p[rhiNOI$p > 0.1] <- NA
rhiNOI$sym <- rhiNOI$p
rhiNOI$sym[rhiNOI$p <= 0.05] <- "*"
rhiNOI$sym[rhiNOI$p > 0.05] <- "+"
rhiNOI$sym[rhiNOI$p > 0.1] <- NA

hiMEI <- read.csv(paste0(datap,"HI_regress_MEI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
hiMEI$Lag <- as.factor(hiMEI$Lag)
hiMEI$Type <- factor(hiMEI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rhiMEI <- hiMEI[,c("Lag","Type","coef","p")]
rhiMEI$p[rhiMEI$p > 0.1] <- NA
rhiMEI$sym <- rhiMEI$p
rhiMEI$sym[rhiMEI$p <= 0.05] <- "*"
rhiMEI$sym[rhiMEI$p > 0.05] <- "+"
rhiMEI$sym[rhiMEI$p > 0.1] <- NA


## CHK
chkSOI <- read.csv(paste0(datap,"CHK_regress_SOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
chkSOI$Lag <- as.factor(chkSOI$Lag)
chkSOI$Type <- factor(chkSOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rchkSOI <- chkSOI[,c("Lag","Type","coef","p")]
rchkSOI$p[rchkSOI$p > 0.1] <- NA
rchkSOI$sym <- rchkSOI$p
rchkSOI$sym[rchkSOI$p <= 0.05] <- "*"
rchkSOI$sym[rchkSOI$p > 0.05] <- "+"
rchkSOI$sym[rchkSOI$p > 0.1] <- NA

chkNino34 <- read.csv(paste0(datap,"CHK_regress_Nino34_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
chkNino34$Lag <- as.factor(chkNino34$Lag)
chkNino34$Type <- factor(chkNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rchkNino34 <- chkNino34[,c("Lag","Type","coef","p")]
rchkNino34$p[rchkNino34$p > 0.1] <- NA
rchkNino34$sym <- rchkNino34$p
rchkNino34$sym[rchkNino34$p <= 0.05] <- "*"
rchkNino34$sym[rchkNino34$p > 0.05] <- "+"
rchkNino34$sym[rchkNino34$p > 0.1] <- NA

chkNino12 <- read.csv(paste0(datap,"CHK_regress_Nino12_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
chkNino12$Lag <- as.factor(chkNino12$Lag)
chkNino12$Type <- factor(chkNino12$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rchkNino12 <- chkNino12[,c("Lag","Type","coef","p")]
rchkNino12$p[rchkNino12$p > 0.1] <- NA
rchkNino12$sym <- rchkNino12$p
rchkNino12$sym[rchkNino12$p <= 0.05] <- "*"
rchkNino12$sym[rchkNino12$p > 0.05] <- "+"
rchkNino12$sym[rchkNino12$p > 0.1] <- NA

chkNOI <- read.csv(paste0(datap,"CHK_regress_NOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
chkNOI$Lag <- as.factor(chkNOI$Lag)
chkNOI$Type <- factor(chkNOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rchkNOI <- chkNOI[,c("Lag","Type","coef","p")]
rchkNOI$p[rchkNOI$p > 0.1] <- NA
rchkNOI$sym <- rchkNOI$p
rchkNOI$sym[rchkNOI$p <= 0.05] <- "*"
rchkNOI$sym[rchkNOI$p > 0.05] <- "+"
rchkNOI$sym[rchkNOI$p > 0.1] <- NA

chkMEI <- read.csv(paste0(datap,"CHK_regress_MEI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
chkMEI$Lag <- as.factor(chkMEI$Lag)
chkMEI$Type <- factor(chkMEI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rchkMEI <- chkMEI[,c("Lag","Type","coef","p")]
rchkMEI$p[rchkMEI$p > 0.1] <- NA
rchkMEI$sym <- rchkMEI$p
rchkMEI$sym[rchkMEI$p <= 0.05] <- "*"
rchkMEI$sym[rchkMEI$p > 0.05] <- "+"
rchkMEI$sym[rchkMEI$p > 0.1] <- NA


##GMX
gmxSOI <- read.csv(paste0(datap,"MX_regress_SOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gmxSOI$Lag <- as.factor(gmxSOI$Lag)
gmxSOI$Type <- factor(gmxSOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rgmxSOI <- gmxSOI[,c("Lag","Type","coef","p")]
rgmxSOI$p[rgmxSOI$p > 0.1] <- NA
rgmxSOI$sym <- rgmxSOI$p
rgmxSOI$sym[rgmxSOI$p <= 0.05] <- "*"
rgmxSOI$sym[rgmxSOI$p > 0.05] <- "+"
rgmxSOI$sym[rgmxSOI$p > 0.1] <- NA

gmxNino34 <- read.csv(paste0(datap,"MX_regress_Nino34_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gmxNino34$Lag <- as.factor(gmxNino34$Lag)
gmxNino34$Type <- factor(gmxNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rgmxNino34 <- gmxNino34[,c("Lag","Type","coef","p")]
rgmxNino34$p[rgmxNino34$p > 0.1] <- NA
rgmxNino34$sym <- rgmxNino34$p
rgmxNino34$sym[rgmxNino34$p <= 0.05] <- "*"
rgmxNino34$sym[rgmxNino34$p > 0.05] <- "+"
rgmxNino34$sym[rgmxNino34$p > 0.1] <- NA

gmxNino12 <- read.csv(paste0(datap,"MX_regress_Nino12_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gmxNino12$Lag <- as.factor(gmxNino12$Lag)
gmxNino12$Type <- factor(gmxNino12$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rgmxNino12 <- gmxNino12[,c("Lag","Type","coef","p")]
rgmxNino12$p[rgmxNino12$p > 0.1] <- NA
rgmxNino12$sym <- rgmxNino12$p
rgmxNino12$sym[rgmxNino12$p <= 0.05] <- "*"
rgmxNino12$sym[rgmxNino12$p > 0.05] <- "+"
rgmxNino12$sym[rgmxNino12$p > 0.1] <- NA

gmxNOI <- read.csv(paste0(datap,"MX_regress_NOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gmxNOI$Lag <- as.factor(gmxNOI$Lag)
gmxNOI$Type <- factor(gmxNOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rgmxNOI <- gmxNOI[,c("Lag","Type","coef","p")]
rgmxNOI$p[rgmxNOI$p > 0.1] <- NA
rgmxNOI$sym <- rgmxNOI$p
rgmxNOI$sym[rgmxNOI$p <= 0.05] <- "*"
rgmxNOI$sym[rgmxNOI$p > 0.05] <- "+"
rgmxNOI$sym[rgmxNOI$p > 0.1] <- NA

gmxMEI <- read.csv(paste0(datap,"MX_regress_MEI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
gmxMEI$Lag <- as.factor(gmxMEI$Lag)
gmxMEI$Type <- factor(gmxMEI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rgmxMEI <- gmxMEI[,c("Lag","Type","coef","p")]
rgmxMEI$p[rgmxMEI$p > 0.1] <- NA
rgmxMEI$sym <- rgmxMEI$p
rgmxMEI$sym[rgmxMEI$p <= 0.05] <- "*"
rgmxMEI$sym[rgmxMEI$p > 0.05] <- "+"
rgmxMEI$sym[rgmxMEI$p > 0.1] <- NA


## NE
neSOI <- read.csv(paste0(datap,"NE_regress_SOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
neSOI$Lag <- as.factor(neSOI$Lag)
neSOI$Type <- factor(neSOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rneSOI <- neSOI[,c("Lag","Type","coef","p")]
rneSOI$p[rneSOI$p > 0.1] <- NA
rneSOI$sym <- rneSOI$p
rneSOI$sym[rneSOI$p <= 0.05] <- "*"
rneSOI$sym[rneSOI$p > 0.05] <- "+"
rneSOI$sym[rneSOI$p > 0.1] <- NA

neNino34 <- read.csv(paste0(datap,"NE_regress_Nino34_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
neNino34$Lag <- as.factor(neNino34$Lag)
neNino34$Type <- factor(neNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rneNino34 <- neNino34[,c("Lag","Type","coef","p")]
rneNino34$p[rneNino34$p > 0.1] <- NA
rneNino34$sym <- rneNino34$p
rneNino34$sym[rneNino34$p <= 0.05] <- "*"
rneNino34$sym[rneNino34$p > 0.05] <- "+"
rneNino34$sym[rneNino34$p > 0.1] <- NA

neNino12 <- read.csv(paste0(datap,"NE_regress_Nino12_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
neNino12$Lag <- as.factor(neNino12$Lag)
neNino12$Type <- factor(neNino12$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rneNino12 <- neNino12[,c("Lag","Type","coef","p")]
rneNino12$p[rneNino12$p > 0.1] <- NA
rneNino12$sym <- rneNino12$p
rneNino12$sym[rneNino12$p <= 0.05] <- "*"
rneNino12$sym[rneNino12$p > 0.05] <- "+"
rneNino12$sym[rneNino12$p > 0.1] <- NA

neNOI <- read.csv(paste0(datap,"NE_regress_NOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
neNOI$Lag <- as.factor(neNOI$Lag)
neNOI$Type <- factor(neNOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rneNOI <- neNOI[,c("Lag","Type","coef","p")]
rneNOI$p[rneNOI$p > 0.1] <- NA
rneNOI$sym <- rneNOI$p
rneNOI$sym[rneNOI$p <= 0.05] <- "*"
rneNOI$sym[rneNOI$p > 0.05] <- "+"
rneNOI$sym[rneNOI$p > 0.1] <- NA

neMEI <- read.csv(paste0(datap,"NE_regress_MEI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
neMEI$Lag <- as.factor(neMEI$Lag)
neMEI$Type <- factor(neMEI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rneMEI <- neMEI[,c("Lag","Type","coef","p")]
rneMEI$p[rneMEI$p > 0.1] <- NA
rneMEI$sym <- rneMEI$p
rneMEI$sym[rneMEI$p <= 0.05] <- "*"
rneMEI$sym[rneMEI$p > 0.05] <- "+"
rneMEI$sym[rneMEI$p > 0.1] <- NA


##SE
seSOI <- read.csv(paste0(datap,"SE_regress_SOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
seSOI$Lag <- as.factor(seSOI$Lag)
seSOI$Type <- factor(seSOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rseSOI <- seSOI[,c("Lag","Type","coef","p")]
rseSOI$p[rseSOI$p > 0.1] <- NA
rseSOI$sym <- rseSOI$p
rseSOI$sym[rseSOI$p <= 0.05] <- "*"
rseSOI$sym[rseSOI$p > 0.05] <- "+"
rseSOI$sym[rseSOI$p > 0.1] <- NA

seNino34 <- read.csv(paste0(datap,"SE_regress_Nino34_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
seNino34$Lag <- as.factor(seNino34$Lag)
seNino34$Type <- factor(seNino34$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rseNino34 <- seNino34[,c("Lag","Type","coef","p")]
rseNino34$p[rseNino34$p > 0.1] <- NA
rseNino34$sym <- rseNino34$p
rseNino34$sym[rseNino34$p <= 0.05] <- "*"
rseNino34$sym[rseNino34$p > 0.05] <- "+"
rseNino34$sym[rseNino34$p > 0.1] <- NA

seNino12 <- read.csv(paste0(datap,"SE_regress_Nino12_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
seNino12$Lag <- as.factor(seNino12$Lag)
seNino12$Type <- factor(seNino12$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rseNino12 <- seNino12[,c("Lag","Type","coef","p")]
rseNino12$p[rseNino12$p > 0.1] <- NA
rseNino12$sym <- rseNino12$p
rseNino12$sym[rseNino12$p <= 0.05] <- "*"
rseNino12$sym[rseNino12$p > 0.05] <- "+"
rseNino12$sym[rseNino12$p > 0.1] <- NA

seNOI <- read.csv(paste0(datap,"SE_regress_NOI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
seNOI$Lag <- as.factor(seNOI$Lag)
seNOI$Type <- factor(seNOI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rseNOI <- seNOI[,c("Lag","Type","coef","p")]
rseNOI$p[rseNOI$p > 0.1] <- NA
rseNOI$sym <- rseNOI$p
rseNOI$sym[rseNOI$p <= 0.05] <- "*"
rseNOI$sym[rseNOI$p > 0.05] <- "+"
rseNOI$sym[rseNOI$p > 0.1] <- NA

seMEI <- read.csv(paste0(datap,"SE_regress_MEI_div2SD_melt.csv"),sep=",",header = T,stringsAsFactors = F)
seMEI$Lag <- as.factor(seMEI$Lag)
seMEI$Type <- factor(seMEI$Type, 
                      levels = c("Tb","Det","B","D","A","P","F","L","M","S","LZloss","LZbiom","Tp"))
rseMEI <- seMEI[,c("Lag","Type","coef","p")]
rseMEI$p[rseMEI$p > 0.1] <- NA
rseMEI$sym <- rseMEI$p
rseMEI$sym[rseMEI$p <= 0.05] <- "*"
rseMEI$sym[rseMEI$p > 0.05] <- "+"
rseMEI$sym[rseMEI$p > 0.1] <- NA


###====== Plots ===================================================================
## CCE
c5 <- ggplot(data = rcceSOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.80,0.80),  
                       name="SOI\nCoef") +
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

c1 <- ggplot(data = rcceNino12, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.80,0.80),  
                       name="Nino12\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)

c2 <- ggplot(data = rcceNOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.80,0.80),  
                       name="NOI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("CCE") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3) 
  

c3 <- ggplot(data = rcceMEI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.80,0.80),  
                       name="MEI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_CCE_coef_climate2_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c1,c2,c3,c4,c5,
           nrow = 2, ncol = 3)
dev.off()


## EBS
e5 <- ggplot(data = rebsSOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.05,1.05),  
                       name="SOI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

e4 <- ggplot(data = rebsNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.05,1.05),  
                       name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

e1 <- ggplot(data = rebsNino12, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.05,1.05),  
                       name="Nino12\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

e2 <- ggplot(data = rebsNOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.05,1.05),  
                       name="NOI\nCoef") +
  theme_minimal() + ggtitle("EBS") + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

e3 <- ggplot(data = rebsMEI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.05,1.05),  
                       name="MEI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_EBS_coef_climate2_inputs_fish.png'), 
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
g5 <- ggplot(data = rgakSOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="SOI\nCoef") +
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

g1 <- ggplot(data = rgakNino12, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino12\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

g2 <- ggplot(data = rgakNOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("GAK") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

g3 <- ggplot(data = rgakMEI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="MEI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_GAK_coef_climate2_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( g1,g2,g3,g4,g5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## HI
h5 <- ggplot(data = rhiSOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.50,0.50),  
                       name="SOI\nCoef") +
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

h1 <- ggplot(data = rhiNino12, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.50,0.50),  
                       name="Nino12\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

h2 <- ggplot(data = rhiNOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.50,0.50),  
                       name="NOI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("HI") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

h3 <- ggplot(data = rhiMEI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.50,0.50),  
                       name="MEI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_HI_coef_climate2_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( h1,h2,h3,h4,h5,
           nrow = 2, ncol = 3)#,
           #align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## CHK
k5 <- ggplot(data = rchkSOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.2,1.2),  
                       name="SOI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

k4 <- ggplot(data = rchkNino34, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.2,1.2),  
                       name="Nino34\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

k1 <- ggplot(data = rchkNino12, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.91,0.91),  
                       name="Nino12\nCoef") +
  theme_minimal() + ggtitle("") + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

k2 <- ggplot(data = rchkNOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.2,1.2),  
                       name="NOI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("CHK") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

k3 <- ggplot(data = rchkMEI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.2,1.2),  
                       name="MEI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_CHK_coef_climate2_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( k1,k2,k3,k4,k5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## GMX
m5 <- ggplot(data = rgmxSOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.45,0.45),  
                       name="SOI\nCoef") +
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

m1 <- ggplot(data = rgmxNino12, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.45,0.45),  
                       name="Nino12\nCoef") +
  theme_minimal() + ggtitle("") + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

m2 <- ggplot(data = rgmxNOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.99,0.99),  
                       name="NOI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("GMX") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

m3 <- ggplot(data = rgmxMEI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.45,0.45),  
                       name="MEI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_GMX_coef_climate2_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( m1,m2,m3,m4,m5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## SE
s5 <- ggplot(data = rseSOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.55,0.55),  
                       name="SOI\nCoef") +
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

s1 <- ggplot(data = rseNino12, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.55,0.55),  
                       name="Nino12\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

s2 <- ggplot(data = rseNOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.55,0.55),  
                       name="NOI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("SEUS") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

s3 <- ggplot(data = rseMEI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.55,0.55),  
                       name="MEI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_SE_coef_climate2_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( s1,s2,s3,s4,s5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## NE
n5 <- ggplot(data = rneSOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.25,1.25),  
                       name="SOI\nCoef") +
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

n1 <- ggplot(data = rneNino12, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.25,1.25),  
                       name="Nino12\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

n2 <- ggplot(data = rneNOI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.25,1.25),  
                       name="NOI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("NEUS") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

n3 <- ggplot(data = rneMEI, aes(y=Type, x=Lag, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-1.25,1.25),  
                       name="MEI\nCoef") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 
#geom_text(aes(Lag, Type, label = signif(p,digits = 1)), color = "black", size = 3)  

png(paste0(figp,'Heatmaps_NE_coef_climate2_inputs_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( n1,n2,n3,n4,n5,
           nrow = 2, ncol = 3)#, labels = "auto", label_size = 12, hjust = -4)
dev.off()






