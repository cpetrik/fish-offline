# Mult linear regression of forcing on fish biomass
# Max lag is 3yrs
# Run once (7/22/2023) & save output
# Create separate code for plotting results
# Include different lags in MLR, let AICc select

rm(list=ls())

library("MuMIn")


source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/corrs/"

### ------------------------------------------------------------
# load data
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/FOSI/"
#datap <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/"

datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"


fcoef1 <- read.csv(paste0(datar,"LME21_F_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef1 <- read.csv(paste0(datar,"LME21_P_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef1 <- read.csv(paste0(datar,"LME21_D_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
acoef1 <- read.csv(paste0(datar,"LME21_A_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bcoef1 <- read.csv(paste0(datar,"LME21_B_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)

fcoef2 <- read.csv(paste0(datar,"LME50_F_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef2 <- read.csv(paste0(datar,"LME50_P_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef2 <- read.csv(paste0(datar,"LME50_D_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
acoef2 <- read.csv(paste0(datar,"LME50_A_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bcoef2 <- read.csv(paste0(datar,"LME50_B_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)

fcoef3 <- read.csv(paste0(datar,"LME350_F_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef3 <- read.csv(paste0(datar,"LME350_P_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef3 <- read.csv(paste0(datar,"LME350_D_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
acoef3 <- read.csv(paste0(datar,"LME350_A_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bcoef3 <- read.csv(paste0(datar,"LME350_B_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)

fpval1 <- read.csv(paste0(datar,"LME21_F_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
ppval1 <- read.csv(paste0(datar,"LME21_P_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dpval1 <- read.csv(paste0(datar,"LME21_D_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
apval1 <- read.csv(paste0(datar,"LME21_A_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bpval1 <- read.csv(paste0(datar,"LME21_B_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)

fpval2 <- read.csv(paste0(datar,"LME50_F_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
ppval2 <- read.csv(paste0(datar,"LME50_P_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dpval2 <- read.csv(paste0(datar,"LME50_D_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
apval2 <- read.csv(paste0(datar,"LME50_A_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bpval2 <- read.csv(paste0(datar,"LME50_B_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)

fpval3 <- read.csv(paste0(datar,"LME350_F_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
ppval3 <- read.csv(paste0(datar,"LME350_P_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dpval3 <- read.csv(paste0(datar,"LME350_D_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
apval3 <- read.csv(paste0(datar,"LME350_A_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bpval3 <- read.csv(paste0(datar,"LME350_B_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)


# 1 = 1:21, 2 = 50-61, 3 = 35-50
fcoef1 <- fcoef1[1:21,]
pcoef1 <- pcoef1[1:21,]
dcoef1 <- dcoef1[1:21,]
acoef1 <- acoef1[1:21,]
bcoef1 <- bcoef1[1:21,]

fcoef2 <- fcoef2[50:62,]
pcoef2 <- pcoef2[50:62,]
dcoef2 <- dcoef2[50:62,]
acoef2 <- acoef2[50:62,]
bcoef2 <- bcoef2[50:62,]

fcoef3 <- fcoef3[35:49,]
pcoef3 <- pcoef3[35:49,]
dcoef3 <- dcoef3[35:49,]
acoef3 <- acoef3[35:49,]
bcoef3 <- bcoef3[35:49,]

fpval1 <- fpval1[1:21,]
ppval1 <- ppval1[1:21,]
dpval1 <- dpval1[1:21,]
apval1 <- apval1[1:21,]
bpval1 <- bpval1[1:21,]

fpval2 <- fpval2[50:62,]
ppval2 <- ppval2[50:62,]
dpval2 <- dpval2[50:62,]
apval2 <- apval2[50:62,]
bpval2 <- bpval2[50:62,]

fpval3 <- fpval3[35:49,]
ppval3 <- ppval3[35:49,]
dpval3 <- dpval3[35:49,]
apval3 <- apval3[35:49,]
bpval3 <- bpval3[35:49,]


### Combine
fcoef <- rbind(fcoef1,fcoef3,fcoef2)
pcoef <- rbind(pcoef1,pcoef3,pcoef2)
dcoef <- rbind(dcoef1,dcoef3,dcoef2)
acoef <- rbind(acoef1,acoef3,acoef2)
bcoef <- rbind(bcoef1,bcoef3,bcoef2)

fpval <- rbind(fpval1,fpval3,fpval2)
ppval <- rbind(ppval1,ppval3,ppval2)
dpval <- rbind(dpval1,dpval3,dpval2)
apval <- rbind(apval1,apval3,apval2)
bpval <- rbind(bpval1,bpval3,bpval2)


### MLR of drivers ---------------------------------------------------------

# Missing 22-34 & 62

write.table(fcoef,paste0(datar,"LMEmost_F_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datar,"LMEmost_P_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datar,"LMEmost_D_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datar,"LMEmost_A_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(bcoef,paste0(datar,"LMEmost_B_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)

write.table(fpval,paste0(datar,"LMEmost_F_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(ppval,paste0(datar,"LMEmost_P_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(dpval,paste0(datar,"LMEmost_D_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(apval,paste0(datar,"LMEmost_A_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(bpval,paste0(datar,"LMEmost_B_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)

