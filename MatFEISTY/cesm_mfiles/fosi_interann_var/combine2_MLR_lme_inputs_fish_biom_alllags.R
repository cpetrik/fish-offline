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


fcoef1 <- read.csv(paste0(datar,"LME2125_F_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef1 <- read.csv(paste0(datar,"LME2125_P_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef1 <- read.csv(paste0(datar,"LME2125_D_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
acoef1 <- read.csv(paste0(datar,"LME2125_A_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bcoef1 <- read.csv(paste0(datar,"LME2125_B_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)

fcoef2 <- read.csv(paste0(datar,"LMEmost_F_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef2 <- read.csv(paste0(datar,"LMEmost_P_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef2 <- read.csv(paste0(datar,"LMEmost_D_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
acoef2 <- read.csv(paste0(datar,"LMEmost_A_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bcoef2 <- read.csv(paste0(datar,"LMEmost_B_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)

fpval1 <- read.csv(paste0(datar,"LME2125_F_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
ppval1 <- read.csv(paste0(datar,"LME2125_P_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dpval1 <- read.csv(paste0(datar,"LME2125_D_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
apval1 <- read.csv(paste0(datar,"LME2125_A_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bpval1 <- read.csv(paste0(datar,"LME2125_B_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)

fpval2 <- read.csv(paste0(datar,"LMEmost_F_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
ppval2 <- read.csv(paste0(datar,"LMEmost_P_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dpval2 <- read.csv(paste0(datar,"LMEmost_D_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
apval2 <- read.csv(paste0(datar,"LMEmost_A_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bpval2 <- read.csv(paste0(datar,"LMEmost_B_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)


# 1 = 1:24, 2 = all except 21-35,62
fcoef1 <- fcoef1[1:24,]
pcoef1 <- pcoef1[1:24,]
dcoef1 <- dcoef1[1:24,]
acoef1 <- acoef1[1:24,]
bcoef1 <- bcoef1[1:24,]

fcoef2 <- fcoef2[22:49,]
pcoef2 <- pcoef2[22:49,]
dcoef2 <- dcoef2[22:49,]
acoef2 <- acoef2[22:49,]
bcoef2 <- bcoef2[22:49,]

fpval1 <- fpval1[1:24,]
ppval1 <- ppval1[1:24,]
dpval1 <- dpval1[1:24,]
apval1 <- apval1[1:24,]
bpval1 <- bpval1[1:24,]

fpval2 <- fpval2[22:49,]
ppval2 <- ppval2[22:49,]
dpval2 <- dpval2[22:49,]
apval2 <- apval2[22:49,]
bpval2 <- bpval2[22:49,]



### Combine
fcoef <- rbind(fcoef1,fcoef2)
pcoef <- rbind(pcoef1,pcoef2)
dcoef <- rbind(dcoef1,dcoef2)
acoef <- rbind(acoef1,acoef2)
bcoef <- rbind(bcoef1,bcoef2)

fpval <- rbind(fpval1,fpval2)
ppval <- rbind(ppval1,ppval2)
dpval <- rbind(dpval1,dpval2)
apval <- rbind(apval1,apval2)
bpval <- rbind(bpval1,bpval2)


### MLR of drivers ---------------------------------------------------------

# Missing 25-34 & 62

write.table(fcoef,paste0(datar,"LMEmost2_F_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datar,"LMEmost2_P_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datar,"LMEmost2_D_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datar,"LMEmost2_A_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(bcoef,paste0(datar,"LMEmost2_B_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)

write.table(fpval,paste0(datar,"LMEmost2_F_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(ppval,paste0(datar,"LMEmost2_P_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(dpval,paste0(datar,"LMEmost2_D_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(apval,paste0(datar,"LMEmost2_A_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(bpval,paste0(datar,"LMEmost2_B_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)

