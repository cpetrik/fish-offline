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

#datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"
datar <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/"

fcoef1 <- read.csv(paste0(datar,"LME25.35.63.64_F_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef1 <- read.csv(paste0(datar,"LME25.35.63.64_P_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef1 <- read.csv(paste0(datar,"LME25.35.63.64_D_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
acoef1 <- read.csv(paste0(datar,"LME25.35.63.64_A_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bcoef1 <- read.csv(paste0(datar,"LME25.35.63.64_B_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)

fcoef2 <- read.csv(paste0(datar,"LMEmost2_F_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef2 <- read.csv(paste0(datar,"LMEmost2_P_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef2 <- read.csv(paste0(datar,"LMEmost2_D_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
acoef2 <- read.csv(paste0(datar,"LMEmost2_A_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bcoef2 <- read.csv(paste0(datar,"LMEmost2_B_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)

fpval1 <- read.csv(paste0(datar,"LME25.35.63.64_F_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
ppval1 <- read.csv(paste0(datar,"LME25.35.63.64_P_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dpval1 <- read.csv(paste0(datar,"LME25.35.63.64_D_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
apval1 <- read.csv(paste0(datar,"LME25.35.63.64_A_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bpval1 <- read.csv(paste0(datar,"LME25.35.63.64_B_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)

fpval2 <- read.csv(paste0(datar,"LMEmost2_F_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
ppval2 <- read.csv(paste0(datar,"LMEmost2_P_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
dpval2 <- read.csv(paste0(datar,"LMEmost2_D_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
apval2 <- read.csv(paste0(datar,"LMEmost2_A_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)
bpval2 <- read.csv(paste0(datar,"LMEmost2_B_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",header = T,stringsAsFactors = F)


# 1 = 24:34; 62:63, 2 = 1:23; 25:51
fcoef12 <- fcoef1[24:34,]
pcoef12 <- pcoef1[24:34,]
dcoef12 <- dcoef1[24:34,]
acoef12 <- acoef1[24:34,]
bcoef12 <- bcoef1[24:34,]

fcoef11 <- fcoef2[1:23,]
pcoef11 <- pcoef2[1:23,]
dcoef11 <- dcoef2[1:23,]
acoef11 <- acoef2[1:23,]
bcoef11 <- bcoef2[1:23,]

fcoef13 <- fcoef2[25:51,]
pcoef13 <- pcoef2[25:51,]
dcoef13 <- dcoef2[25:51,]
acoef13 <- acoef2[25:51,]
bcoef13 <- bcoef2[25:51,]

fcoef14 <- fcoef1[62:63,]
pcoef14 <- pcoef1[62:63,]
dcoef14 <- dcoef1[62:63,]
acoef14 <- acoef1[62:63,]
bcoef14 <- bcoef1[62:63,]


fpval12 <- fpval1[24:34,]
ppval12 <- ppval1[24:34,]
dpval12 <- dpval1[24:34,]
apval12 <- apval1[24:34,]
bpval12 <- bpval1[24:34,]

fpval11 <- fpval2[1:23,]
ppval11 <- ppval2[1:23,]
dpval11 <- dpval2[1:23,]
apval11 <- apval2[1:23,]
bpval11 <- bpval2[1:23,]

fpval13 <- fpval2[25:51,]
ppval13 <- ppval2[25:51,]
dpval13 <- dpval2[25:51,]
apval13 <- apval2[25:51,]
bpval13 <- bpval2[25:51,]

fpval14 <- fpval1[62:63,]
ppval14 <- ppval1[62:63,]
dpval14 <- dpval1[62:63,]
apval14 <- apval1[62:63,]
bpval14 <- bpval1[62:63,]



### Combine
fcoef <- rbind(fcoef11,fcoef12,fcoef13,fcoef14)
pcoef <- rbind(pcoef11,pcoef12,pcoef13,pcoef14)
dcoef <- rbind(dcoef11,dcoef12,dcoef13,dcoef14)
acoef <- rbind(acoef11,acoef12,acoef13,acoef14)
bcoef <- rbind(bcoef11,bcoef12,bcoef13,bcoef14)

fpval <- rbind(fpval11,fpval12,fpval13,fpval14)
ppval <- rbind(ppval11,ppval12,ppval13,ppval14)
dpval <- rbind(dpval11,dpval12,dpval13,dpval14)
apval <- rbind(apval11,apval12,apval13,apval14)
bpval <- rbind(bpval11,bpval12,bpval13,bpval14)


### MLR of drivers ---------------------------------------------------------

# Missing 25-34 & 62

write.table(fcoef,paste0(datar,"LME_F_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datar,"LME_P_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datar,"LME_D_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datar,"LME_A_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(bcoef,paste0(datar,"LME_B_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)

write.table(fpval,paste0(datar,"LME_F_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(ppval,paste0(datar,"LME_P_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(dpval,paste0(datar,"LME_D_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(apval,paste0(datar,"LME_A_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(bpval,paste0(datar,"LME_B_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)

