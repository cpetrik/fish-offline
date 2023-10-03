# Mult linear regression of forcing on fish biomass
# Run once (10/1/2023) & save output
# Create separate code for plotting results
# Include different lags in MLR, let AICc select
# Max lag is 2yrs
# Remove zmeso biom
# Save R2
# Remove interactions (should I incr time lags since can have more covariates?)

rm(list=ls())

library("MuMIn")


source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/corrs/"

### ------------------------------------------------------------
# load data
datap <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/"
TP <- read.csv(paste0(datap,"LMEs_anoms_TP_trans.csv"),sep=",",header = T,stringsAsFactors = F)
TB <- read.csv(paste0(datap,"LMEs_anoms_TB_trans.csv"),sep=",",header = T,stringsAsFactors = F)
Det <- read.csv(paste0(datap,"LMEs_anoms_Det_trans.csv"),sep=",",header = T,stringsAsFactors = F)
ZmLoss <- read.csv(paste0(datap,"LMEs_anoms_ZmLoss_trans.csv"),sep=",",header = T,stringsAsFactors = F)

#dataf <- "/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/"
dataf <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/Fish-MIP/Phase3/fishing/"
F2 <- read.csv(paste0(dataf,"LMEs_FishMIP_cpue_anoms_div2SD_F.csv"),sep=",",header = T,stringsAsFactors = F)
P2 <- read.csv(paste0(dataf,"LMEs_FishMIP_cpue_anoms_div2SD_P.csv"),sep=",",header = T,stringsAsFactors = F)
D2 <- read.csv(paste0(dataf,"LMEs_FishMIP_cpue_anoms_div2SD_D.csv"),sep=",",header = T,stringsAsFactors = F)
A2 <- read.csv(paste0(dataf,"LMEs_FishMIP_cpue_anoms_div2SD_A.csv"),sep=",",header = T,stringsAsFactors = F)

datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"

cname <- names(Det)
nlme <- length(cname)
lname <- names(F2)[2:64]

### Divide by 2SD 
TP2 <- TP
TB2 <- TB
Det2 <- Det
ZmLoss2 <- ZmLoss

for (i in 1:nlme) {
  TP2[,i] <- TP[,i]/(2*sd(TP[,i]))
  TB2[,i] <- TB[,i]/(2*sd(TB[,i]))
  Det2[,i] <- Det[,i]/(2*sd(Det[,i]))
  ZmLoss2[,i] <- ZmLoss[,i]/(2*sd(ZmLoss[,i]))
}

names(TP2) <- lname
names(TB2) <- lname
names(Det2) <- lname
names(ZmLoss2) <- lname

TP2$Year <- c(1948:2015)
TB2$Year <- c(1948:2015)
Det2$Year <- c(1948:2015)
ZmLoss2$Year <- c(1948:2015)

TP2 <- subset.data.frame(TP2, Year >= 1961 & Year <= 2010)
TB2 <- subset.data.frame(TB2, Year >= 1961 & Year <= 2010)
Det2 <- subset.data.frame(Det2, Year >= 1961 & Year <= 2010)
ZmLoss2 <- subset.data.frame(ZmLoss2, Year >= 1961 & Year <= 2010)

#LME 64 of P all NaNs
P2$LME64 <- 0


### MLR of drivers ---------------------------------------------------------
dname <- c('TP','TB','ZmLoss','Det')

#Loop over LMEs  1:nlme
for (i in 1:nlme) {
  cn <- lname[i]
  
  ### Drivers
  maxl <- 3
  yst <- 1+maxl
  yen <- 50+maxl
  
  drive <- data.frame(matrix(ncol = (4*3), nrow = yen))
  names(drive) <- c('TP0','TP1','TP2',
                    'TB0','TB1','TB2',
                    'Det0','Det1','Det2',
                    'ZmLoss0','ZmLoss1','ZmLoss2')
  
  drive$TP0[(yst-3):(yen-3)] <- TP2[,cn]
  drive$TP1[(yst-2):(yen-2)] <- TP2[,cn]
  drive$TP2[(yst-1):(yen-1)] <- TP2[,cn]
  
  drive$TB0[(yst-3):(yen-3)] <- TB2[,cn]
  drive$TB1[(yst-2):(yen-2)] <- TB2[,cn]
  drive$TB2[(yst-1):(yen-1)] <- TB2[,cn]
  
  drive$Det0[(yst-3):(yen-3)] <- Det2[,cn]
  drive$Det1[(yst-2):(yen-2)] <- Det2[,cn]
  drive$Det2[(yst-1):(yen-1)] <- Det2[,cn]
  
  drive$ZmLoss0[(yst-3):(yen-3)] <- ZmLoss2[,cn]
  drive$ZmLoss1[(yst-2):(yen-2)] <- ZmLoss2[,cn]
  drive$ZmLoss2[(yst-1):(yen-1)] <- ZmLoss2[,cn]
  
  
  
  ### F ----------------------------------------------------------------
  ffish <- drive
  ffish$Fish <- NA
  ffish$Fish[(yst-3):(yen-3)] <- F2[,cn]
  ffish <- na.omit(ffish)
  
  ## Test mod for creating array with results
  # Create arrays for coefficients & p-vals
  options(na.action = "na.fail")
  if (i==1) {
    tmod <- lm(Fish ~ TP0+TP1+TP2 + TB0+TB1+TB2 + Det0+Det1+Det2 + 
                 ZmLoss0+ZmLoss1+ZmLoss2 , data=ffish)
    tcombo <- dredge(tmod, extra = c("R^2", F = function(x)
      summary(x)$fstatistic[[1]]))
    fcoef <- data.frame(matrix(ncol = length(tcombo), nrow = nlme))
    names(fcoef) <- names(tcombo)
    fpval <- fcoef[,1:15]
    fcoef$LME <- cname
    fpval$LME <- cname
    pcoef <- fcoef
    ppval <- fpval
    dcoef <- fcoef
    dpval <- fpval
    acoef <- fcoef
    apval <- fpval
  }
  
  # F actual
  fmod <- lm(Fish ~ TP0+TP1+TP2 + TB0+TB1+TB2 + Det0+Det1+Det2 + 
               ZmLoss0+ZmLoss1+ZmLoss2 , data=ffish)
  fcombo <- dredge(fmod, extra = c("R^2", F = function(x)
    summary(x)$fstatistic[[1]]))
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (fcombo)[1]
  fcoef[i, match(names(xx), colnames(fcoef))] = xx
  y <- summary(get.models(fcombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  fpval[i, match(names(yy), colnames(fpval))] = yy
  
  ### P ----------------------------------------------------------------
  pfish <- drive
  pfish$Fish <- NA
  pfish$Fish[(yst-3):(yen-3)] <- P2[,cn]
  pfish <- na.omit(pfish)
  
  pmod <- lm(Fish ~ TP0+TP1+TP2 + TB0+TB1+TB2 + Det0+Det1+Det2 + 
               ZmLoss0+ZmLoss1+ZmLoss2 , data=pfish)
  pcombo <- dredge(pmod, extra = c("R^2", F = function(x)
    summary(x)$fstatistic[[1]]))
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (pcombo)[1]
  pcoef[i, match(names(xx), colnames(pcoef))] = xx
  y <- summary(get.models(pcombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  ppval[i, match(names(yy), colnames(ppval))] = yy
  
  
  ### D ----------------------------------------------------------------
  dfish <- drive
  dfish$Fish <- NA
  dfish$Fish[(yst-3):(yen-3)] <- D2[,cn]
  dfish <- na.omit(dfish)
  
  dmod <- lm(Fish ~ TP0+TP1+TP2 + TB0+TB1+TB2 + Det0+Det1+Det2 + 
               ZmLoss0+ZmLoss1+ZmLoss2 , data=dfish)
  dcombo <- dredge(dmod, extra = c("R^2", F = function(x)
    summary(x)$fstatistic[[1]]))
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (dcombo)[1]
  dcoef[i, match(names(xx), colnames(dcoef))] = xx
  y <- summary(get.models(dcombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  dpval[i, match(names(yy), colnames(dpval))] = yy
  
  ### A ----------------------------------------------------------------
  afish <- drive
  afish$Fish <- NA
  afish$Fish[(yst-3):(yen-3)] <- A2[,cn]
  afish <- na.omit(afish)
  
  amod <- lm(Fish ~ TP0+TP1+TP2 + TB0+TB1+TB2 + Det0+Det1+Det2 + 
               ZmLoss0+ZmLoss1+ZmLoss2 , data=afish)
  acombo <- dredge(amod, extra = c("R^2", F = function(x)
    summary(x)$fstatistic[[1]]))
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (acombo)[1]
  acoef[i, match(names(xx), colnames(acoef))] = xx
  y <- summary(get.models(acombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  apval[i, match(names(yy), colnames(apval))] = yy
  
  
} # LMEs

write.table(fcoef,paste0(dataf,"LME_FishMIP_cpue_F_mlr_coeffs_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(dataf,"LME_FishMIP_cpue_P_mlr_coeffs_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(dataf,"LME_FishMIP_cpue_D_mlr_coeffs_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",row.names=F)
write.table(acoef,paste0(dataf,"LME_FishMIP_cpue_A_mlr_coeffs_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",row.names=F)

write.table(fpval,paste0(dataf,"LME_FishMIP_cpue_F_mlr_pvals_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",row.names=F)
write.table(ppval,paste0(dataf,"LME_FishMIP_cpue_P_mlr_pvals_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",row.names=F)
write.table(dpval,paste0(dataf,"LME_FishMIP_cpue_D_mlr_pvals_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",row.names=F)
write.table(apval,paste0(dataf,"LME_FishMIP_cpue_A_mlr_pvals_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",row.names=F)

write.table(fcoef,paste0(datap,"LME_FishMIP_cpue_F_mlr_coeffs_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datap,"LME_FishMIP_cpue_P_mlr_coeffs_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datap,"LME_FishMIP_cpue_D_mlr_coeffs_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datap,"LME_FishMIP_cpue_A_mlr_coeffs_ALLdiv2SD_alllags_v2_noint.csv"),sep=",",row.names=F)

