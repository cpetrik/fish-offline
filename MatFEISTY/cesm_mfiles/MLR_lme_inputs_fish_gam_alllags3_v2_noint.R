# Mult linear regression of forcing on fish "recruitment" (gamma)
# Run once (2/7/2024) & save output
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
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/FOSI/"
TP2 <- read.csv(paste0(datap,"LMEs_anoms_div2SD_TP.csv"),sep=",",header = T,stringsAsFactors = F)
TB2 <- read.csv(paste0(datap,"LMEs_anoms_div2SD_TB.csv"),sep=",",header = T,stringsAsFactors = F)
Det2 <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Det.csv"),sep=",",header = T,stringsAsFactors = F)
ZmLoss2 <- read.csv(paste0(datap,"LMEs_anoms_div2SD_ZmLoss.csv"),sep=",",header = T,stringsAsFactors = F)

FF <- read.csv(paste0(datap,"LMEs_gam_anoms_div2SD_F.csv"),sep=",",header = T,stringsAsFactors = F)
P <- read.csv(paste0(datap,"LMEs_gam_anoms_div2SD_P.csv"),sep=",",header = T,stringsAsFactors = F)
D <- read.csv(paste0(datap,"LMEs_gam_anoms_div2SD_D.csv"),sep=",",header = T,stringsAsFactors = F)
A <- read.csv(paste0(datap,"LMEs_gam_anoms_div2SD_A.csv"),sep=",",header = T,stringsAsFactors = F)

datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"
#datar <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/"


cname <- names(A)
nlme <- length(cname)

TP2$Year <- c(1:68)
TB2$Year <- c(1:68)
Det2$Year <- c(1:68)
ZmLoss2$Year <- c(1:68)
A$Year <- c(1:68)
D$Year <- c(1:68)
FF$Year <- c(1:68)
P$Year <- c(1:68)

### MLR of drivers ---------------------------------------------------------
dname <- c('TP','TB','ZmLoss','Det')


#Loop over LMEs 
for (i in 1:nlme) {
  cn <- cname[i]
  
  ### Drivers
  maxl <- 3
  yst <- 1+maxl
  yen <- 68+maxl
  
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
  ffish$Fish[(yst-3):(yen-3)] <- FF[,cn]
  ffish <- na.omit(ffish)
  
  options(na.action = "na.fail")
  fmod <- lm(Fish ~ TP0+TP1+TP2 + TB0+TB1+TB2 + Det0+Det1+Det2 + 
               ZmLoss0+ZmLoss1+ZmLoss2, data=ffish)
  fcombo <- dredge(fmod, extra = c("R^2", F = function(x)
    summary(x)$fstatistic[[1]]))
  
  ## Create arrays for coefficients & p-vals
  if (i==1) {
    fcoef <- data.frame(matrix(ncol = length(fcombo), nrow = nlme))
    names(fcoef) <- names(fcombo)
    fpval <- fcoef[,1:15]
    fcoef$LME <- cname
    fpval$LME <- cname
    pcoef <- fcoef
    ppval <- fpval
    dcoef <- fcoef
    dpval <- fpval
    acoef <- fcoef
    apval <- fpval
    bcoef <- fcoef
    bpval <- fpval
  }
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (fcombo)[1]
  fcoef[i, match(names(xx), colnames(fcoef))] = xx
  y <- summary(get.models(fcombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  fpval[i, match(names(yy), colnames(fpval))] = yy
  
  ### P ----------------------------------------------------------------
  pfish <- drive
  pfish$Fish <- NA
  pfish$Fish[(yst-3):(yen-3)] <- P[,cn]
  pfish <- na.omit(pfish)
  
  pmod <- lm(Fish ~ TP0+TP1+TP2 + TB0+TB1+TB2 + Det0+Det1+Det2 + 
               ZmLoss0+ZmLoss1+ZmLoss2, data=pfish)
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
  dfish$Fish[(yst-3):(yen-3)] <- D[,cn]
  dfish <- na.omit(dfish)
  
  dmod <- lm(Fish ~ TP0+TP1+TP2 + TB0+TB1+TB2 + Det0+Det1+Det2 + 
               ZmLoss0+ZmLoss1+ZmLoss2, data=dfish)
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
  afish$Fish[(yst-3):(yen-3)] <- A[,cn]
  afish <- na.omit(afish)
  
  amod <- lm(Fish ~ TP0+TP1+TP2 + TB0+TB1+TB2 + Det0+Det1+Det2 + 
               ZmLoss0+ZmLoss1+ZmLoss2, data=afish)
  acombo <- dredge(amod, extra = c("R^2", F = function(x)
    summary(x)$fstatistic[[1]]))
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (acombo)[1]
  acoef[i, match(names(xx), colnames(acoef))] = xx
  y <- summary(get.models(acombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  apval[i, match(names(yy), colnames(apval))] = yy
  
} # LMEs

write.table(fcoef,paste0(datar,"LME_F_gam_mlr_coeffs_ALLdiv2SD_alllags3_v2_noint.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datar,"LME_P_gam_mlr_coeffs_ALLdiv2SD_alllags3_v2_noint.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datar,"LME_D_gam_mlr_coeffs_ALLdiv2SD_alllags3_v2_noint.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datar,"LME_A_gam_mlr_coeffs_ALLdiv2SD_alllags3_v2_noint.csv"),sep=",",row.names=F)

write.table(fpval,paste0(datar,"LME_F_gam_mlr_pvals_ALLdiv2SD_alllags3_v2_noint.csv"),sep=",",row.names=F)
write.table(ppval,paste0(datar,"LME_P_gam_mlr_pvals_ALLdiv2SD_alllags3_v2_noint.csv"),sep=",",row.names=F)
write.table(dpval,paste0(datar,"LME_D_gam_mlr_pvals_ALLdiv2SD_alllags3_v2_noint.csv"),sep=",",row.names=F)
write.table(apval,paste0(datar,"LME_A_gam_mlr_pvals_ALLdiv2SD_alllags3_v2_noint.csv"),sep=",",row.names=F)

