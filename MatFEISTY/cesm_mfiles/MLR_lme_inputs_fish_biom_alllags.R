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
#datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/FOSI/"
# TP <- read.csv(paste0(datap,"LMEs_anoms_div2SD_TP.csv"),sep=",",header = T,stringsAsFactors = F)
# TB <- read.csv(paste0(datap,"LMEs_anoms_div2SD_TB.csv"),sep=",",header = T,stringsAsFactors = F)
# Det <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Det.csv"),sep=",",header = T,stringsAsFactors = F)
# Zmeso <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Zmeso.csv"),sep=",",header = T,stringsAsFactors = F)
# ZmLoss <- read.csv(paste0(datap,"LMEs_anoms_div2SD_ZmLoss.csv"),sep=",",header = T,stringsAsFactors = F)
# S <- read.csv(paste0(datap,"LMEs_anoms_div2SD_S.csv"),sep=",",header = T,stringsAsFactors = F)
# M <- read.csv(paste0(datap,"LMEs_anoms_div2SD_M.csv"),sep=",",header = T,stringsAsFactors = F)
# L <- read.csv(paste0(datap,"LMEs_anoms_div2SD_L.csv"),sep=",",header = T,stringsAsFactors = F)
# FF <- read.csv(paste0(datap,"LMEs_anoms_div2SD_F.csv"),sep=",",header = T,stringsAsFactors = F)
# P <- read.csv(paste0(datap,"LMEs_anoms_div2SD_P.csv"),sep=",",header = T,stringsAsFactors = F)
# D <- read.csv(paste0(datap,"LMEs_anoms_div2SD_D.csv"),sep=",",header = T,stringsAsFactors = F)
# A <- read.csv(paste0(datap,"LMEs_anoms_div2SD_A.csv"),sep=",",header = T,stringsAsFactors = F)
# B <- read.csv(paste0(datap,"LMEs_anoms_div2SD_B.csv"),sep=",",header = T,stringsAsFactors = F)

datap <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/"
TP <- read.csv(paste0(datap,"LMEs_anoms_TP_trans.csv"),sep=",",header = T,stringsAsFactors = F)
TB <- read.csv(paste0(datap,"LMEs_anoms_TB_trans.csv"),sep=",",header = T,stringsAsFactors = F)
Det <- read.csv(paste0(datap,"LMEs_anoms_Det_trans.csv"),sep=",",header = T,stringsAsFactors = F)
Zmeso <- read.csv(paste0(datap,"LMEs_anoms_Zmeso_trans.csv"),sep=",",header = T,stringsAsFactors = F)
ZmLoss <- read.csv(paste0(datap,"LMEs_anoms_ZmLoss_trans.csv"),sep=",",header = T,stringsAsFactors = F)
FF <- read.csv(paste0(datap,"LMEs_anoms_F_trans.csv"),sep=",",header = T,stringsAsFactors = F)
P <- read.csv(paste0(datap,"LMEs_anoms_P_trans.csv"),sep=",",header = T,stringsAsFactors = F)
D <- read.csv(paste0(datap,"LMEs_anoms_D_trans.csv"),sep=",",header = T,stringsAsFactors = F)
A <- read.csv(paste0(datap,"LMEs_anoms_A_trans.csv"),sep=",",header = T,stringsAsFactors = F)
B <- read.csv(paste0(datap,"LMEs_anoms_B_trans.csv"),sep=",",header = T,stringsAsFactors = F)

datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"


cname <- names(A)
nlme <- length(cname)

### Divide by 2SD 
TP2 <- TP
TB2 <- TB
Det2 <- Det
Zmeso2 <- Zmeso
ZmLoss2 <- ZmLoss
A2 <- A
B2 <- B
D2 <- D
F2 <- FF
P2 <- P

for (i in 1:nlme) {
  TP2[,i] <- TP[,i]/(2*sd(TP[,i]))
  TB2[,i] <- TB[,i]/(2*sd(TB[,i]))
  Det2[,i] <- Det[,i]/(2*sd(Det[,i]))
  Zmeso2[,i] <- Zmeso[,i]/(2*sd(Zmeso[,i]))
  ZmLoss2[,i] <- ZmLoss[,i]/(2*sd(ZmLoss[,i]))
  A2[,i] <- A[,i]/(2*sd(A[,i]))
  B2[,i] <- B[,i]/(2*sd(B[,i]))
  D2[,i] <- D[,i]/(2*sd(D[,i]))
  F2[,i] <- FF[,i]/(2*sd(FF[,i]))
  P2[,i] <- P[,i]/(2*sd(P[,i]))
}

TP2$Year <- c(1:68)
TB2$Year <- c(1:68)
Det2$Year <- c(1:68)
Zmeso2$Year <- c(1:68)
ZmLoss2$Year <- c(1:68)
A2$Year <- c(1:68)
B2$Year <- c(1:68)
D2$Year <- c(1:68)
F2$Year <- c(1:68)
P2$Year <- c(1:68)

### MLR of drivers ---------------------------------------------------------
dname <- c('TP','TB','Zmeso','ZmLoss','Det')

#Loop over LMEs  1:nlme
for (i in 24:34) {
  cn <- cname[i]
  
  ### Drivers
  maxl <- 3
  yst <- 1+maxl
  yen <- 68+maxl
  
  drive <- data.frame(matrix(ncol = (5*4), nrow = yen))
  names(drive) <- c('TP0','TP1','TP2','TP3',
                    'TB0','TB1','TB2','TB3',
                    'Det0','Det1','Det2','Det3',
                    'Zmeso0','Zmeso1','Zmeso2','Zmeso3',
                    'ZmLoss0','ZmLoss1','ZmLoss2','ZmLoss3')
  
  drive$TP0[(yst-3):(yen-3)] <- TP2[,cn]
  drive$TP1[(yst-2):(yen-2)] <- TP2[,cn]
  drive$TP2[(yst-1):(yen-1)] <- TP2[,cn]
  drive$TP3[(yst-0):(yen-0)] <- TP2[,cn]
  
  drive$TB0[(yst-3):(yen-3)] <- TB2[,cn]
  drive$TB1[(yst-2):(yen-2)] <- TB2[,cn]
  drive$TB2[(yst-1):(yen-1)] <- TB2[,cn]
  drive$TB3[(yst-0):(yen-0)] <- TB2[,cn]
  
  drive$Det0[(yst-3):(yen-3)] <- Det2[,cn]
  drive$Det1[(yst-2):(yen-2)] <- Det2[,cn]
  drive$Det2[(yst-1):(yen-1)] <- Det2[,cn]
  drive$Det3[(yst-0):(yen-0)] <- Det2[,cn]
  
  drive$Zmeso0[(yst-3):(yen-3)] <- Zmeso2[,cn]
  drive$Zmeso1[(yst-2):(yen-2)] <- Zmeso2[,cn]
  drive$Zmeso2[(yst-1):(yen-1)] <- Zmeso2[,cn]
  drive$Zmeso3[(yst-0):(yen-0)] <- Zmeso2[,cn]
  
  drive$ZmLoss0[(yst-3):(yen-3)] <- ZmLoss2[,cn]
  drive$ZmLoss1[(yst-2):(yen-2)] <- ZmLoss2[,cn]
  drive$ZmLoss2[(yst-1):(yen-1)] <- ZmLoss2[,cn]
  drive$ZmLoss3[(yst-0):(yen-0)] <- ZmLoss2[,cn]
  
  
  ### F ----------------------------------------------------------------
  ffish <- drive
  ffish$Fish <- NA
  ffish$Fish[(yst-3):(yen-3)] <- F2[,cn]
  ffish <- na.omit(ffish)
  
  options(na.action = "na.fail")
  # fmod <- lm(Fish ~ TP0+TP1+TP2+TP3 + TB0+TB1+TB2+TB3 + Det0+Det1+Det2+Det3 + 
  #              Zmeso0+Zmeso1+Zmeso2+Zmeso3 + ZmLoss0+ZmLoss1+ZmLoss2+ZmLoss3 + 
  #              TP0*Zmeso0 + TP1*Zmeso1 + TP2*Zmeso2 + TP3*Zmeso3 + 
  #              TP0*ZmLoss0 + TP1*ZmLoss1 + TP2*ZmLoss2 + TP3*ZmLoss3 + 
  #              TB0*Det0 + TB1*Det1 + TB2*Det2 + TB3*Det3, data=ffish)
  fmod <- lm(Fish ~ TP0+TP1+TP2 + TB0+TB1+TB2 + Det0+Det1+Det2 + 
               Zmeso0+Zmeso1+Zmeso2 + ZmLoss0+ZmLoss1+ZmLoss2 + 
               TP0*Zmeso0 + TP1*Zmeso1 + TP2*Zmeso2 + 
               TP0*ZmLoss0 + TP1*ZmLoss1 + TP2*ZmLoss2 + 
               TB0*Det0 + TB1*Det1 + TB2*Det2, data=ffish)
  summary(fmod)
  fcombo <- dredge(fmod)
  
  ## Create arrays for coefficients & p-vals
  if (i==1) {
    fcoef <- data.frame(matrix(ncol = length(fcombo), nrow = nlme))
    names(fcoef) <- names(fcombo)
    fpval <- fcoef[,1:25]
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
  pfish$Fish[(yst-3):(yen-3)] <- P2[,cn]
  pfish <- na.omit(pfish)
  
  pmod <- lm(Fish ~ TP0+TP1+TP2 + TB0+TB1+TB2 + Det0+Det1+Det2 + 
               Zmeso0+Zmeso1+Zmeso2 + ZmLoss0+ZmLoss1+ZmLoss2 + 
               TP0*Zmeso0 + TP1*Zmeso1 + TP2*Zmeso2 + 
               TP0*ZmLoss0 + TP1*ZmLoss1 + TP2*ZmLoss2 + 
               TB0*Det0 + TB1*Det1 + TB2*Det2, data=pfish)
  pcombo <- dredge(pmod)
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
               Zmeso0+Zmeso1+Zmeso2 + ZmLoss0+ZmLoss1+ZmLoss2 + 
               TP0*Zmeso0 + TP1*Zmeso1 + TP2*Zmeso2 + 
               TP0*ZmLoss0 + TP1*ZmLoss1 + TP2*ZmLoss2 + 
               TB0*Det0 + TB1*Det1 + TB2*Det2, data=dfish)
  dcombo <- dredge(dmod)
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
               Zmeso0+Zmeso1+Zmeso2 + ZmLoss0+ZmLoss1+ZmLoss2 + 
               TP0*Zmeso0 + TP1*Zmeso1 + TP2*Zmeso2 + 
               TP0*ZmLoss0 + TP1*ZmLoss1 + TP2*ZmLoss2 + 
               TB0*Det0 + TB1*Det1 + TB2*Det2, data=afish)
  acombo <- dredge(amod)
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (acombo)[1]
  acoef[i, match(names(xx), colnames(acoef))] = xx
  y <- summary(get.models(acombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  apval[i, match(names(yy), colnames(apval))] = yy
  
  
  ### B ----------------------------------------------------------------
  bent <- drive
  bent$Bent <- NA
  bent$Bent[(yst-3):(yen-3)] <- B2[,cn]
  bent <- na.omit(bent)
  
  bmod <- lm(Bent ~ TP0+TP1+TP2 + TB0+TB1+TB2 + Det0+Det1+Det2 + 
               Zmeso0+Zmeso1+Zmeso2 + ZmLoss0+ZmLoss1+ZmLoss2 + 
               TP0*Zmeso0 + TP1*Zmeso1 + TP2*Zmeso2 + 
               TP0*ZmLoss0 + TP1*ZmLoss1 + TP2*ZmLoss2 + 
               TB0*Det0 + TB1*Det1 + TB2*Det2, data=bent)
  bcombo <- dredge(bmod)
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (bcombo)[1]
  bcoef[i, match(names(xx), colnames(bcoef))] = xx
  y <- summary(get.models(bcombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  bpval[i, match(names(yy), colnames(bpval))] = yy
  
  
} # LMEs

write.table(fcoef,paste0(datar,"LME2434_F_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datar,"LME2434_P_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datar,"LME2434_D_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datar,"LME2434_A_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(bcoef,paste0(datar,"LME2434_B_mlr_coeffs_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)

write.table(fpval,paste0(datar,"LME2434_F_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(ppval,paste0(datar,"LME2434_P_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(dpval,paste0(datar,"LME2434_D_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(apval,paste0(datar,"LME2434_A_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)
write.table(bpval,paste0(datar,"LME2434_B_mlr_pvals_ALLdiv2SD_alllags.csv"),sep=",",row.names=F)

