# Evaluate R2 or frac var explained by best MLRs

rm(list=ls())

library("MuMIn")

source(file = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/gammaGLM_fns20200504.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/corrs/"

### ------------------------------------------------------------
# load data
#datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/FOSI/"
datap <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/"
TP <- read.csv(paste0(datap,"LMEs_anoms_TP_trans.csv"),sep=",",header = T,stringsAsFactors = F)
TB <- read.csv(paste0(datap,"LMEs_anoms_TB_trans.csv"),sep=",",header = T,stringsAsFactors = F)
Det <- read.csv(paste0(datap,"LMEs_anoms_Det_trans.csv"),sep=",",header = T,stringsAsFactors = F)
ZmLoss <- read.csv(paste0(datap,"LMEs_anoms_ZmLoss_trans.csv"),sep=",",header = T,stringsAsFactors = F)

F2 <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Fnu.csv"),sep=",",header = T,stringsAsFactors = F)
P2 <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Pnu.csv"),sep=",",header = T,stringsAsFactors = F)
D2 <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Dnu.csv"),sep=",",header = T,stringsAsFactors = F)
A2 <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Anu.csv"),sep=",",header = T,stringsAsFactors = F)

#datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"
datar <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/"

fcoef <- read.csv(paste0(datar,"LME_Fnu_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef <- read.csv(paste0(datar,"LME_Pnu_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef <- read.csv(paste0(datar,"LME_Dnu_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
acoef <- read.csv(paste0(datar,"LME_Anu_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)

fcoef$Rcorr <- NA
pcoef$Rcorr <- NA
dcoef$Rcorr <- NA
acoef$Rcorr <- NA

cname <- names(A2)
nlme <- length(cname)-1

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
names(TP2) <- names(A2)[2:nlme]
names(TB2) <- names(A2)[2:nlme]
names(Det2) <- names(A2)[2:nlme]
names(ZmLoss2) <- names(A2)[2:nlme]

TP2$Year <- c(1:68)
TB2$Year <- c(1:68)
Det2$Year <- c(1:68)
ZmLoss2$Year <- c(1:68)


### MLR of drivers ---------------------------------------------------------
dname <- c('TP','TB','ZmLoss','Det')

#Loop over LMEs  1:nlme
for (i in 1:nlme) {
  cn <- cname[i]
  
  ### Drivers
  maxl <- 3
  yst <- 1+maxl
  yen <- 68+maxl
  
  drive <- data.frame(matrix(ncol = (4*2), nrow = yen))
  names(drive) <- c('TP0','TP1',
                    'TB0','TB1',
                    'Det0','Det1',
                    'ZmLoss0','ZmLoss1')
  
  drive$TP0[(yst-3):(yen-3)] <- TP2[,cn]
  drive$TP1[(yst-2):(yen-2)] <- TP2[,cn]
  
  drive$TB0[(yst-3):(yen-3)] <- TB2[,cn]
  drive$TB1[(yst-2):(yen-2)] <- TB2[,cn]
  
  drive$Det0[(yst-3):(yen-3)] <- Det2[,cn]
  drive$Det1[(yst-2):(yen-2)] <- Det2[,cn]
  
  drive$ZmLoss0[(yst-3):(yen-3)] <- ZmLoss2[,cn]
  drive$ZmLoss1[(yst-2):(yen-2)] <- ZmLoss2[,cn]
  
  drive$Det0.TB0 <- drive$Det0 * drive$TB0
  drive$Det0.ZmLoss0 <- drive$Det0 * drive$ZmLoss0
  drive$Det1.TB1 <- drive$Det1 * drive$TB1
  drive$Det1.ZmLoss1 <- drive$Det1 * drive$ZmLoss1
  drive$TP0.ZmLoss0 <- drive$TP0 * drive$ZmLoss0
  drive$TP1.ZmLoss1 <- drive$TP1 * drive$ZmLoss1
  drive$X.Intercept. <- 1
  
  ### F ----------------------------------------------------------------
  # best model
  fmod <- fcoef[i,1:15]
  
  fcols <- names(fmod)
  fdrive <- drive
  fdrive <- na.omit(fdrive)
  fdrive <- fdrive[,fcols]
  
  ffish <- drive
  ffish$Fish <- NA
  ffish$Fish[(yst-3):(yen-3)] <- F2[,cn]
  ffish <- na.omit(ffish)
  
  fmod2 <- data.frame(matrix(ncol = length(fdrive), nrow = nrow(fdrive)))
  for (j in 1:nrow(fdrive)) {
    fmod2[j,] <- fmod
  }
  fpred <- fmod2 * fdrive
  fpred2 <- rowSums(fpred,na.rm=T)
  
  fcoef$Rcorr[i] <- cor(ffish$Fish,fpred2)
  
  ### P ----------------------------------------------------------------
  # best model
  pmod <- pcoef[i,1:15]
  
  pcols <- names(pmod)
  pdrive <- drive
  pdrive <- na.omit(pdrive)
  pdrive <- pdrive[,pcols]
  
  pfish <- drive
  pfish$Fish <- NA
  pfish$Fish[(yst-3):(yen-3)] <- P2[,cn]
  pfish <- na.omit(pfish)
  
  pmod2 <- data.frame(matrix(ncol = length(pdrive), nrow = nrow(pdrive)))
  for (j in 1:nrow(pdrive)) {
    pmod2[j,] <- pmod
  }
  ppred <- pmod2 * pdrive
  ppred2 <- rowSums(ppred,na.rm=T)
  
  pcoef$Rcorr[i] <- cor(ffish$Fish,ppred2)
  
  
  ### D ----------------------------------------------------------------
  # best model
  dmod <- dcoef[i,1:15]
  
  dcols <- names(dmod)
  ddrive <- drive
  ddrive <- na.omit(ddrive)
  ddrive <- ddrive[,dcols]
  
  dfish <- drive
  dfish$Fish <- NA
  dfish$Fish[(yst-3):(yen-3)] <- D2[,cn]
  dfish <- na.omit(dfish)
  
  dmod2 <- data.frame(matrix(ncol = length(ddrive), nrow = nrow(ddrive)))
  for (j in 1:nrow(ddrive)) {
    dmod2[j,] <- dmod
  }
  dpred <- dmod2 * ddrive
  dpred2 <- rowSums(dpred,na.rm=T)
  
  dcoef$Rcorr[i] <- cor(ffish$Fish,dpred2)
  
  ### A ----------------------------------------------------------------
  # best model
  amod <- acoef[i,1:15]
  
  acols <- names(amod)
  adrive <- drive
  adrive <- na.omit(adrive)
  adrive <- adrive[,acols]
  
  afish <- drive
  afish$Fish <- NA
  afish$Fish[(yst-3):(yen-3)] <- A2[,cn]
  afish <- na.omit(afish)
  
  amod2 <- data.frame(matrix(ncol = length(adrive), nrow = nrow(adrive)))
  for (j in 1:nrow(adrive)) {
    amod2[j,] <- amod
  }
  apred <- amod2 * adrive
  apred2 <- rowSums(apred,na.rm=T)
  
  acoef$Rcorr[i] <- cor(ffish$Fish,apred2)
  
  
} # LMEs

### One table with all corrs
Rcorr <- as.data.frame(fcoef$LME)
names(Rcorr) <- "LME"
Rcorr$Fcorr <- fcoef$Rcorr
Rcorr$Pcorr <- pcoef$Rcorr
Rcorr$Dcorr <- dcoef$Rcorr
Rcorr$Acorr <- acoef$Rcorr

write.table(Rcorr,paste0(datap,"LME_mlr_nu_Rcorrs_ALLdiv2SD_alllags_v2.csv"),sep=",",row.names=F)

write.table(fcoef,paste0(datap,"LME_Fnu_mlr_coeffs_Rcorr_ALLdiv2SD_alllags_v2.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datap,"LME_Pnu_mlr_coeffs_Rcorr_ALLdiv2SD_alllags_v2.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datap,"LME_Dnu_mlr_coeffs_Rcorr_ALLdiv2SD_alllags_v2.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datap,"LME_Anu_mlr_coeffs_Rcorr_ALLdiv2SD_alllags_v2.csv"),sep=",",row.names=F)

