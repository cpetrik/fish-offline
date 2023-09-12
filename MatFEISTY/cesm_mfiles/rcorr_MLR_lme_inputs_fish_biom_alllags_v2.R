# Evaluate R2 or frac var explained by best MLRs

rm(list=ls())

library("MuMIn")

source(file = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/gammaGLM_fns20200504.r")
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
ZmLoss <- read.csv(paste0(datap,"LMEs_anoms_ZmLoss_trans.csv"),sep=",",header = T,stringsAsFactors = F)
FF <- read.csv(paste0(datap,"LMEs_anoms_F_trans.csv"),sep=",",header = T,stringsAsFactors = F)
P <- read.csv(paste0(datap,"LMEs_anoms_P_trans.csv"),sep=",",header = T,stringsAsFactors = F)
D <- read.csv(paste0(datap,"LMEs_anoms_D_trans.csv"),sep=",",header = T,stringsAsFactors = F)
A <- read.csv(paste0(datap,"LMEs_anoms_A_trans.csv"),sep=",",header = T,stringsAsFactors = F)
B <- read.csv(paste0(datap,"LMEs_anoms_B_trans.csv"),sep=",",header = T,stringsAsFactors = F)

datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"

fcoef <- read.csv(paste0(datap,"LME_F_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
pcoef <- read.csv(paste0(datap,"LME_P_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
dcoef <- read.csv(paste0(datap,"LME_D_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
acoef <- read.csv(paste0(datap,"LME_A_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)
bcoef <- read.csv(paste0(datap,"LME_B_mlr_coeffs_ALLdiv2SD_alllags_v2.csv"),sep=",",header = T,stringsAsFactors = F)

fcoef$Rcorr <- NA
pcoef$Rcorr <- NA
dcoef$Rcorr <- NA
acoef$Rcorr <- NA
bcoef$Rcorr <- NA

cname <- names(A)
nlme <- length(cname)

### Divide by 2SD 
TP2 <- TP
TB2 <- TB
Det2 <- Det
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
ZmLoss2$Year <- c(1:68)
A2$Year <- c(1:68)
B2$Year <- c(1:68)
D2$Year <- c(1:68)
F2$Year <- c(1:68)
P2$Year <- c(1:68)

### MLR of drivers ---------------------------------------------------------
dname <- c('TP','TB','ZmLoss','Det')

#Loop over LMEs  1:nlme
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
  
  drive$Det0.TB0 <- drive$Det0 * drive$TB0
  drive$Det0.ZmLoss0 <- drive$Det0 * drive$ZmLoss0
  drive$Det1.TB1 <- drive$Det1 * drive$TB1
  drive$Det1.ZmLoss1 <- drive$Det1 * drive$ZmLoss1
  drive$Det2.TB2 <- drive$Det2 * drive$TB2
  drive$Det2.ZmLoss2 <- drive$Det2 * drive$ZmLoss2
  drive$TP0.ZmLoss0 <- drive$TP0 * drive$ZmLoss0
  drive$TP1.ZmLoss1 <- drive$TP1 * drive$ZmLoss1
  drive$TP2.ZmLoss2 <- drive$TP2 * drive$ZmLoss2
  drive$X.Intercept. <- 1
  
  ### F ----------------------------------------------------------------
  # best model
  fmod <- fcoef[i,1:22]
  
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
  pmod <- pcoef[i,1:22]
  
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
  dmod <- dcoef[i,1:22]
  
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
  amod <- acoef[i,1:22]
  
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
  
  ### B ----------------------------------------------------------------
  # best model
  bmod <- bcoef[i,1:22]
  
  bcols <- names(bmod)
  bdrive <- drive
  bdrive <- na.omit(bdrive)
  bdrive <- bdrive[,bcols]
  
  bent <- drive
  bent$Bent <- NA
  bent$Bent[(yst-3):(yen-3)] <- B2[,cn]
  bent <- na.omit(bent)
  
  bmod2 <- data.frame(matrix(ncol = length(bdrive), nrow = nrow(bdrive)))
  for (j in 1:nrow(bdrive)) {
    bmod2[j,] <- bmod
  }
  bpred <- bmod2 * bdrive
  bpred2 <- rowSums(bpred,na.rm=T)
  
  bcoef$Rcorr[i] <- cor(ffish$Fish,bpred2)
  
} # LMEs

### One table with all corrs
Rcorr <- as.data.frame(fcoef$LME)
names(Rcorr) <- "LME"
Rcorr$Fcorr <- fcoef$Rcorr
Rcorr$Pcorr <- pcoef$Rcorr
Rcorr$Dcorr <- dcoef$Rcorr
Rcorr$Acorr <- acoef$Rcorr
Rcorr$Bcorr <- bcoef$Rcorr

write.table(Rcorr,paste0(datap,"LME_mlr_biom_Rcorrs_ALLdiv2SD_alllags_v2.csv"),sep=",",row.names=F)

write.table(fcoef,paste0(datap,"LME_F_mlr_coeffs_Rcorr_ALLdiv2SD_alllags_v2.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datap,"LME_P_mlr_coeffs_Rcorr_ALLdiv2SD_alllags_v2.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datap,"LME_D_mlr_coeffs_Rcorr_ALLdiv2SD_alllags_v2.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datap,"LME_A_mlr_coeffs_Rcorr_ALLdiv2SD_alllags_v2.csv"),sep=",",row.names=F)
write.table(bcoef,paste0(datap,"LME_B_mlr_coeffs_Rcorr_ALLdiv2SD_alllags_v2.csv"),sep=",",row.names=F)

