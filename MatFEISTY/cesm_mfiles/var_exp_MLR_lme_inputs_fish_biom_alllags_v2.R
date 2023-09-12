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
  
  
  ### F ----------------------------------------------------------------
  fdrive <- drive
  fdrive <- na.omit(fdrive)
  
  ffish <- drive
  ffish$Fish <- NA
  ffish$Fish[(yst-3):(yen-3)] <- F2[,cn]
  ffish <- na.omit(ffish)
  
  # best model
  fmod <- fcoef[i,1:22]
  fmod2 <- fmod[,!is.na(fmod)]
  m <- glm( vs ~ LH + cvF + LH:cvF, family = Gamma( link = 'log' ), data = dds )
  # m <- glm( I( 1000 * vs ) ~ LH + cvF + LH:cvF, family = Gamma( link = 'log' ), data = dds ) # altering the scale does not affect the estimates [beyond the intercept]
  
  # deviances and r squares 
  s <- summary( m, dispersion = 1 / gamma.shape( m )$alpha )
  
  #! check this summary
  mAlt <- glm( vs  ~ 0 + LH + LH:cvF, family = Gamma( link = 'log' ), data = dds )
  sAlt <- summary( mAlt, dispersion = 1 / gamma.shape( mAlt )$alpha )
  
  print( sAlt, digits = 3 )
  xtable( sAlt, digits = 2)
  
  # fraction of deviance explained
  fdevexpl( m )
  fdevexpl( mAlt ) # same
  
  # pseudo r square - we need to rescale the dependent vairable
  m <- glm( I( vs * 1000 ) ~ LH + cvF + LH:cvF, family = Gamma( link = 'log' ), data = dds )
  mAlt <- glm( I( vs * 1000 )  ~ 0 + LH + LH:cvF, family = Gamma( link = 'log' ), data = dds )
  fgNagelkerke( m )
  fgNagelkerke( mAlt )# bit different but minor ...
  
  
  
  ### P ----------------------------------------------------------------
  pfish <- drive
  pfish$Fish <- NA
  pfish$Fish[(yst-3):(yen-3)] <- P2[,cn]
  pfish <- na.omit(pfish)
  
  
  
  
  ### D ----------------------------------------------------------------
  dfish <- drive
  dfish$Fish <- NA
  dfish$Fish[(yst-3):(yen-3)] <- D2[,cn]
  dfish <- na.omit(dfish)
  
  
  ### A ----------------------------------------------------------------
  afish <- drive
  afish$Fish <- NA
  afish$Fish[(yst-3):(yen-3)] <- A2[,cn]
  afish <- na.omit(afish)
  
  
  
  ### B ----------------------------------------------------------------
  bent <- drive
  bent$Bent <- NA
  bent$Bent[(yst-3):(yen-3)] <- B2[,cn]
  bent <- na.omit(bent)
  
  
  
} # LMEs

