# Mult linear regression of forcing on fish biomass

rm(list=ls())

library("MuMIn")


source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/corrs/"

### ------------------------------------------------------------
# load data
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/FOSI/"
TP <- read.csv(paste0(datap,"LMEs_anoms_div2SD_TP.csv"),sep=",",header = T,stringsAsFactors = F)
TB <- read.csv(paste0(datap,"LMEs_anoms_div2SD_TB.csv"),sep=",",header = T,stringsAsFactors = F)
Det <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Det.csv"),sep=",",header = T,stringsAsFactors = F)
Zmeso <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Zmeso.csv"),sep=",",header = T,stringsAsFactors = F)
ZmLoss <- read.csv(paste0(datap,"LMEs_anoms_div2SD_ZmLoss.csv"),sep=",",header = T,stringsAsFactors = F)
# S <- read.csv(paste0(datap,"LMEs_anoms_div2SD_S.csv"),sep=",",header = T,stringsAsFactors = F)
# M <- read.csv(paste0(datap,"LMEs_anoms_div2SD_M.csv"),sep=",",header = T,stringsAsFactors = F)
# L <- read.csv(paste0(datap,"LMEs_anoms_div2SD_L.csv"),sep=",",header = T,stringsAsFactors = F)
FF <- read.csv(paste0(datap,"LMEs_anoms_div2SD_F.csv"),sep=",",header = T,stringsAsFactors = F)
P <- read.csv(paste0(datap,"LMEs_anoms_div2SD_P.csv"),sep=",",header = T,stringsAsFactors = F)
D <- read.csv(paste0(datap,"LMEs_anoms_div2SD_D.csv"),sep=",",header = T,stringsAsFactors = F)
A <- read.csv(paste0(datap,"LMEs_anoms_div2SD_A.csv"),sep=",",header = T,stringsAsFactors = F)
B <- read.csv(paste0(datap,"LMEs_anoms_div2SD_B.csv"),sep=",",header = T,stringsAsFactors = F)

datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"
Flag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_F.csv"),sep=",",header = T,stringsAsFactors = F)
Plag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_P.csv"),sep=",",header = T,stringsAsFactors = F)
Dlag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_D.csv"),sep=",",header = T,stringsAsFactors = F)
Alag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_A.csv"),sep=",",header = T,stringsAsFactors = F)
Blag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_B.csv"),sep=",",header = T,stringsAsFactors = F)

### MLR of drivers ---------------------------------------------------------
dname <- c('TP','TB','Zmeso','ZmLoss','Det')

#Start with US LMEs
cname <- names(A)[2:64]
nlme <- length(cname)

#Loop over LMEs 
for (i in 1:nlme) {
  cn <- cname[i]
  
  ### Drivers
  drive <- TP[,c('Year',cn)]
  names(drive) <- c('Year','TP')
  drive$TB <- TB[,cn]
  drive$Det <- Det[,cn]
  drive$Zmeso <- Zmeso[,cn]
  drive$ZmLoss <- ZmLoss[,cn]
  
  
  ### F ----------------------------------------------------------------
  maxf <- max(Flag[i,])
  yst <- 1+maxf
  yen <- 68+maxf
  ffish <- data.frame(matrix(ncol = 6, nrow = yen))
  names(ffish) <- c('Fish','TP','TB','Det','Zmeso','ZmLoss')
  
  ffish$TP[(yst-Flag$TP[i]):(yen-Flag$TP[i])] <- drive$TP 
  ffish$TB[(yst-Flag$TB[i]):(yen-Flag$TB[i])] <- drive$TB 
  ffish$Det[(yst-Flag$Det[i]):(yen-Flag$Det[i])] <- drive$Det
  ffish$Zmeso[(yst-Flag$Zmeso[i]):(yen-Flag$Zmeso[i])] <- drive$Zmeso 
  ffish$ZmLoss[(yst-Flag$ZmLoss[i]):(yen-Flag$ZmLoss[i])] <- drive$ZmLoss
  ffish$Fish[(yst):(yen)] <- FF[,cn]
  ffish <- na.omit(ffish)
  
  options(na.action = "na.fail")
  fmod <- lm(Fish ~ TP + TB + Det + Zmeso + ZmLoss + TP*Zmeso + TP*ZmLoss + TB*Det, data=ffish)
  summary(fmod)
  fcombo <- dredge(fmod)
  ## Create arrays for coefficients & p-vals
  if (i==1) {
    fcoef <- data.frame(matrix(ncol = length(fcombo), nrow = nlme))
    names(fcoef) <- names(fcombo)
    fpval <- fcoef[,1:9]
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
  
  