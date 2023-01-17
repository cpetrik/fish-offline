# Find patterns in forcing-fish correlations

rm(list=ls())


source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/"

### ------------------------------------------------------------
# load data
datam <- "/Volumes/MIP/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/"
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"
CH <- read.csv(paste0(datap,"CHK_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
BS <- read.csv(paste0(datap,"EBS_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
AK <- read.csv(paste0(datap,"GAK_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
HI <- read.csv(paste0(datap,"HI_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
CC <- read.csv(paste0(datap,"CCE_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
MX <- read.csv(paste0(datap,"GMX_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
SE <- read.csv(paste0(datap,"SE_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
NE <- read.csv(paste0(datap,"NE_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)


### LM of forcing ---------------------------------------------------------
cnam <- c('Lag','Type','coef', 'se', 'p','r2')

#Loop over forcing and climate indices
#input forcing
for (j in 2:6) {
  #fish vars
  driver <- names(AK)[j]
  
  AKtab <- data.frame(matrix(ncol = 6, nrow = 48))
  colnames(AKtab) <- cnam
  BStab <- data.frame(matrix(ncol = 6, nrow = 48))
  colnames(BStab) <- cnam
  CCtab <- data.frame(matrix(ncol = 6, nrow = 48))
  colnames(CCtab) <- cnam
  CHtab <- data.frame(matrix(ncol = 6, nrow = 48))
  colnames(CHtab) <- cnam
  HItab <- data.frame(matrix(ncol = 6, nrow = 48))
  colnames(HItab) <- cnam
  MXtab <- data.frame(matrix(ncol = 6, nrow = 48))
  colnames(MXtab) <- cnam
  NEtab <- data.frame(matrix(ncol = 6, nrow = 48))
  colnames(NEtab) <- cnam
  SEtab <- data.frame(matrix(ncol = 6, nrow = 48))
  colnames(SEtab) <- cnam
  
  n <- 0
  for (i in 7:14) {
    fish <- names(AK)[i]
    
    ### AK
    rdat <- as.data.frame(AK[,j])
    names(rdat) <- "inp"
    rdat$inp <- rdat$inp/(2*sd(rdat$inp, na.rm = T))
    adat <- rdat
    adat$out <- AK[,i]
    adat <- na.omit(adat)
    anlen <- dim(adat)[1]
    
    ### BS
    rdat <- as.data.frame(BS[,j])
    names(rdat) <- "inp"
    rdat$inp <- rdat$inp/(2*sd(rdat$inp, na.rm = T))
    bdat <- rdat
    bdat$out <- BS[,i]
    bdat <- na.omit(bdat)
    bnlen <- dim(bdat)[1]
    
    ### CCE
    rdat <- as.data.frame(CC[,j])
    names(rdat) <- "inp"
    rdat$inp <- rdat$inp/(2*sd(rdat$inp, na.rm = T))
    cdat <- rdat
    cdat$out <- CC[,i]
    cdat <- na.omit(cdat)
    cnlen <- dim(cdat)[1]
    
    ##CHK
    rdat <- as.data.frame(CH[,j])
    names(rdat) <- "inp"
    rdat$inp <- rdat$inp/(2*sd(rdat$inp, na.rm = T))
    kdat <- rdat
    kdat$out <- CH[,i]
    kdat <- na.omit(kdat)
    knlen <- dim(kdat)[1]
    
    ### HI
    rdat <- as.data.frame(HI[,j])
    names(rdat) <- "inp"
    rdat$inp <- rdat$inp/(2*sd(rdat$inp, na.rm = T))
    hdat <- rdat
    hdat$out <- HI[,i]
    hdat <- na.omit(hdat)
    hnlen <- dim(hdat)[1]
    
    ### MX
    rdat <- as.data.frame(MX[,j])
    names(rdat) <- "inp"
    rdat$inp <- rdat$inp/(2*sd(rdat$inp, na.rm = T))
    mdat <- rdat
    mdat$out <- MX[,i]
    mdat <- na.omit(mdat)
    mnlen <- dim(mdat)[1]
    
    ### NE
    rdat <- as.data.frame(NE[,j])
    names(rdat) <- "inp"
    rdat$inp <- rdat$inp/(2*sd(rdat$inp, na.rm = T))
    ndat <- rdat
    ndat$out <- NE[,i]
    ndat <- na.omit(ndat)
    nnlen <- dim(ndat)[1]
    
    ### SE
    rdat <- as.data.frame(SE[,j])
    names(rdat) <- "inp"
    rdat$inp <- rdat$inp/(2*sd(rdat$inp, na.rm = T))
    sdat <- rdat
    sdat$out <- SE[,i]
    sdat <- na.omit(sdat)
    snlen <- dim(sdat)[1]
    
    for (k in 0:5) {
      n <- n+1
      
      AKlm <- lm(out[1:(anlen-k)] ~ inp[(k+1):anlen],adat)
      Alag   = summary(AKlm)
      AKtab[n,1] <- k
      AKtab[n,2] <- fish
      AKtab[n,3] <- Alag$coefficients[2,1]
      AKtab[n,4] <- Alag$coefficients[2,2]
      AKtab[n,5] <- Alag$coefficients[2,4]
      AKtab[n,6] <- Alag$r.squared
    
      BSlm <- lm(out[1:(bnlen-k)] ~ inp[(k+1):bnlen],bdat)
      Blag   = summary(BSlm)
      BStab[n,1] <- k
      BStab[n,2] <- fish
      BStab[n,3] <- Blag$coefficients[2,1]
      BStab[n,4] <- Blag$coefficients[2,2]
      BStab[n,5] <- Blag$coefficients[2,4]
      BStab[n,6] <- Blag$r.squared
    
      CClm <- lm(out[1:(cnlen-k)] ~ inp[(k+1):cnlen],cdat)
      Clag   = summary(CClm)
      CCtab[n,1] <- k
      CCtab[n,2] <- fish
      CCtab[n,3] <- Clag$coefficients[2,1]
      CCtab[n,4] <- Clag$coefficients[2,2]
      CCtab[n,5] <- Clag$coefficients[2,4]
      CCtab[n,6] <- Clag$r.squared
    
      CHlm <- lm(out[1:(knlen-k)] ~ inp[(k+1):knlen],kdat)
      Klag   = summary(CHlm)
      CHtab[n,1] <- k
      CHtab[n,2] <- fish
      CHtab[n,3] <- Klag$coefficients[2,1]
      CHtab[n,4] <- Klag$coefficients[2,2]
      CHtab[n,5] <- Klag$coefficients[2,4]
      CHtab[n,6] <- Klag$r.squared
    
      HIlm <- lm(out[1:(hnlen-k)] ~ inp[(k+1):hnlen],hdat)
      Hlag   = summary(HIlm)
      HItab[n,1] <- k
      HItab[n,2] <- fish
      HItab[n,3] <- Hlag$coefficients[2,1]
      HItab[n,4] <- Hlag$coefficients[2,2]
      HItab[n,5] <- Hlag$coefficients[2,4]
      HItab[n,6] <- Hlag$r.squared
    
      MXlm <- lm(out[1:(mnlen-k)] ~ inp[(k+1):mnlen],mdat)
      Mlag   = summary(MXlm)
      MXtab[n,1] <- k
      MXtab[n,2] <- fish
      MXtab[n,3] <- Mlag$coefficients[2,1]
      MXtab[n,4] <- Mlag$coefficients[2,2]
      MXtab[n,5] <- Mlag$coefficients[2,4]
      MXtab[n,6] <- Mlag$r.squared
    
      NElm <- lm(out[1:(nnlen-k)] ~ inp[(k+1):nnlen],ndat)
      Nlag   = summary(NElm)
      NEtab[n,1] <- k
      NEtab[n,2] <- fish
      NEtab[n,3] <- Nlag$coefficients[2,1]
      NEtab[n,4] <- Nlag$coefficients[2,2]
      NEtab[n,5] <- Nlag$coefficients[2,4]
      NEtab[n,6] <- Nlag$r.squared
    
      SElm <- lm(out[1:(snlen-k)] ~ inp[(k+1):snlen],sdat)
      Slag   = summary(SElm)
      SEtab[n,1] <- k
      SEtab[n,2] <- fish
      SEtab[n,3] <- Slag$coefficients[2,1]
      SEtab[n,4] <- Slag$coefficients[2,2]
      SEtab[n,5] <- Slag$coefficients[2,4]
      SEtab[n,6] <- Slag$r.squared
    }
    
  } # inputs & fish
  write.table(AKtab,paste0(datap,"GAK_regress_",driver,"_div2SD_melt.csv"),sep=",",row.names=T)
  write.table(BStab,paste0(datap,"EBS_regress_",driver,"_div2SD_melt.csv"),sep=",",row.names=T)
  write.table(CCtab,paste0(datap,"CCE_regress_",driver,"_div2SD_melt.csv"),sep=",",row.names=T)
  write.table(CHtab,paste0(datap,"CHK_regress_",driver,"_div2SD_melt.csv"),sep=",",row.names=T)
  write.table(HItab,paste0(datap,"HI_regress_",driver,"_div2SD_melt.csv"),sep=",",row.names=T)
  write.table(MXtab,paste0(datap,"MX_regress_",driver,"_div2SD_melt.csv"),sep=",",row.names=T)
  write.table(NEtab,paste0(datap,"NE_regress_",driver,"_div2SD_melt.csv"),sep=",",row.names=T)
  write.table(SEtab,paste0(datap,"SE_regress_",driver,"_div2SD_melt.csv"),sep=",",row.names=T)
  
}


