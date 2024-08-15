# Find patterns in climate correlations

rm(list=ls())

library(modelsummary)


source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/"

### ------------------------------------------------------------
# load data
apath <- '/Users/cpetrik/Dropbox/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
C <- read.csv(paste0(apath, "Climate_anomalies_annual_means.csv"),sep=",",header = T,stringsAsFactors = F)

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


### LM of climate ---------------------------------------------------------
cnam <- c('L0coef', 'L0se', 'L0p','L0r2',
          'L1coef', 'L1se', 'L1p','L1r2',
          'L2coef', 'L2se', 'L2p','L2r2',
          'L3coef', 'L3se', 'L3p','L3r2',
          'L4coef', 'L4se', 'L4p','L4r2',
          'L5coef', 'L5se', 'L5p','L5r2')
#Loop over forcing and climate indices
#climate
for (j in 1:11) {
  #forcing and fish vars
  clim <- names(C)[j]
  
  AKtab <- data.frame(matrix(ncol = 24, nrow = 13))
  colnames(AKtab) <- cnam
  rownames(AKtab) <- names(AK)[2:14]
  BStab <- data.frame(matrix(ncol = 24, nrow = 13))
  colnames(BStab) <- cnam
  rownames(BStab) <- names(BS)[2:14]
  CCtab <- data.frame(matrix(ncol = 24, nrow = 13))
  colnames(CCtab) <- cnam
  rownames(CCtab) <- names(CC)[2:14]
  CHtab <- data.frame(matrix(ncol = 24, nrow = 13))
  colnames(CHtab) <- cnam
  rownames(CHtab) <- names(CH)[2:14]
  HItab <- data.frame(matrix(ncol = 24, nrow = 13))
  colnames(HItab) <- cnam
  rownames(HItab) <- names(HI)[2:14]
  MXtab <- data.frame(matrix(ncol = 24, nrow = 13))
  colnames(MXtab) <- cnam
  rownames(MXtab) <- names(MX)[2:14]
  NEtab <- data.frame(matrix(ncol = 24, nrow = 13))
  colnames(NEtab) <- cnam
  rownames(NEtab) <- names(NE)[2:14]
  SEtab <- data.frame(matrix(ncol = 24, nrow = 13))
  colnames(SEtab) <- cnam
  rownames(SEtab) <- names(SE)[2:14]
  
  for (i in 2:14) {
    driver <- names(AK)[i]
    rdat <- as.data.frame(C[,j])
    names(rdat) <- "inp"
    rdat$inp <- rdat$inp/(2*sd(rdat$inp, na.rm = T))
    
    ### AK
    adat <- rdat
    adat$out <- AK[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    for (k in 0:5) {
      AKlm <- lm(out[1:(nlen-k)] ~ inp[(k+1):nlen],adat)
      
      Alag   = summary(AKlm)
      
      # Want a table for an LME and driver, rows are forcing and fish, cols are lag
      AKtab[i-1,(4*k)+1] <- Alag$coefficients[2,1]
      AKtab[i-1,(4*k)+2] <- Alag$coefficients[2,2]
      AKtab[i-1,(4*k)+3] <- Alag$coefficients[2,4]
      AKtab[i-1,(4*k)+4] <- Alag$r.squared
    }   
    
    ### BS
    adat <- rdat
    adat$out <- BS[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    for (k in 0:5) {
      BSlm <- lm(out[1:(nlen-k)] ~ inp[(k+1):nlen],adat)
      
      Blag   = summary(BSlm)
      
      BStab[i-1,(4*k)+1] <- Blag$coefficients[2,1]
      BStab[i-1,(4*k)+2] <- Blag$coefficients[2,2]
      BStab[i-1,(4*k)+3] <- Blag$coefficients[2,4]
      BStab[i-1,(4*k)+4] <- Blag$r.squared
    }   
    
    ### CC
    adat <- rdat
    adat$out <- CC[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    for (k in 0:5) {
      CClm <- lm(out[1:(nlen-k)] ~ inp[(k+1):nlen],adat)
      
      Clag   = summary(CClm)
      
      CCtab[i-1,(4*k)+1] <- Clag$coefficients[2,1]
      CCtab[i-1,(4*k)+2] <- Clag$coefficients[2,2]
      CCtab[i-1,(4*k)+3] <- Clag$coefficients[2,4]
      CCtab[i-1,(4*k)+4] <- Clag$r.squared
    }  
    
    ### CH
    adat <- rdat
    adat$out <- CH[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    for (k in 0:5) {
      CHlm <- lm(out[1:(nlen-k)] ~ inp[(k+1):nlen],adat)
      
      Klag   = summary(CHlm)
      
      CHtab[i-1,(4*k)+1] <- Klag$coefficients[2,1]
      CHtab[i-1,(4*k)+2] <- Klag$coefficients[2,2]
      CHtab[i-1,(4*k)+3] <- Klag$coefficients[2,4]
      CHtab[i-1,(4*k)+4] <- Klag$r.squared
    }   
    
    ### HI
    adat <- rdat
    adat$out <- HI[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    for (k in 0:5) {
      HIlm <- lm(out[1:(nlen-k)] ~ inp[(k+1):nlen],adat)
      
      Hlag   = summary(HIlm)
      
      HItab[i-1,(4*k)+1] <- Hlag$coefficients[2,1]
      HItab[i-1,(4*k)+2] <- Hlag$coefficients[2,2]
      HItab[i-1,(4*k)+3] <- Hlag$coefficients[2,4]
      HItab[i-1,(4*k)+4] <- Hlag$r.squared
    }  
    
    ### MX
    adat <- rdat
    adat$out <- MX[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    for (k in 0:5) {
      MXlm <- lm(out[1:(nlen-k)] ~ inp[(k+1):nlen],adat)
      
      Mlag   = summary(MXlm)
      
      MXtab[i-1,(4*k)+1] <- Mlag$coefficients[2,1]
      MXtab[i-1,(4*k)+2] <- Mlag$coefficients[2,2]
      MXtab[i-1,(4*k)+3] <- Mlag$coefficients[2,4]
      MXtab[i-1,(4*k)+4] <- Mlag$r.squared
    }  
    
    ### NE
    adat <- rdat
    adat$out <- NE[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    for (k in 0:5) {
      NElm <- lm(out[1:(nlen-k)] ~ inp[(k+1):nlen],adat)
      
      Nlag   = summary(NElm)
      
      NEtab[i-1,(4*k)+1] <- Nlag$coefficients[2,1]
      NEtab[i-1,(4*k)+2] <- Nlag$coefficients[2,2]
      NEtab[i-1,(4*k)+3] <- Nlag$coefficients[2,4]
      NEtab[i-1,(4*k)+4] <- Nlag$r.squared
    }  
    
    ### SE
    adat <- rdat
    adat$out <- SE[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    for (k in 0:5) {
      SElm <- lm(out[1:(nlen-k)] ~ inp[(k+1):nlen],adat)
      
      Slag   = summary(SElm)
      
      SEtab[i-1,(4*k)+1] <- Slag$coefficients[2,1]
      SEtab[i-1,(4*k)+2] <- Slag$coefficients[2,2]
      SEtab[i-1,(4*k)+3] <- Slag$coefficients[2,4]
      SEtab[i-1,(4*k)+4] <- Slag$r.squared
    }
    
    # AKlm0 <- lm(out[1:(nlen-0)] ~ inp[1:nlen],adat)
    # AKlm1 <- lm(out[1:(nlen-1)] ~ inp[2:nlen],adat)
    # AKlm2 <- lm(out[1:(nlen-2)] ~ inp[3:nlen],adat)
    # AKlm3 <- lm(out[1:(nlen-3)] ~ inp[4:nlen],adat)
    # AKlm4 <- lm(out[1:(nlen-4)] ~ inp[5:nlen],adat)
    # AKlm5 <- lm(out[1:(nlen-5)] ~ inp[6:nlen],adat)
    # modelsummary(Amodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"AK_regress_",clim,"_",driver,"_div2SD.txt"))
    
    
  } # inputs & fish
  write.table(AKtab,paste0(datap,"GAK_regress_",clim,"_div2SD.csv"),sep=",",row.names=T)
  write.table(BStab,paste0(datap,"EBS_regress_",clim,"_div2SD.csv"),sep=",",row.names=T)
  write.table(CCtab,paste0(datap,"CCE_regress_",clim,"_div2SD.csv"),sep=",",row.names=T)
  write.table(CHtab,paste0(datap,"CHK_regress_",clim,"_div2SD.csv"),sep=",",row.names=T)
  write.table(HItab,paste0(datap,"HI_regress_",clim,"_div2SD.csv"),sep=",",row.names=T)
  write.table(MXtab,paste0(datap,"MX_regress_",clim,"_div2SD.csv"),sep=",",row.names=T)
  write.table(NEtab,paste0(datap,"NE_regress_",clim,"_div2SD.csv"),sep=",",row.names=T)
  write.table(SEtab,paste0(datap,"SE_regress_",clim,"_div2SD.csv"),sep=",",row.names=T)
  
}


