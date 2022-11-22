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
#Loop over forcing and climate indices
for (i in 2:6) {
  for (j in 1:11) {
    fish <- names(C)[j]
    rdat <- as.data.frame(C[,j])
    names(rdat) <- "inp"
    rdat$inp <- rdat$inp/(2*sd(rdat$inp, na.rm = T))
    
    ### AK
    adat <- rdat
    adat$out <- AK[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    AKlm0 <- lm(out ~ inp,adat)
    AKlm1 <- lm(out[1:(nlen-1)] ~ inp[2:nlen],adat)
    AKlm2 <- lm(out[1:(nlen-2)] ~ inp[3:nlen],adat)
    AKlm3 <- lm(out[1:(nlen-3)] ~ inp[4:nlen],adat)
    AKlm4 <- lm(out[1:(nlen-4)] ~ inp[5:nlen],adat)
    AKlm5 <- lm(out[1:(nlen-5)] ~ inp[6:nlen],adat)
    
    Amodels <- list(
      "lag0"    = AKlm0,
      "lag1"    = AKlm1,
      "lag2"    = AKlm2,
      "lag3"    = AKlm3,
      "lag4"    = AKlm4,
      "lag5"    = AKlm5
    )
    driver <- names(AK)[i]
    modelsummary(Amodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"AK_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Amodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"AK_regress_",fish,"_",driver,"_div2SD.txt"))
    
    ### BS
    adat <- rdat
    adat$out <- BS[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    BSlm0 <- lm(out ~ inp,adat)
    BSlm1 <- lm(out[1:(nlen-1)] ~ inp[2:nlen],adat)
    BSlm2 <- lm(out[1:(nlen-2)] ~ inp[3:nlen],adat)
    BSlm3 <- lm(out[1:(nlen-3)] ~ inp[4:nlen],adat)
    BSlm4 <- lm(out[1:(nlen-4)] ~ inp[5:nlen],adat)
    BSlm5 <- lm(out[1:(nlen-5)] ~ inp[6:nlen],adat)
    
    Bmodels <- list(
      "lag0"    = BSlm0,
      "lag1"    = BSlm1,
      "lag2"    = BSlm2,
      "lag3"    = BSlm3,
      "lag4"    = BSlm4,
      "lag5"    = BSlm5
    )
    driver <- names(BS)[i]
    modelsummary(Bmodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"BS_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Bmodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"BS_regress_",fish,"_",driver,"_div2SD.txt"))
    
    
    ### CC
    adat <- rdat
    adat$out <- CC[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    CClm0 <- lm(out ~ inp,adat)
    CClm1 <- lm(out[1:(nlen-1)] ~ inp[2:nlen],adat)
    CClm2 <- lm(out[1:(nlen-2)] ~ inp[3:nlen],adat)
    CClm3 <- lm(out[1:(nlen-3)] ~ inp[4:nlen],adat)
    CClm4 <- lm(out[1:(nlen-4)] ~ inp[5:nlen],adat)
    CClm5 <- lm(out[1:(nlen-5)] ~ inp[6:nlen],adat)
    
    Emodels <- list(
      "lag0"    = CClm0,
      "lag1"    = CClm1,
      "lag2"    = CClm2,
      "lag3"    = CClm3,
      "lag4"    = CClm4,
      "lag5"    = CClm5
    )
    driver <- names(CC)[i]
    modelsummary(Emodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"CC_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Emodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"CC_regress_",fish,"_",driver,"_div2SD.txt"))
    
    
    ### CH
    adat <- rdat
    adat$out <- CH[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    CHlm0 <- lm(out ~ inp,adat)
    CHlm1 <- lm(out[1:(nlen-1)] ~ inp[2:nlen],adat)
    CHlm2 <- lm(out[1:(nlen-2)] ~ inp[3:nlen],adat)
    CHlm3 <- lm(out[1:(nlen-3)] ~ inp[4:nlen],adat)
    CHlm4 <- lm(out[1:(nlen-4)] ~ inp[5:nlen],adat)
    CHlm5 <- lm(out[1:(nlen-5)] ~ inp[6:nlen],adat)
    
    Cmodels <- list(
      "lag0"    = CHlm0,
      "lag1"    = CHlm1,
      "lag2"    = CHlm2,
      "lag3"    = CHlm3,
      "lag4"    = CHlm4,
      "lag5"    = CHlm5
    )
    driver <- names(CH)[i]
    modelsummary(Cmodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"CH_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Cmodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"CH_regress_",fish,"_",driver,"_div2SD.txt"))
    
    
    ### HI
    adat <- rdat
    adat$out <- HI[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    HIlm0 <- lm(out ~ inp,adat)
    HIlm1 <- lm(out[1:(nlen-1)] ~ inp[2:nlen],adat)
    HIlm2 <- lm(out[1:(nlen-2)] ~ inp[3:nlen],adat)
    HIlm3 <- lm(out[1:(nlen-3)] ~ inp[4:nlen],adat)
    HIlm4 <- lm(out[1:(nlen-4)] ~ inp[5:nlen],adat)
    HIlm5 <- lm(out[1:(nlen-5)] ~ inp[6:nlen],adat)
    
    Hmodels <- list(
      "lag0"    = HIlm0,
      "lag1"    = HIlm1,
      "lag2"    = HIlm2,
      "lag3"    = HIlm3,
      "lag4"    = HIlm4,
      "lag5"    = HIlm5
    )
    driver <- names(HI)[i]
    modelsummary(Hmodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"HI_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Hmodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"HI_regress_",fish,"_",driver,"_div2SD.txt"))
    
    
    ### MX
    adat <- rdat
    adat$out <- MX[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    MXlm0 <- lm(out ~ inp,adat)
    MXlm1 <- lm(out[1:(nlen-1)] ~ inp[2:nlen],adat)
    MXlm2 <- lm(out[1:(nlen-2)] ~ inp[3:nlen],adat)
    MXlm3 <- lm(out[1:(nlen-3)] ~ inp[4:nlen],adat)
    MXlm4 <- lm(out[1:(nlen-4)] ~ inp[5:nlen],adat)
    MXlm5 <- lm(out[1:(nlen-5)] ~ inp[6:nlen],adat)
    
    Mmodels <- list(
      "lag0"    = MXlm0,
      "lag1"    = MXlm1,
      "lag2"    = MXlm2,
      "lag3"    = MXlm3,
      "lag4"    = MXlm4,
      "lag5"    = MXlm5
    )
    driver <- names(MX)[i]
    modelsummary(Mmodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"MX_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Mmodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"MX_regress_",fish,"_",driver,"_div2SD.txt"))
    
    
    ### NE
    adat <- rdat
    adat$out <- NE[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    NElm0 <- lm(out ~ inp,adat)
    NElm1 <- lm(out[1:(nlen-1)] ~ inp[2:nlen],adat)
    NElm2 <- lm(out[1:(nlen-2)] ~ inp[3:nlen],adat)
    NElm3 <- lm(out[1:(nlen-3)] ~ inp[4:nlen],adat)
    NElm4 <- lm(out[1:(nlen-4)] ~ inp[5:nlen],adat)
    NElm5 <- lm(out[1:(nlen-5)] ~ inp[6:nlen],adat)
    
    Nmodels <- list(
      "lag0"    = NElm0,
      "lag1"    = NElm1,
      "lag2"    = NElm2,
      "lag3"    = NElm3,
      "lag4"    = NElm4,
      "lag5"    = NElm5
    )
    driver <- names(NE)[i]
    modelsummary(Nmodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"NE_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Nmodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"NE_regress_",fish,"_",driver,"_div2SD.txt"))
    
    ### SE
    adat <- rdat
    adat$out <- SE[,i]
    adat <- na.omit(adat)
    nlen <- dim(adat)[1]
    
    SElm0 <- lm(out ~ inp,adat)
    SElm1 <- lm(out[1:(nlen-1)] ~ inp[2:nlen],adat)
    SElm2 <- lm(out[1:(nlen-2)] ~ inp[3:nlen],adat)
    SElm3 <- lm(out[1:(nlen-3)] ~ inp[4:nlen],adat)
    SElm4 <- lm(out[1:(nlen-4)] ~ inp[5:nlen],adat)
    SElm5 <- lm(out[1:(nlen-5)] ~ inp[6:nlen],adat)
    
    Smodels <- list(
      "lag0"    = SElm0,
      "lag1"    = SElm1,
      "lag2"    = SElm2,
      "lag3"    = SElm3,
      "lag4"    = SElm4,
      "lag5"    = SElm5
    )
    driver <- names(SE)[i]
    modelsummary(Smodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"SE_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Smodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"SE_regress_",fish,"_",driver,"_div2SD.txt"))
    
  }
}

modelsummary(Smodels, fmt = "%.3e", statistic = 'p.value')

