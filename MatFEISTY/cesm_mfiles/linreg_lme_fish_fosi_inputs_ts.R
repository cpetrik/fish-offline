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
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/"
CH <- read.csv(paste0(datap,"CHK_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
BS <- read.csv(paste0(datap,"EBS_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
AK <- read.csv(paste0(datap,"GAK_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
HI <- read.csv(paste0(datap,"HI_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
CC <- read.csv(paste0(datap,"CCE_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
MX <- read.csv(paste0(datap,"GMX_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
SE <- read.csv(paste0(datap,"SE_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)
NE <- read.csv(paste0(datap,"NE_fosi_inputs_fish_annual_means_anomalies.csv"),sep=",",header = T,stringsAsFactors = F)


### LM of climate ---------------------------------------------------------
#Loop over forcing and fish types
for (i in 2:6) {
  for (j in 7:14) {
    
    ### AK
    inp <- AK[,i]/(2*sd(AK[,i]))
    out <- AK[,j]
    
    AKlm0 <- lm(out ~ inp)
    AKlm1 <- lm(out[1:(68-1)] ~ inp[2:68])
    AKlm2 <- lm(out[1:(68-2)] ~ inp[3:68])
    AKlm3 <- lm(out[1:(68-3)] ~ inp[4:68])
    AKlm4 <- lm(out[1:(68-4)] ~ inp[5:68])
    AKlm5 <- lm(out[1:(68-5)] ~ inp[6:68])
    
    Amodels <- list(
      "lag0"    = AKlm0,
      "lag1"    = AKlm1,
      "lag2"    = AKlm2,
      "lag3"    = AKlm3,
      "lag4"    = AKlm4,
      "lag5"    = AKlm5
    )
    driver <- names(AK)[i]
    fish <- names(AK)[j]
    modelsummary(Amodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"AK_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Amodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"AK_regress_",fish,"_",driver,"_div2SD.txt"))
    
    ### BS
    inp <- BS[,i]/(2*sd(BS[,i]))
    out <- BS[,j]
    
    BSlm0 <- lm(out ~ inp)
    BSlm1 <- lm(out[1:(68-1)] ~ inp[2:68])
    BSlm2 <- lm(out[1:(68-2)] ~ inp[3:68])
    BSlm3 <- lm(out[1:(68-3)] ~ inp[4:68])
    BSlm4 <- lm(out[1:(68-4)] ~ inp[5:68])
    BSlm5 <- lm(out[1:(68-5)] ~ inp[6:68])
    
    Bmodels <- list(
      "lag0"    = BSlm0,
      "lag1"    = BSlm1,
      "lag2"    = BSlm2,
      "lag3"    = BSlm3,
      "lag4"    = BSlm4,
      "lag5"    = BSlm5
    )
    driver <- names(BS)[i]
    fish <- names(BS)[j]
    modelsummary(Bmodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"BS_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Bmodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"BS_regress_",fish,"_",driver,"_div2SD.txt"))
    
    
    ### CC
    inp <- CC[,i]/(2*sd(CC[,i]))
    out <- CC[,j]
    
    CClm0 <- lm(out ~ inp)
    CClm1 <- lm(out[1:(68-1)] ~ inp[2:68])
    CClm2 <- lm(out[1:(68-2)] ~ inp[3:68])
    CClm3 <- lm(out[1:(68-3)] ~ inp[4:68])
    CClm4 <- lm(out[1:(68-4)] ~ inp[5:68])
    CClm5 <- lm(out[1:(68-5)] ~ inp[6:68])
    
    Emodels <- list(
      "lag0"    = CClm0,
      "lag1"    = CClm1,
      "lag2"    = CClm2,
      "lag3"    = CClm3,
      "lag4"    = CClm4,
      "lag5"    = CClm5
    )
    driver <- names(CC)[i]
    fish <- names(CC)[j]
    modelsummary(Emodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"CC_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Emodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"CC_regress_",fish,"_",driver,"_div2SD.txt"))
    
    
    ### CH
    inp <- CH[,i]/(2*sd(CH[,i]))
    out <- CH[,j]
    
    CHlm0 <- lm(out ~ inp)
    CHlm1 <- lm(out[1:(68-1)] ~ inp[2:68])
    CHlm2 <- lm(out[1:(68-2)] ~ inp[3:68])
    CHlm3 <- lm(out[1:(68-3)] ~ inp[4:68])
    CHlm4 <- lm(out[1:(68-4)] ~ inp[5:68])
    CHlm5 <- lm(out[1:(68-5)] ~ inp[6:68])
    
    Cmodels <- list(
      "lag0"    = CHlm0,
      "lag1"    = CHlm1,
      "lag2"    = CHlm2,
      "lag3"    = CHlm3,
      "lag4"    = CHlm4,
      "lag5"    = CHlm5
    )
    driver <- names(CH)[i]
    fish <- names(CH)[j]
    modelsummary(Cmodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"CH_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Cmodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"CH_regress_",fish,"_",driver,"_div2SD.txt"))
    
    
    ### HI
    inp <- HI[,i]/(2*sd(HI[,i]))
    out <- HI[,j]
    
    HIlm0 <- lm(out ~ inp)
    HIlm1 <- lm(out[1:(68-1)] ~ inp[2:68])
    HIlm2 <- lm(out[1:(68-2)] ~ inp[3:68])
    HIlm3 <- lm(out[1:(68-3)] ~ inp[4:68])
    HIlm4 <- lm(out[1:(68-4)] ~ inp[5:68])
    HIlm5 <- lm(out[1:(68-5)] ~ inp[6:68])
    
    Hmodels <- list(
      "lag0"    = HIlm0,
      "lag1"    = HIlm1,
      "lag2"    = HIlm2,
      "lag3"    = HIlm3,
      "lag4"    = HIlm4,
      "lag5"    = HIlm5
    )
    driver <- names(HI)[i]
    fish <- names(HI)[j]
    modelsummary(Hmodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"HI_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Hmodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"HI_regress_",fish,"_",driver,"_div2SD.txt"))
    
    
    ### MX
    inp <- MX[,i]/(2*sd(MX[,i]))
    out <- MX[,j]
    
    MXlm0 <- lm(out ~ inp)
    MXlm1 <- lm(out[1:(68-1)] ~ inp[2:68])
    MXlm2 <- lm(out[1:(68-2)] ~ inp[3:68])
    MXlm3 <- lm(out[1:(68-3)] ~ inp[4:68])
    MXlm4 <- lm(out[1:(68-4)] ~ inp[5:68])
    MXlm5 <- lm(out[1:(68-5)] ~ inp[6:68])
    
    Mmodels <- list(
      "lag0"    = MXlm0,
      "lag1"    = MXlm1,
      "lag2"    = MXlm2,
      "lag3"    = MXlm3,
      "lag4"    = MXlm4,
      "lag5"    = MXlm5
    )
    driver <- names(MX)[i]
    fish <- names(MX)[j]
    modelsummary(Mmodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"MX_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Mmodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"MX_regress_",fish,"_",driver,"_div2SD.txt"))
    
    
    ### NE
    inp <- NE[,i]/(2*sd(NE[,i]))
    out <- NE[,j]
    
    NElm0 <- lm(out ~ inp)
    NElm1 <- lm(out[1:(68-1)] ~ inp[2:68])
    NElm2 <- lm(out[1:(68-2)] ~ inp[3:68])
    NElm3 <- lm(out[1:(68-3)] ~ inp[4:68])
    NElm4 <- lm(out[1:(68-4)] ~ inp[5:68])
    NElm5 <- lm(out[1:(68-5)] ~ inp[6:68])
    
    Nmodels <- list(
      "lag0"    = NElm0,
      "lag1"    = NElm1,
      "lag2"    = NElm2,
      "lag3"    = NElm3,
      "lag4"    = NElm4,
      "lag5"    = NElm5
    )
    driver <- names(NE)[i]
    fish <- names(NE)[j]
    modelsummary(Nmodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"NE_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Nmodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"NE_regress_",fish,"_",driver,"_div2SD.txt"))
    
    ### SE
    inp <- SE[,i]/(2*sd(SE[,i]))
    out <- SE[,j]
    
    SElm0 <- lm(out ~ inp)
    SElm1 <- lm(out[1:(68-1)] ~ inp[2:68])
    SElm2 <- lm(out[1:(68-2)] ~ inp[3:68])
    SElm3 <- lm(out[1:(68-3)] ~ inp[4:68])
    SElm4 <- lm(out[1:(68-4)] ~ inp[5:68])
    SElm5 <- lm(out[1:(68-5)] ~ inp[6:68])
    
    Smodels <- list(
      "lag0"    = SElm0,
      "lag1"    = SElm1,
      "lag2"    = SElm2,
      "lag3"    = SElm3,
      "lag4"    = SElm4,
      "lag5"    = SElm5
    )
    driver <- names(SE)[i]
    fish <- names(SE)[j]
    modelsummary(Smodels, fmt = "%.3e", estimate = "stars", output = paste0(datap,"SE_regress_",fish,"_",driver,"_div2SD.txt"))
    modelsummary(Smodels, fmt = "%.3e", estimate = "stars", output = paste0(datam,"SE_regress_",fish,"_",driver,"_div2SD.txt"))
    
  }
}

modelsummary(Smodels, fmt = "%.3e", statistic = 'p.value')

