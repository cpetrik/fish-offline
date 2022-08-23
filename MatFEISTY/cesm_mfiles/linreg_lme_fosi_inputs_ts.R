# Find patterns in climate correlations

rm(list=ls())

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/"

### ------------------------------------------------------------
# load data
apath <- '/Users/cpetrik/Dropbox/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
C <- read.csv(paste0(apath, "Climate_anomalies_annual_means.csv"),sep=",",header = T,stringsAsFactors = F)

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
#CAN
Clm0 <- lm(Lzmeso~Lchl, data=CM)
Clm1 <- lm(Lzmeso~Lchl, data=subset(CM, biome==1))
Clm2 <- lm(Lzmeso~Lchl, data=subset(CM, biome==2))
Clm3 <- lm(Lzmeso~Lchl, data=subset(CM, biome==3))


### Save all summary stats
## Diff table for each biome
#GLobal
Gmodels <- list(
  "CAN"     = Clm0,
  "CMCC"    = Mlm0,
  "CNRM"    = Nlm0,
  "GFDL"    = Glm0,
  "IPSL"    = Ilm0,
  "UK"      = Ulm0,
  "obsGLMM" = Olm0,
  "obsSM"   = Slm0,
  "obsMO-S" = CSlm0,
  "obsMO-M" = CMlm0
)
modelsummary(Gmodels, fmt=2) #estimate = "stars")
modelsummary(Gmodels, fmt=2, output = paste0(datap,"regress_global_table_mgC.docx"))

#LC
Lmodels <- list(
  "CAN"     = Clm1,
  "CMCC"    = Mlm1,
  "CNRM"    = Nlm1,
  "GFDL"    = Glm1,
  "IPSL"    = Ilm1,
  "UK"      = Ulm1,
  "obsGLMM" = Olm1,
  "obsSM"   = Slm1,
  "obsMO-S" = CSlm1,
  "obsMO-M" = CMlm1
)
modelsummary(Lmodels, fmt=2, output = paste0(datap,"regress_LC_table_mgC.docx"))

#HCSS
Smodels <- list(
  "CAN"     = Clm2,
  "CMCC"    = Mlm2,
  "CNRM"    = Nlm2,
  "GFDL"    = Glm2,
  "IPSL"    = Ilm2,
  "UK"      = Ulm2,
  "obsGLMM" = Olm2,
  "obsSM"   = Slm2,
  "obsMO-S" = CSlm2,
  "obsMO-M" = CMlm2
)
modelsummary(Smodels, fmt=2, output = paste0(datap,"regress_HCSS_table_mgC.docx"))

#HCPS
Pmodels <- list(
  "CAN"     = Clm3,
  "CMCC"    = Mlm3,
  "CNRM"    = Nlm3,
  "GFDL"    = Glm3,
  "IPSL"    = Ilm3,
  "UK"      = Ulm3,
  "obsGLMM" = Olm3,
  "obsSM"   = Slm3,
  "obsMO-S" = CSlm3,
  "obsMO-M" = CMlm3
)
modelsummary(Pmodels, fmt=2, output = paste0(datap,"regress_HCPS_table_mgC.docx"))

