# Mult linear regression of forcing on fish biomass

rm(list=ls())

library(dplyr)

setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/cpue_analysis/")
figp <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/corrs/"

### ------------------------------------------------------------
# load data
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regress_cpue/"
##CPUE
#obsfish
Feo <- read.csv(paste0(datap,"LMEs_corr_cpue_chlyrs15_driver_obsfish2015_maxcorr_posfood_F.csv"),sep=",",header = T,stringsAsFactors = F)
Peo <- read.csv(paste0(datap,"LMEs_corr_cpue_chlyrs15_driver_obsfish2015_maxcorr_posfood_P.csv"),sep=",",header = T,stringsAsFactors = F)
Deo <- read.csv(paste0(datap,"LMEs_corr_cpue_chlyrs15_driver_obsfish2015_maxcorr_posfood_D.csv"),sep=",",header = T,stringsAsFactors = F)
Aeo <- read.csv(paste0(datap,"LMEs_corr_cpue_chlyrs15_driver_obsfish2015_maxcorr_posfood_A.csv"),sep=",",header = T,stringsAsFactors = F)
#ctrlfish
Fec <- read.csv(paste0(datap,"LMEs_corr_cpue_chlyrs15_inputs_feisty_maxcorr_posfood_F.csv"),sep=",",header = T,stringsAsFactors = F)
Pec <- read.csv(paste0(datap,"LMEs_corr_cpue_chlyrs15_inputs_feisty_maxcorr_posfood_P.csv"),sep=",",header = T,stringsAsFactors = F)
Dec <- read.csv(paste0(datap,"LMEs_corr_cpue_chlyrs15_inputs_feisty_maxcorr_posfood_D.csv"),sep=",",header = T,stringsAsFactors = F)
Aec <- read.csv(paste0(datap,"LMEs_corr_cpue_chlyrs15_inputs_feisty_maxcorr_posfood_A.csv"),sep=",",header = T,stringsAsFactors = F)
##Catch
#obsfish
Fco <- read.csv(paste0(datap,"LMEs_corr_catch_chlyrs15_inputs_obsfish2015_maxcorr_posfood_F.csv"),sep=",",header = T,stringsAsFactors = F)
Pco <- read.csv(paste0(datap,"LMEs_corr_catch_chlyrs15_inputs_obsfish2015_maxcorr_posfood_P.csv"),sep=",",header = T,stringsAsFactors = F)
Dco <- read.csv(paste0(datap,"LMEs_corr_catch_chlyrs15_inputs_obsfish2015_maxcorr_posfood_D.csv"),sep=",",header = T,stringsAsFactors = F)
Aco <- read.csv(paste0(datap,"LMEs_corr_catch_chlyrs15_inputs_obsfish2015_maxcorr_posfood_A.csv"),sep=",",header = T,stringsAsFactors = F)
#ctrlfish
Fcc <- read.csv(paste0(datap,"LMEs_corr_catch_chlyrs15_inputs_feisty_maxcorr_posfood_F.csv"),sep=",",header = T,stringsAsFactors = F)
Pcc <- read.csv(paste0(datap,"LMEs_corr_catch_chlyrs15_inputs_feisty_maxcorr_posfood_P.csv"),sep=",",header = T,stringsAsFactors = F)
Dcc <- read.csv(paste0(datap,"LMEs_corr_catch_chlyrs15_inputs_feisty_maxcorr_posfood_D.csv"),sep=",",header = T,stringsAsFactors = F)
Acc <- read.csv(paste0(datap,"LMEs_corr_catch_chlyrs15_inputs_feisty_maxcorr_posfood_A.csv"),sep=",",header = T,stringsAsFactors = F)


### Summarize R2 ---------------------------------------------------------
## CPUE
#ctrlfish
cAec <- summary((Aec$corr)^2)
tab_cpueC <- data.frame(4,4)
tab_cpueC[1,1] <- cAec[1]
tab_cpueC[1,2] <- cAec[4]
tab_cpueC[1,3] <- cAec[6]

cFec <- summary((Fec$corr)^2)
tab_cpueC[2,1] <- cFec[1]
tab_cpueC[2,2] <- cFec[4]
tab_cpueC[2,3] <- cFec[6]

cPec <- summary((Pec$corr)^2)
tab_cpueC[3,1] <- cPec[1]
tab_cpueC[3,2] <- cPec[4]
tab_cpueC[3,3] <- cPec[6]

cDec <- summary((Dec$corr)^2)
tab_cpueC[4,1] <- cDec[1]
tab_cpueC[4,2] <- cDec[4]
tab_cpueC[4,3] <- cDec[6]

tab_cpueC[,4] <- "Const"
tab_cpueC[,5] <- "CPUE"
names(tab_cpueC) <- c("Min","Mean","Max","Effort","Cvar")
row.names(tab_cpueC) <- c("All","Forage","LgPel","Dem")

#obsfish
cAeo <- summary((Aeo$corr)^2)
tab_cpueO <- data.frame(4,4)
tab_cpueO[1,1] <- cAeo[1]
tab_cpueO[1,2] <- cAeo[4]
tab_cpueO[1,3] <- cAeo[6]

cFeo <- summary((Feo$corr)^2)
tab_cpueO[2,1] <- cFeo[1]
tab_cpueO[2,2] <- cFeo[4]
tab_cpueO[2,3] <- cFeo[6]

cPeo <- summary((Peo$corr)^2)
tab_cpueO[3,1] <- cPeo[1]
tab_cpueO[3,2] <- cPeo[4]
tab_cpueO[3,3] <- cPeo[6]

cDeo <- summary((Deo$corr)^2)
tab_cpueO[4,1] <- cDeo[1]
tab_cpueO[4,2] <- cDeo[4]
tab_cpueO[4,3] <- cDeo[6]

tab_cpueO[,4] <- "Obs"
tab_cpueO[,5] <- "CPUE"
names(tab_cpueO) <- c("Min","Mean","Max","Effort","Cvar")
row.names(tab_cpueO) <- c("All","Forage","LgPel","Dem")

## Catch
#ctrlfish
cAcc <- summary((Acc$corr)^2)
tab_catchC <- data.frame(4,4)
tab_catchC[1,1] <- cAcc[1]
tab_catchC[1,2] <- cAcc[4]
tab_catchC[1,3] <- cAcc[6]

cFcc <- summary((Fcc$corr)^2)
tab_catchC[2,1] <- cFcc[1]
tab_catchC[2,2] <- cFcc[4]
tab_catchC[2,3] <- cFcc[6]

cPcc <- summary((Pcc$corr)^2)
tab_catchC[3,1] <- cPcc[1]
tab_catchC[3,2] <- cPcc[4]
tab_catchC[3,3] <- cPcc[6]

cDcc <- summary((Dcc$corr)^2)
tab_catchC[4,1] <- cDcc[1]
tab_catchC[4,2] <- cDcc[4]
tab_catchC[4,3] <- cDcc[6]

tab_catchC[,4] <- "Const"
tab_catchC[,5] <- "Catch"
names(tab_catchC) <- c("Min","Mean","Max","Effort","Cvar")
row.names(tab_catchC) <- c("All","Forage","LgPel","Dem")

#obsfish
cAco <- summary((Aco$corr)^2)
tab_catchO <- data.frame(4,4)
tab_catchO[1,1] <- cAco[1]
tab_catchO[1,2] <- cAco[4]
tab_catchO[1,3] <- cAco[6]

cFco <- summary((Fco$corr)^2)
tab_catchO[2,1] <- cFco[1]
tab_catchO[2,2] <- cFco[4]
tab_catchO[2,3] <- cFco[6]

cPco <- summary((Pco$corr)^2)
tab_catchO[3,1] <- cPco[1]
tab_catchO[3,2] <- cPco[4]
tab_catchO[3,3] <- cPco[6]

cDco <- summary((Dco$corr)^2)
tab_catchO[4,1] <- cDco[1]
tab_catchO[4,2] <- cDco[4]
tab_catchO[4,3] <- cDco[6]

tab_catchO[,4] <- "Obs"
tab_catchO[,5] <- "Catch"
names(tab_catchO) <- c("Min","Mean","Max","Effort","Cvar")
row.names(tab_catchO) <- c("All","Forage","LgPel","Dem")


# write.table(tab_cpueC,paste0(datap,"Range_LMEs_R2_cpue_chlyrs15.csv"),sep=",",row.names=F)
# write.table(tab_cpueO,paste0(datap,"Range_LMEs_R2_cpue_chlyrs15_obsfish2015.csv"),sep=",",row.names=F)
# write.table(tab_catchC,paste0(datap,"Range_LMEs_R2_catch_chlyrs15.csv"),sep=",",row.names=F)
# write.table(tab_catchO,paste0(datap,"Range_LMEs_R2_catch_chlyrs15_obsfish2015.csv"),sep=",",row.names=F)

##Stack into one
tab_all <- bind_rows(tab_cpueC,tab_cpueO,tab_catchC,tab_catchO)

write.table(tab_all,paste0(datap,"Range_mean_LMEs_R2_chlyrs15_fntypes.csv"),sep=",",row.names=F)


  