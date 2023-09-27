################################################################################

# Aggregate F, P, and D again
# Only 1948-2015 saved for comp to FOSI

################################################################################

rm(list=ls())

library( plyr )

setwd("/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
fpath <- "/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/"

## Load Yannick's data
saup <- read.csv(paste0(fpath,"catch_histsoc_1869_2017_EEZ_addFAO.csv"),sep=",",header = T,stringsAsFactors = T)


## Define subsets in FGroup
levels(saup$FGroup)

fntypes <- c("bathydemersal>=90cm","bathydemersal30-90cm",
             "bathypelagic<30cm","bathypelagic>=90cm","bathypelagic30-90cm",   
             "benthopelagic<30cm","benthopelagic>=90cm","benthopelagic30-90cm",  
             "demersal>=90cm","demersal30-90cm",      
             "cephalopods","flatfish<90cm","flatfish>=90cm",          
             "pelagic<30cm","pelagic>=90cm","pelagic30-90cm",       
             "rays>=90cm","shark>=90cm",
             "reef-associated>=90cm","reef-associated30-90cm")

ffgn <- c("bathypelagic<30cm","cephalopods","pelagic<30cm")
pfgn <- c("bathypelagic>=90cm","bathypelagic30-90cm","pelagic>=90cm","pelagic30-90cm")
dfgn <- c("bathydemersal>=90cm","bathydemersal30-90cm","benthopelagic>=90cm","benthopelagic30-90cm",
          "demersal>=90cm","demersal30-90cm","flatfish>=90cm","reef-associated>=90cm","reef-associated30-90cm")
pdgn <- c("shark>=90cm","rays>=90cm")

## Subset
feisty1 <- subset.data.frame(saup, Year >= 1948)
feisty1 <- subset.data.frame(feisty1, Year <= 2015)
feisty1 <- subset.data.frame(feisty1, FGroup %in% fntypes)
feisty1$FGroup <- droplevels(feisty1$FGroup)
summary(feisty1)

## Assign to a FEISTY fn type
feisty1$FF <- 0
feisty1$LP <- 0
feisty1$DF <- 0

feisty1$FF[feisty1$FGroup %in% ffgn] <- 1
feisty1$LP[feisty1$FGroup %in% pfgn] <- 1
feisty1$DF[feisty1$FGroup %in% dfgn] <- 1
feisty1$LP[feisty1$FGroup %in% pdgn] <- 0.5
feisty1$DF[feisty1$FGroup %in% pdgn] <- 0.5

## Calc catch with weights
feisty1$Fwcatch <- feisty1$Reported * feisty1$FF
feisty1$Pwcatch <- feisty1$Reported * feisty1$LP
feisty1$Dwcatch <- feisty1$Reported * feisty1$DF

## Calc totals by Year & LME
lme <- ddply( feisty1, .( Year, LME ), summarize, 
              FwCatch = sum( as.numeric(Fwcatch), na.rm = TRUE ),
              PwCatch = sum( as.numeric(Pwcatch), na.rm = TRUE ),
              DwCatch = sum( as.numeric(Dwcatch), na.rm = TRUE ))

write.table(lme,paste0(fpath,"FishMIP_Phase3a_LME_Catch_annual_1948-2015.csv"),sep=",",row.names=F)



