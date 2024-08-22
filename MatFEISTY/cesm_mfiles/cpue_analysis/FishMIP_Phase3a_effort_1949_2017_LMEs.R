################################################################################

# Aggregate F, P, and D again
# Only 1948-2015 saved for comp to FOSI

################################################################################

rm(list=ls())

library( plyr )

setwd("/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
spath <- "/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/"
fpath <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Fish-MIP/Phase3/fishing/"

## Load Yannick's data
#saup <- read.csv(paste0(fpath,"catch_histsoc_1869_2017_EEZ_addFAO.csv"),sep=",",header = T,stringsAsFactors = T)
# 1841-2010 effort only has LMEs 0-32
effort0 <- read.csv(paste0(spath,"effort_isimip3a_histsoc_1841_2010.csv"),sep=",",header = T,stringsAsFactors = T)
effort1 <- read.csv(paste0(spath,"effort_histsoc_1841_2017_EEZ_addFAO.csv"),sep=",",header = T,stringsAsFactors = T)
#effort <- read.csv(paste0(fpath,"effort_histsoc_1961_2010.csv"),sep=",",header = T,stringsAsFactors = T)

#Check new yrs
summary(effort0) #50585449
# Year             Sector            fao_area          LME           eez_country_name         SAUP      
# Min.   :1841   Artisanal :10144734   Min.   :18.00   Min.   : 0.00   High Seas  : 3339732   Min.   :  8.0  
# 1st Qu.:1884   Industrial:40440715   1st Qu.:27.00   1st Qu.:11.00   Russian Fed: 2324603   1st Qu.:250.0  
# Median :1927                         Median :37.00   Median :26.00   Indonesia  : 1733079   Median :410.0  
# Mean   :1928                         Mean   :44.76   Mean   :24.83   USA        : 1677043   Mean   :460.2  
# 3rd Qu.:1975                         3rd Qu.:61.00   3rd Qu.:36.00   Canada     : 1619055   3rd Qu.:643.0  
# Max.   :2010                         Max.   :88.00   Max.   :66.00   US (Alaska): 1600346   Max.   :892.0  
#                                                                      (Other)    :38291591                  

# Gear                              FGroup                         NomActive           Phase         
# Trawl_Midwater_or_Unsp: 6000391   pelagic30-90cm     : 5096287   Min.   :        0   experiment:17015493  
# Gillnets              : 5985173   demersal<30cm      : 5045565   1st Qu.:        0   transition:33569956  
# Lines_Longlines       : 5834719   pelagic>=90cm      : 4973619   Median :        3                        
# Others_Others         : 5357098   demersal30-90cm    : 4215150   Mean   :    16717                        
# Trawl_Bottom          : 5138961   benthopelagic>=90cm: 3368950   3rd Qu.:      142                        
# Seine_Purse_Seine     : 4386060   pelagic<30cm       : 2978346   Max.   :233495502                        
# (Other)               :17883047   (Other)            :24907532                                  

summary(effort1) #53324778
# Year         eez_country_name         SAUP                           Gear                          FGroup        
# Min.   :1841   High Seas  : 3645502   Min.   :  8.0   Trawl_Midwater_or_Unsp: 6291467   pelagic30-90cm     : 5332127  
# 1st Qu.:1886   Russian Fed: 2411741   1st Qu.:250.0   Gillnets              : 6290780   pelagic>=90cm      : 5311188  
# Median :1932   Indonesia  : 1793676   Median :410.0   Lines_Longlines       : 6139193   demersal<30cm      : 5307582  
# Mean   :1932   USA        : 1763229   Mean   :461.5   Others_Others         : 5615941   demersal30-90cm    : 4387668  
# 3rd Qu.:1980   Canada     : 1661037   3rd Qu.:643.0   Trawl_Bottom          : 5374664   benthopelagic>=90cm: 3513976  
# Max.   :2017   US (Alaska): 1656784   Max.   :892.0   Seine_Purse_Seine     : 4684404   pelagic<30cm       : 3116514  
#                (Other)    :40392809                   (Other)               :18928329   (Other)            :26355723  

# Sector                LME             fao_area        NomActive           Phase         
# Artisanal :10750206   Min.   : 0.00   Min.   :18.00   Min.   :        0   experiment:17015493  
# Industrial:42574572   1st Qu.:11.00   1st Qu.:27.00   1st Qu.:        0   spin-up   : 5829800  
#                       Median :26.00   Median :37.00   Median :        3   transition:27740156  
#                       Mean   :24.64   Mean   :44.84   Mean   :    18698   validation: 2739329  
#                       3rd Qu.:36.00   3rd Qu.:61.00   3rd Qu.:      157                        
#                       Max.   :66.00   Max.   :88.00   Max.   :233495502                        


## Define subsets in FGroup
levels(effort1$FGroup)

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
feisty0 <- subset.data.frame(effort0, Year >= 1948)
feisty0 <- subset.data.frame(feisty0, Year <= 2015)
feisty0 <- subset.data.frame(feisty0, FGroup %in% fntypes)
feisty0$FGroup <- droplevels(feisty0$FGroup)
summary(feisty0)

feisty1 <- subset.data.frame(effort1, Year >= 1948)
feisty1 <- subset.data.frame(feisty1, FGroup %in% fntypes)
feisty1$FGroup <- droplevels(feisty1$FGroup)
summary(feisty1)

## Assign to a FEISTY fn type
feisty0$FF <- 0
feisty0$LP <- 0
feisty0$DF <- 0

feisty1$FF <- 0
feisty1$LP <- 0
feisty1$DF <- 0

feisty0$FF[feisty0$FGroup %in% ffgn] <- 1
feisty0$LP[feisty0$FGroup %in% pfgn] <- 1
feisty0$DF[feisty0$FGroup %in% dfgn] <- 1
feisty0$LP[feisty0$FGroup %in% pdgn] <- 0.5
feisty0$DF[feisty0$FGroup %in% pdgn] <- 0.5

feisty1$FF[feisty1$FGroup %in% ffgn] <- 1
feisty1$LP[feisty1$FGroup %in% pfgn] <- 1
feisty1$DF[feisty1$FGroup %in% dfgn] <- 1
feisty1$LP[feisty1$FGroup %in% pdgn] <- 0.5
feisty1$DF[feisty1$FGroup %in% pdgn] <- 0.5

## Calc catch with weights
feisty0$Fweffort <- feisty0$NomActive * feisty0$FF
feisty0$Pweffort <- feisty0$NomActive * feisty0$LP
feisty0$Dweffort <- feisty0$NomActive * feisty0$DF

feisty1$Fweffort <- feisty1$NomActive * feisty1$FF
feisty1$Pweffort <- feisty1$NomActive * feisty1$LP
feisty1$Dweffort <- feisty1$NomActive * feisty1$DF

## Calc totals by Year & LME
lme0 <- ddply( feisty0, .( Year, LME ), summarize, 
              Fweffort = sum( as.numeric(Fweffort), na.rm = TRUE ),
              Pweffort = sum( as.numeric(Pweffort), na.rm = TRUE ),
              Dweffort = sum( as.numeric(Dweffort), na.rm = TRUE ))

lme1 <- ddply( feisty1, .( Year, LME ), summarize, 
              Fweffort = sum( as.numeric(Fweffort), na.rm = TRUE ),
              Pweffort = sum( as.numeric(Pweffort), na.rm = TRUE ),
              Dweffort = sum( as.numeric(Dweffort), na.rm = TRUE ))

## Save in both places
write.table(lme1,paste0(fpath,"FishMIP_Phase3a_LME_Effort_annual_1948-2017.csv"),sep=",",row.names=F)
write.table(lme1,paste0(spath,"FishMIP_Phase3a_LME_Effort_annual_1948-2017.csv"),sep=",",row.names=F)



