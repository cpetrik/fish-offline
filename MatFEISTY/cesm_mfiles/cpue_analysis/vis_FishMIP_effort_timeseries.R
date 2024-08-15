################################################################################

# FishMIP Phase3a Effort

################################################################################

rm(list=ls())

dir <- '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

### Load data
efrt <- read.csv(paste0(dir,'effort_isimip3a_histsoc_1841_2010.csv'),sep=",",header = T,stringsAsFactors = T)

efrt <- subset(efrt, Year>1947)

cid <- c("Year","LME","FGroup","NomActive")
efrt <- efrt[,cid]

levels( efrt$FGroup )

#odir <- paste0('/Volumes/petrik-lab/Feisty/NC/NEMURO/',cfile,'/')
#write.table(all,paste0(odir,"FEISTY_NEMURO_annual_outputs_popprod.csv"),sep=",",row.names=F)



