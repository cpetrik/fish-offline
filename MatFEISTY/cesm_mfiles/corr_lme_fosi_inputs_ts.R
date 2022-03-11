# Find patterns in climate correlations

rm(list=ls())

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/"

library(Hmisc)
library(ggplot2)
library(gridExtra)
library(corrgram)
library(scatterplot3d)
library(rpart)
library(tree)
library(pvclust)
library(mclust)
library(dendextend)
library(gplots)
library("corrplot")
library(scico)
library(reshape2)

### -------------------------------- Raw 4rt -------------------------------
# load data
ddir <- "/Volumes/MIP/GCM_DATA/CESM/FOSI/"
sigM <- read.csv(paste0(ddir,"LME_fosi_inputs_sigLzooC.csv"),sep=",",header = F,stringsAsFactors = T)
sigL <- read.csv(paste0(ddir,"LME_fosi_inputs_sigZooLoss.csv"),sep=",",header = F,stringsAsFactors = T)
sigP <- read.csv(paste0(ddir,"LME_fosi_inputs_sigTP.csv"),sep=",",header = F,stringsAsFactors = T)
sigD <- read.csv(paste0(ddir,"LME_fosi_inputs_sigDet.csv"),sep=",",header = F,stringsAsFactors = T)
sigB <- read.csv(paste0(ddir,"LME_fosi_inputs_sigTB.csv"),sep=",",header = F,stringsAsFactors = T)

names(sigM) <- c("LME","climate","lag","rZ","pZ")#,"type")
names(sigL) <- c("LME","climate","lag","rL","pL")#,"type")
names(sigP) <- c("LME","climate","lag","rP","pP")#,"type")
names(sigD) <- c("LME","climate","lag","rD","pD")#,"type")
names(sigB) <- c("LME","climate","lag","rB","pB")#,"type")

sig1 <- merge(sigL,sigM,by=c("LME","climate","lag"),all=T)
sig4 <- merge(sig1,sigP,by=c("LME","climate","lag"),all=T)
sig5 <- merge(sig4,sigB,by=c("LME","climate","lag"),all=T)
sigALL <- merge(sig5,sigD,by=c("LME","climate","lag"),all=T)


sigCHK <- subset(sigALL, LME=="CHK")
sigEBS <- subset(sigALL, LME=="EBS")
sigGAK <- subset(sigALL, LME=="GAK")
sigCCE <- subset(sigALL, LME=="CCE")
sigHI <- subset(sigALL, LME=="HI")
sigGMX <- subset(sigALL, LME=="GMX")
sigSE <- subset(sigALL, LME=="SE")
sigNE <- subset(sigALL, LME=="NE")

write.table(sigCHK,paste0(ddir,"LME_fosi_inputs_sig_climate_corrs_CHK.csv"),sep=",",row.names=F)
write.table(sigEBS,paste0(ddir,"LME_fosi_inputs_sig_climate_corrs_EBS.csv"),sep=",",row.names=F)
write.table(sigGAK,paste0(ddir,"LME_fosi_inputs_sig_climate_corrs_GAK.csv"),sep=",",row.names=F)
write.table(sigCCE,paste0(ddir,"LME_fosi_inputs_sig_climate_corrs_CCE.csv"),sep=",",row.names=F)
write.table(sigHI,paste0(ddir,"LME_fosi_inputs_sig_climate_corrs_HI.csv"),sep=",",row.names=F)
write.table(sigGMX,paste0(ddir,"LME_fosi_inputs_sig_climate_corrs_GMX.csv"),sep=",",row.names=F)
write.table(sigSE,paste0(ddir,"LME_fosi_inputs_sig_climate_corrs_SE.csv"),sep=",",row.names=F)
write.table(sigNE,paste0(ddir,"LME_fosi_inputs_sig_climate_corrs_NE.csv"),sep=",",row.names=F)

### ------------------------------------------------------------------------ ###
### ----------------------------- Heatmaps --------------------------------- ###
### ------------------------------------------------------------------------ ###
library(cowplot) #plot_grid

### Rearrange =====================================================
## CCE
cceNino34 <- subset(sigCCE, climate=="Nino34")
cceNino34 <- cceNino34[,c("lag","rP","rZ","rL","rD","rB")]
names(cceNino34) <- c("lag","TP","LZ","LZm","Det","TB")
cceNino34$lag <- as.factor(cceNino34$lag)
rcceNino34 <- melt(cceNino34)
names(rcceNino34) <- c("Lag","Type","Corr")

cceNino4 <- subset(sigCCE, climate=="Nino4")
cceNino4 <- cceNino4[,c("lag","rP","rZ","rL","rD","rB")]
names(cceNino4) <- c("lag","TP","LZ","LZm","Det","TB")
cceNino4$lag <- as.factor(cceNino4$lag)
rcceNino4 <- melt(cceNino4)
names(rcceNino4) <- c("Lag","Type","Corr")

cceNOI <- subset(sigCCE, climate=="NOI")
cceNOI <- cceNOI[,c("lag","rP","rZ","rL","rD","rB")]
names(cceNOI) <- c("lag","TP","LZ","LZm","Det","TB")
cceNOI$lag <- as.factor(cceNOI$lag)
rcceNOI <- melt(cceNOI)
names(rcceNOI) <- c("Lag","Type","Corr")

ccePDO <- subset(sigCCE, climate=="PDO")
ccePDO <- ccePDO[,c("lag","rP","rZ","rL","rD","rB")]
names(ccePDO) <- c("lag","TP","LZ","LZm","Det","TB")
ccePDO$lag <- as.factor(ccePDO$lag)
rccePDO <- melt(ccePDO)
names(rccePDO) <- c("Lag","Type","Corr")

cceSOI <- subset(sigCCE, climate=="SOI")
cceSOI <- cceSOI[,c("lag","rP","rZ","rL","rD","rB")]
names(cceSOI) <- c("lag","TP","LZ","LZm","Det","TB")
cceSOI$lag <- as.factor(cceSOI$lag)
rcceSOI <- melt(cceSOI)
names(rcceSOI) <- c("Lag","Type","Corr")


## EBS
ebsNino34 <- subset(sigEBS, climate=="Nino34")
ebsNino34 <- ebsNino34[,c("lag","rP","rZ","rL","rD","rB")]
names(ebsNino34) <- c("lag","TP","LZ","LZm","Det","TB")
ebsNino34$lag <- as.factor(ebsNino34$lag)
rebsNino34 <- melt(ebsNino34)
names(rebsNino34) <- c("Lag","Type","Corr")

ebsNino4 <- subset(sigEBS, climate=="Nino4")
ebsNino4 <- ebsNino4[,c("lag","rP","rZ","rL","rD","rB")]
names(ebsNino4) <- c("lag","TP","LZ","LZm","Det","TB")
ebsNino4$lag <- as.factor(ebsNino4$lag)
rebsNino4 <- melt(ebsNino4)
names(rebsNino4) <- c("Lag","Type","Corr")

ebsNOI <- subset(sigEBS, climate=="NOI")
ebsNOI <- ebsNOI[,c("lag","rP","rZ","rL","rD","rB")]
names(ebsNOI) <- c("lag","TP","LZ","LZm","Det","TB")
ebsNOI$lag <- as.factor(ebsNOI$lag)
rebsNOI <- melt(ebsNOI)
names(rebsNOI) <- c("Lag","Type","Corr")

ebsPDO <- subset(sigEBS, climate=="PDO")
ebsPDO <- ebsPDO[,c("lag","rP","rZ","rL","rD","rB")]
names(ebsPDO) <- c("lag","TP","LZ","LZm","Det","TB")
ebsPDO$lag <- as.factor(ebsPDO$lag)
rebsPDO <- melt(ebsPDO)
names(rebsPDO) <- c("Lag","Type","Corr")

ebsSOI <- subset(sigEBS, climate=="SOI")
ebsSOI <- ebsSOI[,c("lag","rP","rZ","rL","rD","rB")]
names(ebsSOI) <- c("lag","TP","LZ","LZm","Det","TB")
ebsSOI$lag <- as.factor(ebsSOI$lag)
rebsSOI <- melt(ebsSOI)
names(rebsSOI) <- c("Lag","Type","Corr")

ebsMEI <- subset(sigEBS, climate=="MEI")
ebsMEI <- ebsMEI[,c("lag","rP","rZ","rL","rD","rB")]
names(ebsMEI) <- c("lag","TP","LZ","LZm","Det","TB")
ebsMEI$lag <- as.factor(ebsMEI$lag)
rebsMEI <- melt(ebsMEI)
names(rebsMEI) <- c("Lag","Type","Corr")


## GAK
gakMEI <- subset(sigGAK, climate=="MEI")
gakMEI <- gakMEI[,c("lag","rP","rZ","rL","rD","rB")]
names(gakMEI) <- c("lag","TP","LZ","LZm","Det","TB")
gakMEI$lag <- as.factor(gakMEI$lag)
rgakMEI <- melt(gakMEI)
names(rgakMEI) <- c("Lag","Type","Corr")

gakNino3 <- subset(sigGAK, climate=="Nino3")
gakNino3 <- gakNino3[,c("lag","rP","rZ","rL","rD","rB")]
names(gakNino3) <- c("lag","TP","LZ","LZm","Det","TB")
gakNino3$lag <- as.factor(gakNino3$lag)
rgakNino3 <- melt(gakNino3)
names(rgakNino3) <- c("Lag","Type","Corr")

gakNino34 <- subset(sigGAK, climate=="Nino34")
gakNino34 <- gakNino34[,c("lag","rP","rZ","rL","rD","rB")]
names(gakNino34) <- c("lag","TP","LZ","LZm","Det","TB")
gakNino34$lag <- as.factor(gakNino34$lag)
rgakNino34 <- melt(gakNino34)
names(rgakNino34) <- c("Lag","Type","Corr")

gakNino4 <- subset(sigGAK, climate=="Nino4")
gakNino4 <- gakNino4[,c("lag","rP","rZ","rL","rD","rB")]
names(gakNino4) <- c("lag","TP","LZ","LZm","Det","TB")
gakNino4$lag <- as.factor(gakNino4$lag)
rgakNino4 <- melt(gakNino4)
names(rgakNino4) <- c("Lag","Type","Corr")

gakNOI <- subset(sigGAK, climate=="NOI")
gakNOI <- gakNOI[,c("lag","rP","rZ","rL","rD","rB")]
names(gakNOI) <- c("lag","TP","LZ","LZm","Det","TB")
gakNOI$lag <- as.factor(gakNOI$lag)
rgakNOI <- melt(gakNOI)
names(rgakNOI) <- c("Lag","Type","Corr")

gakPDO <- subset(sigGAK, climate=="PDO")
gakPDO <- gakPDO[,c("lag","rP","rZ","rL","rD","rB")]
names(gakPDO) <- c("lag","TP","LZ","LZm","Det","TB")
gakPDO$lag <- as.factor(gakPDO$lag)
rgakPDO <- melt(gakPDO)
names(rgakPDO) <- c("Lag","Type","Corr")


## HI
hiNino3 <- subset(sigHI, climate=="Nino3")
hiNino3 <- hiNino3[,c("lag","rP","rZ","rL","rD","rB")]
names(hiNino3) <- c("lag","TP","LZ","LZm","Det","TB")
hiNino3$lag <- as.factor(hiNino3$lag)
rhiNino3 <- melt(hiNino3)
names(rhiNino3) <- c("Lag","Type","Corr")

hiNino34 <- subset(sigHI, climate=="Nino34")
hiNino34 <- hiNino34[,c("lag","rP","rZ","rL","rD","rB")]
names(hiNino34) <- c("lag","TP","LZ","LZm","Det","TB")
hiNino34$lag <- as.factor(hiNino34$lag)
rhiNino34 <- melt(hiNino34)
names(rhiNino34) <- c("Lag","Type","Corr")

hiNino4 <- subset(sigHI, climate=="Nino4")
hiNino4 <- hiNino4[,c("lag","rP","rZ","rL","rD","rB")]
names(hiNino4) <- c("lag","TP","LZ","LZm","Det","TB")
hiNino4$lag <- as.factor(hiNino4$lag)
rhiNino4 <- melt(hiNino4)
names(rhiNino4) <- c("Lag","Type","Corr")

hiNOI <- subset(sigHI, climate=="NOI")
hiNOI <- hiNOI[,c("lag","rP","rZ","rL","rD","rB")]
names(hiNOI) <- c("lag","TP","LZ","LZm","Det","TB")
hiNOI$lag <- as.factor(hiNOI$lag)
rhiNOI <- melt(hiNOI)
names(rhiNOI) <- c("Lag","Type","Corr")

hiPDO <- subset(sigHI, climate=="PDO")
hiPDO <- hiPDO[,c("lag","rP","rZ","rL","rD","rB")]
names(hiPDO) <- c("lag","TP","LZ","LZm","Det","TB")
hiPDO$lag <- as.factor(hiPDO$lag)
rhiPDO <- melt(hiPDO)
names(rhiPDO) <- c("Lag","Type","Corr")


## CHK
chkAMO <- subset(sigCHK, climate=="AMO")
chkAMO <- chkAMO[,c("lag","rP","rZ","rL","rD","rB")]
names(chkAMO) <- c("lag","TP","LZ","LZm","Det","TB")
chkAMO$lag <- as.factor(chkAMO$lag)
rchkAMO <- melt(chkAMO)
names(rchkAMO) <- c("Lag","Type","Corr")

chkMEI <- subset(sigCHK, climate=="MEI")
chkMEI <- chkMEI[,c("lag","rP","rZ","rL","rD","rB")]
names(chkMEI) <- c("lag","TP","LZ","LZm","Det","TB")
chkMEI$lag <- as.factor(chkMEI$lag)
rchkMEI <- melt(chkMEI)
names(rchkMEI) <- c("Lag","Type","Corr")

chkPDO <- subset(sigCHK, climate=="PDO")
chkPDO <- chkPDO[,c("lag","rP","rZ","rL","rD","rB")]
names(chkPDO) <- c("lag","TP","LZ","LZm","Det","TB")
chkPDO$lag <- as.factor(chkPDO$lag)
rchkPDO <- melt(chkPDO)
names(rchkPDO) <- c("Lag","Type","Corr")


##GMX
gmxPDO <- subset(sigGMX, climate=="PDO")
gmxPDO <- gmxPDO[,c("lag","rP","rZ","rL","rD","rB")]
names(gmxPDO) <- c("lag","TP","LZ","LZm","Det","TB")
gmxPDO$lag <- as.factor(gmxPDO$lag)
rgmxPDO <- melt(gmxPDO)
names(rgmxPDO) <- c("Lag","Type","Corr")

gmxNOI <- subset(sigGMX, climate=="NOI")
gmxNOI <- gmxNOI[,c("lag","rP","rZ","rL","rD","rB")]
names(gmxNOI) <- c("lag","TP","LZ","LZm","Det","TB")
gmxNOI$lag <- as.factor(gmxNOI$lag)
rgmxNOI <- melt(gmxNOI)
names(rgmxNOI) <- c("Lag","Type","Corr")

gmxAMO <- subset(sigGMX, climate=="AMO")
gmxAMO <- gmxAMO[,c("lag","rP","rZ","rL","rD","rB")]
names(gmxAMO) <- c("lag","TP","LZ","LZm","Det","TB")
gmxAMO$lag <- as.factor(gmxAMO$lag)
rgmxAMO <- melt(gmxAMO)
names(rgmxAMO) <- c("Lag","Type","Corr")


## NE
neAMO <- subset(sigNE, climate=="AMO")
neAMO <- neAMO[,c("lag","rP","rZ","rL","rD","rB")]
names(neAMO) <- c("lag","TP","LZ","LZm","Det","TB")
neAMO$lag <- as.factor(neAMO$lag)
rneAMO <- melt(neAMO)
names(rneAMO) <- c("Lag","Type","Corr")

neNAO <- subset(sigNE, climate=="NAO")
neNAO <- neNAO[,c("lag","rP","rZ","rL","rD","rB")]
names(neNAO) <- c("lag","TP","LZ","LZm","Det","TB")
neNAO$lag <- as.factor(neNAO$lag)
rneNAO <- melt(neNAO)
names(rneNAO) <- c("Lag","Type","Corr")


##SE
seNAO <- subset(sigSE, climate=="NAO")
seNAO <- seNAO[,c("lag","rP","rZ","rL","rD","rB")]
names(seNAO) <- c("lag","TP","LZ","LZm","Det","TB")
seNAO$lag <- as.factor(seNAO$lag)
rseNAO <- melt(seNAO)
names(rseNAO) <- c("Lag","Type","Corr")

seAMO <- subset(sigSE, climate=="AMO")
seAMO <- seAMO[,c("lag","rP","rZ","rL","rD","rB")]
names(seAMO) <- c("lag","TP","LZ","LZm","Det","TB")
seAMO$lag <- as.factor(seAMO$lag)
rseAMO <- melt(seAMO)
names(rseAMO) <- c("Lag","Type","Corr")

sePDO <- subset(sigSE, climate=="PDO")
sePDO <- sePDO[,c("lag","rP","rZ","rL","rD","rB")]
names(sePDO) <- c("lag","TP","LZ","LZm","Det","TB")
sePDO$lag <- as.factor(sePDO$lag)
rsePDO <- melt(sePDO)
names(rsePDO) <- c("Lag","Type","Corr")




###====== Plots ===================================================================
## CCE
c1 <- ggplot(data = rcceNino34, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino34\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

c2 <- ggplot(data = rcceNino4, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino4\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

c3 <- ggplot(data = rcceNOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

c4 <- ggplot(data = rccePDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

c5 <- ggplot(data = rcceSOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="SOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 


png(paste0(figp,'Heatmaps_CCE_corr_clims_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( c1,c2,c3,c4,c5,
           nrow = 2, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## EBS
e0 <- ggplot(data = rebsMEI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="MEI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e1 <- ggplot(data = rebsNino34, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino34\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e2 <- ggplot(data = rebsNino4, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino4\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e3 <- ggplot(data = rebsNOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e4 <- ggplot(data = rebsPDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

e5 <- ggplot(data = rebsSOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="SOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 


png(paste0(figp,'Heatmaps_EBS_corr_clims_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( e0,e1,e2,e3,e4,e5,
           nrow = 2, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## GAK
g0 <- ggplot(data = rgakNino3, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino3\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g1 <- ggplot(data = rgakNino34, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino34\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g2 <- ggplot(data = rgakNino4, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino4\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g3 <- ggplot(data = rgakNOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g4 <- ggplot(data = rgakPDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

g5 <- ggplot(data = rgakMEI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="MEI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 


png(paste0(figp,'Heatmaps_GAK_corr_clims_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( g0,g1,g2,g3,g4,g5,
           nrow = 2, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## HI
h0 <- ggplot(data = rhiNino3, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino3\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

h1 <- ggplot(data = rhiNino34, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino34\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

h2 <- ggplot(data = rhiNino4, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="Nino4\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

h3 <- ggplot(data = rhiNOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_HI_corr_clims_inputs.png'), 
    width = 6*300,        # 5 x 300 pixels
    height = 6*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( h0,h1,h2,h3,
           nrow = 2, ncol = 2,
           rel_widths = c(1,1), rel_heights = c(1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## CHK
k0 <- ggplot(data = rchkMEI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="MEI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

k1 <- ggplot(data = rchkPDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

k2 <- ggplot(data = rchkAMO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AMO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3)

png(paste0(figp,'Heatmaps_CHK_corr_clims_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 3*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( k0,k1,k2,
           nrow = 1, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## GMX
m1 <- ggplot(data = rgmxAMO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AMO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

m3 <- ggplot(data = rgmxPDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

m2 <- ggplot(data = rgmxNOI, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NOI\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_GMX_corr_clims_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 3*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( m1,m2,m3,
           nrow = 1, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## SE
s1 <- ggplot(data = rseAMO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AMO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

s2 <- ggplot(data = rseNAO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NAO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

s3 <- ggplot(data = rsePDO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="PDO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

png(paste0(figp,'Heatmaps_SE_corr_clims_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 3*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( s1,s2,s3,
           nrow = 1, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()


## NE
n1 <- ggplot(data = rneAMO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="AMO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 

n2 <- ggplot(data = rneNAO, aes(y=Type, x=Lag, fill=Corr)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", limit = c(-0.75,0.75),  
                       name="NAO\nPearson\nCorrelation") +
  theme_minimal() + labs(x="")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + 
  geom_text(aes(Lag, Type, label = signif(Corr,digits = 2)), color = "black", size = 3) 


png(paste0(figp,'Heatmaps_NE_corr_clims_inputs.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 3*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( n1,n2,
           nrow = 1, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()

