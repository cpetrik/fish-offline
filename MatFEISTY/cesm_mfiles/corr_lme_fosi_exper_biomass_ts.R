# Find patterns in climate correlations

rm(list=ls())

source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/"

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
ddir <- "/Volumes/MIP/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/"
sigS <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigS.csv"),sep=",",header = F,stringsAsFactors = T)
sigM <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigM.csv"),sep=",",header = F,stringsAsFactors = T)
sigL <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigL.csv"),sep=",",header = F,stringsAsFactors = T)
sigF <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigF.csv"),sep=",",header = F,stringsAsFactors = T)
sigP <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigP.csv"),sep=",",header = F,stringsAsFactors = T)
sigD <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigD.csv"),sep=",",header = F,stringsAsFactors = T)
sigA <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigA.csv"),sep=",",header = F,stringsAsFactors = T)
sigB <- read.csv(paste0(ddir,"LME_fosi_fished_v14_All_fish03_sigB.csv"),sep=",",header = F,stringsAsFactors = T)

# sigS[6] <- "S"
# sigM[6] <- "M"
# sigL[6] <- "L"
# sigF[6] <- "F"
# sigP[6] <- "P"
# sigD[6] <- "D"
# sigA[6] <- "A"
# sigB[6] <- "B"

names(sigS) <- c("LME","climate","lag","rS","pS")#,"type")
names(sigM) <- c("LME","climate","lag","rM","pM")#,"type")
names(sigL) <- c("LME","climate","lag","rL","pL")#,"type")
names(sigF) <- c("LME","climate","lag","rF","pF")#,"type")
names(sigP) <- c("LME","climate","lag","rP","pP")#,"type")
names(sigD) <- c("LME","climate","lag","rD","pD")#,"type")
names(sigA) <- c("LME","climate","lag","rA","pA")#,"type")
names(sigB) <- c("LME","climate","lag","rB","pB")#,"type")

sig1 <- merge(sigS,sigM,by=c("LME","climate","lag"),all=T)
sig2 <- merge(sig1,sigL,by=c("LME","climate","lag"),all=T)
sig3 <- merge(sig2,sigF,by=c("LME","climate","lag"),all=T)
sig4 <- merge(sig3,sigP,by=c("LME","climate","lag"),all=T)
sig5 <- merge(sig4,sigD,by=c("LME","climate","lag"),all=T)
sig6 <- merge(sig5,sigA,by=c("LME","climate","lag"),all=T)
sigALL <- merge(sig6,sigB,by=c("LME","climate","lag"),all=T)


sigCHK <- subset(sigALL, LME=="CHK")
sigEBS <- subset(sigALL, LME=="EBS")
sigGAK <- subset(sigALL, LME=="GAK")
sigCCE <- subset(sigALL, LME=="CCE")
sigHI <- subset(sigALL, LME=="HI")
sigGMX <- subset(sigALL, LME=="GMX")
sigSE <- subset(sigALL, LME=="SE")
sigNE <- subset(sigALL, LME=="NE")

write.table(sigCHK,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_climate_corrs_CHK.csv"),sep=",",row.names=F)
write.table(sigEBS,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_climate_corrs_EBS.csv"),sep=",",row.names=F)
write.table(sigGAK,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_climate_corrs_GAK.csv"),sep=",",row.names=F)
write.table(sigCCE,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_climate_corrs_CCE.csv"),sep=",",row.names=F)
write.table(sigHI,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_climate_corrs_HI.csv"),sep=",",row.names=F)
write.table(sigGMX,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_climate_corrs_GMX.csv"),sep=",",row.names=F)
write.table(sigSE,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_climate_corrs_SE.csv"),sep=",",row.names=F)
write.table(sigNE,paste0(ddir,"LME_fosi_fished_v14_All_fish03_sig_climate_corrs_NE.csv"),sep=",",row.names=F)

### ------------------------------------------------------------------------ ###
### ----------------------------- Heatmaps --------------------------------- ###
### ------------------------------------------------------------------------ ###
library(cowplot) #plot_grid

### Rearrange =====================================================
## CCE
cceNino34 <- subset(sigCCE, climate=="Nino34")
cceNino34 <- cceNino34[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(cceNino34) <- c("lag","B","S","M","L","F","P","D","A")
cceNino34$lag <- as.factor(cceNino34$lag)
rcceNino34 <- melt(cceNino34)
names(rcceNino34) <- c("Lag","Type","Corr")

cceNino4 <- subset(sigCCE, climate=="Nino4")
cceNino4 <- cceNino4[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(cceNino4) <- c("lag","B","S","M","L","F","P","D","A")
cceNino4$lag <- as.factor(cceNino4$lag)
rcceNino4 <- melt(cceNino4)
names(rcceNino4) <- c("Lag","Type","Corr")

cceNOI <- subset(sigCCE, climate=="NOI")
cceNOI <- cceNOI[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(cceNOI) <- c("lag","B","S","M","L","F","P","D","A")
cceNOI$lag <- as.factor(cceNOI$lag)
rcceNOI <- melt(cceNOI)
names(rcceNOI) <- c("Lag","Type","Corr")

ccePDO <- subset(sigCCE, climate=="PDO")
ccePDO <- ccePDO[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(ccePDO) <- c("lag","B","S","M","L","F","P","D","A")
ccePDO$lag <- as.factor(ccePDO$lag)
rccePDO <- melt(ccePDO)
names(rccePDO) <- c("Lag","Type","Corr")

cceSOI <- subset(sigCCE, climate=="SOI")
cceSOI <- cceSOI[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(cceSOI) <- c("lag","B","S","M","L","F","P","D","A")
cceSOI$lag <- as.factor(cceSOI$lag)
rcceSOI <- melt(cceSOI)
names(rcceSOI) <- c("Lag","Type","Corr")


## EBS
ebsNino34 <- subset(sigEBS, climate=="Nino34")
ebsNino34 <- ebsNino34[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(ebsNino34) <- c("lag","B","S","M","L","F","P","D","A")
ebsNino34$lag <- as.factor(ebsNino34$lag)
rebsNino34 <- melt(ebsNino34)
names(rebsNino34) <- c("Lag","Type","Corr")

ebsNino4 <- subset(sigEBS, climate=="Nino4")
ebsNino4 <- ebsNino4[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(ebsNino4) <- c("lag","B","S","M","L","F","P","D","A")
ebsNino4$lag <- as.factor(ebsNino4$lag)
rebsNino4 <- melt(ebsNino4)
names(rebsNino4) <- c("Lag","Type","Corr")

ebsNOI <- subset(sigEBS, climate=="NOI")
ebsNOI <- ebsNOI[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(ebsNOI) <- c("lag","B","S","M","L","F","P","D","A")
ebsNOI$lag <- as.factor(ebsNOI$lag)
rebsNOI <- melt(ebsNOI)
names(rebsNOI) <- c("Lag","Type","Corr")

ebsPDO <- subset(sigEBS, climate=="PDO")
ebsPDO <- ebsPDO[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(ebsPDO) <- c("lag","B","S","M","L","F","P","D","A")
ebsPDO$lag <- as.factor(ebsPDO$lag)
rebsPDO <- melt(ebsPDO)
names(rebsPDO) <- c("Lag","Type","Corr")

ebsSOI <- subset(sigEBS, climate=="SOI")
ebsSOI <- ebsSOI[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(ebsSOI) <- c("lag","B","S","M","L","F","P","D","A")
ebsSOI$lag <- as.factor(ebsSOI$lag)
rebsSOI <- melt(ebsSOI)
names(rebsSOI) <- c("Lag","Type","Corr")

ebsMEI <- subset(sigEBS, climate=="MEI")
ebsMEI <- ebsMEI[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(ebsMEI) <- c("lag","B","S","M","L","F","P","D","A")
ebsMEI$lag <- as.factor(ebsMEI$lag)
rebsMEI <- melt(ebsMEI)
names(rebsMEI) <- c("Lag","Type","Corr")


## GAK
gakMEI <- subset(sigGAK, climate=="MEI")
gakMEI <- gakMEI[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gakMEI) <- c("lag","B","S","M","L","F","P","D","A")
gakMEI$lag <- as.factor(gakMEI$lag)
rgakMEI <- melt(gakMEI)
names(rgakMEI) <- c("Lag","Type","Corr")

gakNino3 <- subset(sigGAK, climate=="Nino3")
gakNino3 <- gakNino3[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gakNino3) <- c("lag","B","S","M","L","F","P","D","A")
gakNino3$lag <- as.factor(gakNino3$lag)
rgakNino3 <- melt(gakNino3)
names(rgakNino3) <- c("Lag","Type","Corr")

gakNino34 <- subset(sigGAK, climate=="Nino34")
gakNino34 <- gakNino34[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gakNino34) <- c("lag","B","S","M","L","F","P","D","A")
gakNino34$lag <- as.factor(gakNino34$lag)
rgakNino34 <- melt(gakNino34)
names(rgakNino34) <- c("Lag","Type","Corr")

gakNino4 <- subset(sigGAK, climate=="Nino4")
gakNino4 <- gakNino4[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gakNino4) <- c("lag","B","S","M","L","F","P","D","A")
gakNino4$lag <- as.factor(gakNino4$lag)
rgakNino4 <- melt(gakNino4)
names(rgakNino4) <- c("Lag","Type","Corr")

gakNOI <- subset(sigGAK, climate=="NOI")
gakNOI <- gakNOI[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gakNOI) <- c("lag","B","S","M","L","F","P","D","A")
gakNOI$lag <- as.factor(gakNOI$lag)
rgakNOI <- melt(gakNOI)
names(rgakNOI) <- c("Lag","Type","Corr")

gakPDO <- subset(sigGAK, climate=="PDO")
gakPDO <- gakPDO[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gakPDO) <- c("lag","B","S","M","L","F","P","D","A")
gakPDO$lag <- as.factor(gakPDO$lag)
rgakPDO <- melt(gakPDO)
names(rgakPDO) <- c("Lag","Type","Corr")


## HI
hiNino3 <- subset(sigHI, climate=="Nino3")
hiNino3 <- hiNino3[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(hiNino3) <- c("lag","B","S","M","L","F","P","D","A")
hiNino3$lag <- as.factor(hiNino3$lag)
rhiNino3 <- melt(hiNino3)
names(rhiNino3) <- c("Lag","Type","Corr")

hiNino34 <- subset(sigHI, climate=="Nino34")
hiNino34 <- hiNino34[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(hiNino34) <- c("lag","B","S","M","L","F","P","D","A")
hiNino34$lag <- as.factor(hiNino34$lag)
rhiNino34 <- melt(hiNino34)
names(rhiNino34) <- c("Lag","Type","Corr")

hiNino4 <- subset(sigHI, climate=="Nino4")
hiNino4 <- hiNino4[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(hiNino4) <- c("lag","B","S","M","L","F","P","D","A")
hiNino4$lag <- as.factor(hiNino4$lag)
rhiNino4 <- melt(hiNino4)
names(rhiNino4) <- c("Lag","Type","Corr")

hiNOI <- subset(sigHI, climate=="NOI")
hiNOI <- hiNOI[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(hiNOI) <- c("lag","B","S","M","L","F","P","D","A")
hiNOI$lag <- as.factor(hiNOI$lag)
rhiNOI <- melt(hiNOI)
names(rhiNOI) <- c("Lag","Type","Corr")

hiPDO <- subset(sigHI, climate=="PDO")
hiPDO <- hiPDO[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(hiPDO) <- c("lag","B","S","M","L","F","P","D","A")
hiPDO$lag <- as.factor(hiPDO$lag)
rhiPDO <- melt(hiPDO)
names(rhiPDO) <- c("Lag","Type","Corr")


## CHK
chkAMO <- subset(sigCHK, climate=="AMO")
chkAMO <- chkAMO[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(chkAMO) <- c("lag","B","S","M","L","F","P","D","A")
chkAMO$lag <- as.factor(chkAMO$lag)
rchkAMO <- melt(chkAMO)
names(rchkAMO) <- c("Lag","Type","Corr")

chkMEI <- subset(sigCHK, climate=="MEI")
chkMEI <- chkMEI[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(chkMEI) <- c("lag","B","S","M","L","F","P","D","A")
chkMEI$lag <- as.factor(chkMEI$lag)
rchkMEI <- melt(chkMEI)
names(rchkMEI) <- c("Lag","Type","Corr")

chkPDO <- subset(sigCHK, climate=="PDO")
chkPDO <- chkPDO[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(chkPDO) <- c("lag","B","S","M","L","F","P","D","A")
chkPDO$lag <- as.factor(chkPDO$lag)
rchkPDO <- melt(chkPDO)
names(rchkPDO) <- c("Lag","Type","Corr")


##GMX
gmxPDO <- subset(sigGMX, climate=="PDO")
gmxPDO <- gmxPDO[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gmxPDO) <- c("lag","B","S","M","L","F","P","D","A")
gmxPDO$lag <- as.factor(gmxPDO$lag)
rgmxPDO <- melt(gmxPDO)
names(rgmxPDO) <- c("Lag","Type","Corr")

gmxNOI <- subset(sigGMX, climate=="NOI")
gmxNOI <- gmxNOI[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gmxNOI) <- c("lag","B","S","M","L","F","P","D","A")
gmxNOI$lag <- as.factor(gmxNOI$lag)
rgmxNOI <- melt(gmxNOI)
names(rgmxNOI) <- c("Lag","Type","Corr")

gmxAMO <- subset(sigGMX, climate=="AMO")
gmxAMO <- gmxAMO[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(gmxAMO) <- c("lag","B","S","M","L","F","P","D","A")
gmxAMO$lag <- as.factor(gmxAMO$lag)
rgmxAMO <- melt(gmxAMO)
names(rgmxAMO) <- c("Lag","Type","Corr")


## NE
neAMO <- subset(sigNE, climate=="AMO")
neAMO <- neAMO[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(neAMO) <- c("lag","B","S","M","L","F","P","D","A")
neAMO$lag <- as.factor(neAMO$lag)
rneAMO <- melt(neAMO)
names(rneAMO) <- c("Lag","Type","Corr")

neNAO <- subset(sigNE, climate=="NAO")
neNAO <- neNAO[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(neNAO) <- c("lag","B","S","M","L","F","P","D","A")
neNAO$lag <- as.factor(neNAO$lag)
rneNAO <- melt(neNAO)
names(rneNAO) <- c("Lag","Type","Corr")


##SE
seNAO <- subset(sigSE, climate=="NAO")
seNAO <- seNAO[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(seNAO) <- c("lag","B","S","M","L","F","P","D","A")
seNAO$lag <- as.factor(seNAO$lag)
rseNAO <- melt(seNAO)
names(rseNAO) <- c("Lag","Type","Corr")

seAMO <- subset(sigSE, climate=="AMO")
seAMO <- seAMO[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(seAMO) <- c("lag","B","S","M","L","F","P","D","A")
seAMO$lag <- as.factor(seAMO$lag)
rseAMO <- melt(seAMO)
names(rseAMO) <- c("Lag","Type","Corr")

sePDO <- subset(sigSE, climate=="PDO")
sePDO <- sePDO[,c("lag","rB","rS","rM","rL","rF","rP","rD","rA")]
names(sePDO) <- c("lag","B","S","M","L","F","P","D","A")
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


png(paste0(figp,'Heatmaps_CCE_corr_clims_fish.png'), 
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


png(paste0(figp,'Heatmaps_EBS_corr_clims_fish.png'), 
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


png(paste0(figp,'Heatmaps_GAK_corr_clims_fish.png'), 
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

png(paste0(figp,'Heatmaps_HI_corr_clims_fish.png'), 
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

png(paste0(figp,'Heatmaps_CHK_corr_clims_fish.png'), 
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

png(paste0(figp,'Heatmaps_GMX_corr_clims_fish.png'), 
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

png(paste0(figp,'Heatmaps_SE_corr_clims_fish.png'), 
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


png(paste0(figp,'Heatmaps_NE_corr_clims_fish.png'), 
    width = 10*300,        # 5 x 300 pixels
    height = 3*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)
plot_grid( n1,n2,
           nrow = 1, ncol = 3,
           rel_widths = c(1,1), rel_heights = c(1,1,1) ,
           align = 'h' )#, labels = "auto", label_size = 12, hjust = -4)
dev.off()

