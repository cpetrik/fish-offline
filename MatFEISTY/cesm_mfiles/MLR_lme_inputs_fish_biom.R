# Mult linear regression of forcing on fish biomass

rm(list=ls())

library("MuMIn")


source(file = "/Users/cpetrik/Dropbox/UCSC/CVC-LCM/CVC_data/outmigration_model/HighstatLibV6.r")
setwd("/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/cesm_mfiles/")
figp <- "/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/corrs/"

### ------------------------------------------------------------
# load data
datap <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/FOSI/"
TP <- read.csv(paste0(datap,"LMEs_anoms_div2SD_TP.csv"),sep=",",header = T,stringsAsFactors = F)
TB <- read.csv(paste0(datap,"LMEs_anoms_div2SD_TB.csv"),sep=",",header = T,stringsAsFactors = F)
Det <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Det.csv"),sep=",",header = T,stringsAsFactors = F)
Zmeso <- read.csv(paste0(datap,"LMEs_anoms_div2SD_Zmeso.csv"),sep=",",header = T,stringsAsFactors = F)
ZmLoss <- read.csv(paste0(datap,"LMEs_anoms_div2SD_ZmLoss.csv"),sep=",",header = T,stringsAsFactors = F)
# S <- read.csv(paste0(datap,"LMEs_anoms_div2SD_S.csv"),sep=",",header = T,stringsAsFactors = F)
# M <- read.csv(paste0(datap,"LMEs_anoms_div2SD_M.csv"),sep=",",header = T,stringsAsFactors = F)
# L <- read.csv(paste0(datap,"LMEs_anoms_div2SD_L.csv"),sep=",",header = T,stringsAsFactors = F)
FF <- read.csv(paste0(datap,"LMEs_anoms_div2SD_F.csv"),sep=",",header = T,stringsAsFactors = F)
P <- read.csv(paste0(datap,"LMEs_anoms_div2SD_P.csv"),sep=",",header = T,stringsAsFactors = F)
D <- read.csv(paste0(datap,"LMEs_anoms_div2SD_D.csv"),sep=",",header = T,stringsAsFactors = F)
A <- read.csv(paste0(datap,"LMEs_anoms_div2SD_A.csv"),sep=",",header = T,stringsAsFactors = F)
B <- read.csv(paste0(datap,"LMEs_anoms_div2SD_B.csv"),sep=",",header = T,stringsAsFactors = F)

datar <- "/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/regressions/"
Flag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_F.csv"),sep=",",header = T,stringsAsFactors = F)
Plag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_P.csv"),sep=",",header = T,stringsAsFactors = F)
Dlag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_D.csv"),sep=",",header = T,stringsAsFactors = F)
Alag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_A.csv"),sep=",",header = T,stringsAsFactors = F)
Blag <- read.csv(paste0(datar,"LMEs_regress_drivers_ALLdiv2SD_siglag_B.csv"),sep=",",header = T,stringsAsFactors = F)


### MLR of drivers ---------------------------------------------------------
dname <- c('TP','TB','Zmeso','ZmLoss','Det')

#Start with US LMEs
cname <- names(A)[2:64]
nlme <- length(cname)

#Loop over LMEs 
for (i in 1:nlme) {
  cn <- cname[i]
  
  ### Drivers
  drive <- TP[,c('Year',cn)]
  names(drive) <- c('Year','TP')
  drive$TB <- TB[,cn]
  drive$Det <- Det[,cn]
  drive$Zmeso <- Zmeso[,cn]
  drive$ZmLoss <- ZmLoss[,cn]
  
  
  ### F ----------------------------------------------------------------
  maxf <- max(Flag[i,])
  yst <- 1+maxf
  yen <- 68+maxf
  ffish <- data.frame(matrix(ncol = 6, nrow = yen))
  names(ffish) <- c('Fish','TP','TB','Det','Zmeso','ZmLoss')
  
  ffish$TP[(yst-Flag$TP[i]):(yen-Flag$TP[i])] <- drive$TP 
  ffish$TB[(yst-Flag$TB[i]):(yen-Flag$TB[i])] <- drive$TB 
  ffish$Det[(yst-Flag$Det[i]):(yen-Flag$Det[i])] <- drive$Det
  ffish$Zmeso[(yst-Flag$Zmeso[i]):(yen-Flag$Zmeso[i])] <- drive$Zmeso 
  ffish$ZmLoss[(yst-Flag$ZmLoss[i]):(yen-Flag$ZmLoss[i])] <- drive$ZmLoss
  ffish$Fish[(yst):(yen)] <- FF[,cn]
  ffish <- na.omit(ffish)
  
  options(na.action = "na.fail")
  fmod <- lm(Fish ~ TP + TB + Det + Zmeso + ZmLoss + TP*Zmeso + TP*ZmLoss + TB*Det, data=ffish)
  summary(fmod)
  fcombo <- dredge(fmod)
  ## Create arrays for coefficients & p-vals
  if (i==1) {
    fcoef <- data.frame(matrix(ncol = length(fcombo), nrow = nlme))
    names(fcoef) <- names(fcombo)
    fpval <- fcoef[,1:9]
    fcoef$LME <- cname
    fpval$LME <- cname
    pcoef <- fcoef
    ppval <- fpval
    dcoef <- fcoef
    dpval <- fpval
    acoef <- fcoef
    apval <- fpval
    bcoef <- fcoef
    bpval <- fpval
  }
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (fcombo)[1]
  fcoef[i, match(names(xx), colnames(fcoef))] = xx
  y <- summary(get.models(fcombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  fpval[i, match(names(yy), colnames(fpval))] = yy
  
  ### P ----------------------------------------------------------------
  maxp <- max(Plag[i,])
  yst <- 1+maxp
  yen <- 68+maxp
  pfish <- data.frame(matrix(ncol = 6, nrow = yen))
  names(pfish) <- c('Fish','TP','TB','Det','Zmeso','ZmLoss')
  
  pfish$TP[(yst-Plag$TP[i]):(yen-Plag$TP[i])] <- drive$TP 
  pfish$TB[(yst-Plag$TB[i]):(yen-Plag$TB[i])] <- drive$TB 
  pfish$Det[(yst-Plag$Det[i]):(yen-Plag$Det[i])] <- drive$Det
  pfish$Zmeso[(yst-Plag$Zmeso[i]):(yen-Plag$Zmeso[i])] <- drive$Zmeso 
  pfish$ZmLoss[(yst-Plag$ZmLoss[i]):(yen-Plag$ZmLoss[i])] <- drive$ZmLoss
  pfish$Fish[(yst):(yen)] <- P[,cn]
  pfish <- na.omit(pfish)
  
  pmod <- lm(Fish ~ TP + TB + Det + Zmeso + ZmLoss + TP*Zmeso + TP*ZmLoss + TB*Det, data=pfish)
  pcombo <- dredge(pmod)
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (pcombo)[1]
  pcoef[i, match(names(xx), colnames(pcoef))] = xx
  y <- summary(get.models(pcombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  ppval[i, match(names(yy), colnames(ppval))] = yy
  
  
  ### D ----------------------------------------------------------------
  maxd <- max(Dlag[i,])
  yst <- 1+maxd
  yen <- 68+maxd
  dfish <- data.frame(matrix(ncol = 6, nrow = yen))
  names(dfish) <- c('Fish','TP','TB','Det','Zmeso','ZmLoss')
  
  dfish$TP[(yst-Dlag$TP[i]):(yen-Dlag$TP[i])] <- drive$TP 
  dfish$TB[(yst-Dlag$TB[i]):(yen-Dlag$TB[i])] <- drive$TB 
  dfish$Det[(yst-Dlag$Det[i]):(yen-Dlag$Det[i])] <- drive$Det
  dfish$Zmeso[(yst-Dlag$Zmeso[i]):(yen-Dlag$Zmeso[i])] <- drive$Zmeso 
  dfish$ZmLoss[(yst-Dlag$ZmLoss[i]):(yen-Dlag$ZmLoss[i])] <- drive$ZmLoss
  dfish$Fish[(yst):(yen)] <- D[,cn]
  dfish <- na.omit(dfish)
  
  dmod <- lm(Fish ~ TP + TB + Det + Zmeso + ZmLoss + TP*Zmeso + TP*ZmLoss + TB*Det, data=dfish)
  dcombo <- dredge(dmod)
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (dcombo)[1]
  dcoef[i, match(names(xx), colnames(dcoef))] = xx
  y <- summary(get.models(dcombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  dpval[i, match(names(yy), colnames(dpval))] = yy
  
  ### A ----------------------------------------------------------------
  maxa <- max(Alag[i,])
  yst <- 1+maxa
  yen <- 68+maxa
  afish <- data.frame(matrix(ncol = 6, nrow = yen))
  names(afish) <- c('Fish','TP','TB','Det','Zmeso','ZmLoss')
  
  afish$TP[(yst-Alag$TP[i]):(yen-Alag$TP[i])] <- drive$TP 
  afish$TB[(yst-Alag$TB[i]):(yen-Alag$TB[i])] <- drive$TB 
  afish$Det[(yst-Alag$Det[i]):(yen-Alag$Det[i])] <- drive$Det
  afish$Zmeso[(yst-Alag$Zmeso[i]):(yen-Alag$Zmeso[i])] <- drive$Zmeso 
  afish$ZmLoss[(yst-Alag$ZmLoss[i]):(yen-Alag$ZmLoss[i])] <- drive$ZmLoss
  afish$Fish[(yst):(yen)] <- A[,cn]
  afish <- na.omit(afish)
  
  amod <- lm(Fish ~ TP + TB + Det + Zmeso + ZmLoss + TP*Zmeso + TP*ZmLoss + TB*Det, data=afish)
  acombo <- dredge(amod)
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (acombo)[1]
  acoef[i, match(names(xx), colnames(acoef))] = xx
  y <- summary(get.models(acombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  apval[i, match(names(yy), colnames(apval))] = yy
  
  
  ### B ----------------------------------------------------------------
  maxb <- max(Blag[i,])
  yst <- 1+maxb
  yen <- 68+maxb
  bent <- data.frame(matrix(ncol = 6, nrow = yen))
  names(bent) <- c('Bent','TP','TB','Det','Zmeso','ZmLoss')
  
  bent$TP[(yst-Blag$TP[i]):(yen-Blag$TP[i])] <- drive$TP 
  bent$TB[(yst-Blag$TB[i]):(yen-Blag$TB[i])] <- drive$TB 
  bent$Det[(yst-Blag$Det[i]):(yen-Blag$Det[i])] <- drive$Det
  bent$Zmeso[(yst-Blag$Zmeso[i]):(yen-Blag$Zmeso[i])] <- drive$Zmeso 
  bent$ZmLoss[(yst-Blag$ZmLoss[i]):(yen-Blag$ZmLoss[i])] <- drive$ZmLoss
  bent$Bent[(yst):(yen)] <- B[,cn]
  bent <- na.omit(bent)
  
  bmod <- lm(Bent ~ TP + TB + Det + Zmeso + ZmLoss + TP*Zmeso + TP*ZmLoss + TB*Det, data=bent)
  bcombo <- dredge(bmod)
  ## Put coefficients & p-val of best model in arrays for saving
  xx <- (bcombo)[1]
  bcoef[i, match(names(xx), colnames(bcoef))] = xx
  y <- summary(get.models(bcombo, 1)[[1]])
  yy <- (y$coefficients[,4])
  bpval[i, match(names(yy), colnames(bpval))] = yy
  
  
} # LMEs

write.table(fcoef,paste0(datar,"LME_F_mlr_coeffs_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(pcoef,paste0(datar,"LME_P_mlr_coeffs_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(dcoef,paste0(datar,"LME_D_mlr_coeffs_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(acoef,paste0(datar,"LME_A_mlr_coeffs_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(bcoef,paste0(datar,"LME_B_mlr_coeffs_ALLdiv2SD.csv"),sep=",",row.names=F)

write.table(fpval,paste0(datar,"LME_F_mlr_pvals_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(ppval,paste0(datar,"LME_P_mlr_pvals_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(dpval,paste0(datar,"LME_D_mlr_pvals_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(apval,paste0(datar,"LME_A_mlr_pvals_ALLdiv2SD.csv"),sep=",",row.names=F)
write.table(bpval,paste0(datar,"LME_B_mlr_pvals_ALLdiv2SD.csv"),sep=",",row.names=F)



### Would heatmap or something help find patterns?
library(ggplot2)
library(cowplot) #plot_grid
library(gridExtra)
library(ggpubr) #cowplot and gridExtra
library(gplots)

Acoef <- acoef[,2:9]
Acoef$LME <- acoef$LME
test <- melt(Acoef)
names(test) <- c("LME","driver","coef")

a1 <- ggplot(data = test, aes(y=LME, x=driver, fill=coef)) + 
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", name="coef") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 8, hjust = 1)) +
  coord_fixed() + ggtitle("All fish") + 
  geom_text(aes(Lag, Type, label = sym), color = "black", size = 3) 

#scale_fill_distiller(palette = "RdBu", limit = c(-1.05,1.05),
#                     breaks=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5),
#                     labels=c(-1.5,-1.0,-0.5,0,0.5,1.0,1.5), name="Det\ncorr") +


# Cluster LMEs
library(rpart)
library(tree)
library(pvclust)
library(mclust)
library(dendextend)

vmat <- as.matrix(vec[2:39,c(2:6,8)])
rvec <- as.data.frame(t(vmat))
bvec <- rvec[,c(high,20,36)]
tvec <- bvec[1:5,params]
tvec <- as.data.frame(t(tvec))

tvec <- Acoef[,1:8]
rownames(tvec) <- Acoef[,9]
tvec[is.na(tvec)] <- 0

d_fish <- dist(tvec) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish, k=8) 
# We hang the dendrogram a bit:
#dfish <- hang.dendrogram(dfish,hang_height=0.1)
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

# And plot:
pdf(paste0(figp,"Cluster_LME_A_mlr_coeffs_ALLdiv2SD.pdf"))
par(mfrow = c(1,1))
par(mar = c(3,3,3,7))
plot(dfish, 
     main = "hclust complete", 
     horiz =  TRUE,  nodePar = list(cex = .007))
dev.off()

some_col_func <- colorspace::diverge_hcl(25)
col_breaks = c(-15,-5,-2,-1.1, # for low
               seq(-1,1,length=21), # for blue
               1.1,2,5,15) # for high
test <- colorspace::diverge_hcl(28)

hv <- heatmap.2(as.matrix(tvec), 
                trace="none",          
                col = some_col_func)
names(hv)
#color mapping
crang <- hv$colorTable
crang$color <- as.character(crang$color)
# Extract the range associated with white
hv$colorTable[hv$colorTable[,"color"]=="#E2E2E2",]

png(paste0(figp,"hclust_wardD_heatmap_LME_A_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(tvec), 
                  main = "All fish",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
                  #RowSideColors = ccol$col,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  #add.expr=text(x=0.2, y=-0.1, srt=45, xpd=NA, adj=0, labels="Mag"),
                  col = some_col_func #, breaks=col_breaks
)
dev.off()


## Forage
Fcoef <- fcoef[,2:9]
rownames(Fcoef) <- fcoef$LME
Fcoef[is.na(Fcoef)] <- 0
d_fish <- dist(Fcoef) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish, k=8) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_F_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Fcoef), 
                  main = "Forage fish",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()


## Lg Pel
Pcoef <- pcoef[,2:9]
rownames(Pcoef) <- pcoef$LME
Pcoef[is.na(Pcoef)] <- 0
d_fish <- dist(Pcoef) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish, k=8) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_P_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Pcoef), 
                  main = "Large Pelagics",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()


## Dem
Dcoef <- dcoef[,2:9]
rownames(Dcoef) <- dcoef$LME
Dcoef[is.na(Dcoef)] <- 0
d_fish <- dist(Dcoef) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish, k=8) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_D_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Dcoef), 
                  main = "Demersals",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()


## Bent
Bcoef <- bcoef[,2:9]
rownames(Bcoef) <- bcoef$LME
Bcoef[is.na(Bcoef)] <- 0
d_fish <- dist(Bcoef) # method="man" # is a bit better
hc_fish <- hclust(d_fish, method = "ward.D")
dfish <- as.dendrogram(hc_fish)
# Color the branches based on the clusters:
dfish <- color_branches(dfish, k=8) 
# reduce the size of the labels:
dfish <- set(dfish, "labels_cex", 0.8)

png(paste0(figp,"hclust_wardD_heatmap_LME_B_mlr_coeffs_ALLdiv2SD.png"),    # create PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels
    height = 6.5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)  
gplots::heatmap.2(as.matrix(Bcoef), 
                  main = "Benthos",
                  srtCol = 45,
                  dendrogram = "row",
                  Rowv = dfish,
                  Colv = "NA", # this to make sure the columns are not ordered
                  trace="none",          
                  margins =c(3,6),      
                  key.xlab = "Coefficient",
                  denscol = "grey",
                  key.title=NA, # no title
                  density.info = "none",
                  offsetCol=-0.6,
                  col = some_col_func)
dev.off()


