% Pcolor plots instead of R heatmap
% fish biomass-driver linear regression coeffs for color
% p-val as star or dot

clear 
close all

%% Fish data
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

load([fpath,'FOSI_biom_regress_drivers_div2SD.mat']);

ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/'...
    cfile '/corrs/'];

% LMEs
% lid = [54,1:2,65,10,3,5:7]; %ADD 65 = Aleutian Islands
% lname = {'CHK','EBS','GAK','AI','HI','CCE','GMX','SE','NE'};

%% CCE
ymat = repmat([1.5:7.5],9,1);
fmat = repmat([1.5:9.5],7,1);
fmat = fmat';

%coeffs
cpmat = (squeeze(Cmat_CCE(:,:,1)));
[ni,nj] = size(cpmat);
cpmat(:,nj+1) = nan;
cpmat(ni+1,:) = nan;

cbmat = (squeeze(Cmat_CCE(:,:,2)));
cbmat(:,nj+1) = nan;
cbmat(ni+1,:) = nan;

cdmat = (squeeze(Cmat_CCE(:,:,3)));
cdmat(:,nj+1) = nan;
cdmat(ni+1,:) = nan;

czmat = (squeeze(Cmat_CCE(:,:,4)));
czmat(:,nj+1) = nan;
czmat(ni+1,:) = nan;

clmat = (squeeze(Cmat_CCE(:,:,5)));
clmat(:,nj+1) = nan;
clmat(ni+1,:) = nan;

%pvals
cpsig = flipud(squeeze(Pmat_CCE(:,:,1))); %not sure if these should be flipped or not
cpsigy = ymat(cpsig < 0.05);
cpsigf = fmat(cpsig < 0.05);
%psig(0.1 > psig > 0.05) = "."

cbsig = flipud(squeeze(Pmat_CCE(:,:,2)));
cbsigy = ymat(cbsig < 0.05);
cbsigf = fmat(cbsig < 0.05);

czsig = flipud(squeeze(Pmat_CCE(:,:,4)));
czsigy = ymat(czsig < 0.05);
czsigf = fmat(czsig < 0.05);

clsig = flipud(squeeze(Pmat_CCE(:,:,5)));
clsigy = ymat(clsig < 0.05);
clsigf = fmat(clsig < 0.05);

cdsig = flipud(squeeze(Pmat_CCE(:,:,3)));
cdsigy = ymat(cdsig < 0.05);
cdsigf = fmat(cdsig < 0.05);

%% AI

%coeffs
apmat = (squeeze(Cmat_AI(:,:,1)));
[ni,nj] = size(apmat);
apmat(:,nj+1) = nan;
apmat(ni+1,:) = nan;

abmat = (squeeze(Cmat_AI(:,:,2)));
abmat(:,nj+1) = nan;
abmat(ni+1,:) = nan;

azmat = (squeeze(Cmat_AI(:,:,4)));
azmat(:,nj+1) = nan;
azmat(ni+1,:) = nan;

almat = (squeeze(Cmat_AI(:,:,5)));
almat(:,nj+1) = nan;
almat(ni+1,:) = nan;

admat = (squeeze(Cmat_AI(:,:,3)));
admat(:,nj+1) = nan;
admat(ni+1,:) = nan;

%pvals
apsig = flipud(squeeze(Pmat_AI(:,:,1)));
apsigy = ymat(apsig < 0.05);
apsigf = fmat(apsig < 0.05);
%psig(0.1 > psig > 0.05) = "."

absig = flipud(squeeze(Pmat_AI(:,:,2)));
absigy = ymat(absig < 0.05);
absigf = fmat(absig < 0.05);

azsig = flipud(squeeze(Pmat_AI(:,:,4)));
azsigy = ymat(azsig < 0.05);
azsigf = fmat(azsig < 0.05);

alsig = flipud(squeeze(Pmat_AI(:,:,5)));
alsigy = ymat(alsig < 0.05);
alsigf = fmat(alsig < 0.05);

adsig = flipud(squeeze(Pmat_AI(:,:,3)));
adsigy = ymat(adsig < 0.05);
adsigf = fmat(adsig < 0.05);


%% CHK

%coeffs
kpmat = (squeeze(Cmat_CHK(:,:,1)));
[ni,nj] = size(kpmat);
kpmat(:,nj+1) = nan;
kpmat(ni+1,:) = nan;

kbmat = (squeeze(Cmat_CHK(:,:,2)));
kbmat(:,nj+1) = nan;
kbmat(ni+1,:) = nan;

kzmat = (squeeze(Cmat_CHK(:,:,4)));
kzmat(:,nj+1) = nan;
kzmat(ni+1,:) = nan;

klmat = (squeeze(Cmat_CHK(:,:,5)));
klmat(:,nj+1) = nan;
klmat(ni+1,:) = nan;

kdmat = (squeeze(Cmat_CHK(:,:,3)));
kdmat(:,nj+1) = nan;
kdmat(ni+1,:) = nan;

%pvals
kpsig = flipud(squeeze(Pmat_CHK(:,:,1)));
kpsigy = ymat(kpsig < 0.05);
kpsigf = fmat(kpsig < 0.05);
%psig(0.1 > psig > 0.05) = "."

kbsig = flipud(squeeze(Pmat_CHK(:,:,2)));
kbsigy = ymat(kbsig < 0.05);
kbsigf = fmat(kbsig < 0.05);

kzsig = flipud(squeeze(Pmat_CHK(:,:,4)));
kzsigy = ymat(kzsig < 0.05);
kzsigf = fmat(kzsig < 0.05);

klsig = flipud(squeeze(Pmat_CHK(:,:,5)));
klsigy = ymat(klsig < 0.05);
klsigf = fmat(klsig < 0.05);

kdsig = flipud(squeeze(Pmat_CHK(:,:,3)));
kdsigy = ymat(kdsig < 0.05);
kdsigf = fmat(kdsig < 0.05);

%% EBS

%coeffs
epmat = (squeeze(Cmat_EBS(:,:,1)));
[ni,nj] = size(epmat);
epmat(:,nj+1) = nan;
epmat(ni+1,:) = nan;

ebmat = (squeeze(Cmat_EBS(:,:,2)));
ebmat(:,nj+1) = nan;
ebmat(ni+1,:) = nan;

ezmat = (squeeze(Cmat_EBS(:,:,4)));
ezmat(:,nj+1) = nan;
ezmat(ni+1,:) = nan;

elmat = (squeeze(Cmat_EBS(:,:,5)));
elmat(:,nj+1) = nan;
elmat(ni+1,:) = nan;

edmat = (squeeze(Cmat_EBS(:,:,3)));
edmat(:,nj+1) = nan;
edmat(ni+1,:) = nan;

%pvals
epsig = flipud(squeeze(Pmat_EBS(:,:,1)));
epsigy = ymat(epsig < 0.05);
epsigf = fmat(epsig < 0.05);
%psig(0.1 > psig > 0.05) = "."

ebsig = flipud(squeeze(Pmat_EBS(:,:,2)));
ebsigy = ymat(ebsig < 0.05);
ebsigf = fmat(ebsig < 0.05);

ezsig = flipud(squeeze(Pmat_EBS(:,:,4)));
ezsigy = ymat(ezsig < 0.05);
ezsigf = fmat(ezsig < 0.05);

elsig = flipud(squeeze(Pmat_EBS(:,:,5)));
elsigy = ymat(elsig < 0.05);
elsigf = fmat(elsig < 0.05);

edsig = flipud(squeeze(Pmat_EBS(:,:,3)));
edsigy = ymat(edsig < 0.05);
edsigf = fmat(edsig < 0.05);


%% GAK

%coeffs
gpmat = (squeeze(Cmat_GAK(:,:,1)));
[ni,nj] = size(gpmat);
gpmat(:,nj+1) = nan;
gpmat(ni+1,:) = nan;

gbmat = (squeeze(Cmat_GAK(:,:,2)));
gbmat(:,nj+1) = nan;
gbmat(ni+1,:) = nan;

gzmat = (squeeze(Cmat_GAK(:,:,4)));
gzmat(:,nj+1) = nan;
gzmat(ni+1,:) = nan;

glmat = (squeeze(Cmat_GAK(:,:,5)));
glmat(:,nj+1) = nan;
glmat(ni+1,:) = nan;

gdmat = (squeeze(Cmat_GAK(:,:,3)));
gdmat(:,nj+1) = nan;
gdmat(ni+1,:) = nan;

%pvals
gpsig = flipud(squeeze(Pmat_GAK(:,:,1)));
gpsigy = ymat(gpsig < 0.05);
gpsigf = fmat(gpsig < 0.05);
%psig(0.1 > psig > 0.05) = "."

gbsig = flipud(squeeze(Pmat_GAK(:,:,2)));
gbsigy = ymat(gbsig < 0.05);
gbsigf = fmat(gbsig < 0.05);

gzsig = flipud(squeeze(Pmat_GAK(:,:,4)));
gzsigy = ymat(gzsig < 0.05);
gzsigf = fmat(gzsig < 0.05);

glsig = flipud(squeeze(Pmat_GAK(:,:,5)));
glsigy = ymat(glsig < 0.05);
glsigf = fmat(glsig < 0.05);

gdsig = flipud(squeeze(Pmat_GAK(:,:,3)));
gdsigy = ymat(gdsig < 0.05);
gdsigf = fmat(gdsig < 0.05);

%% GMX

%coeffs
mpmat = (squeeze(Cmat_GMX(:,:,1)));
[ni,nj] = size(mpmat);
mpmat(:,nj+1) = nan;
mpmat(ni+1,:) = nan;

mbmat = (squeeze(Cmat_GMX(:,:,2)));
mbmat(:,nj+1) = nan;
mbmat(ni+1,:) = nan;

mzmat = (squeeze(Cmat_GMX(:,:,4)));
mzmat(:,nj+1) = nan;
mzmat(ni+1,:) = nan;

mlmat = (squeeze(Cmat_GMX(:,:,5)));
mlmat(:,nj+1) = nan;
mlmat(ni+1,:) = nan;

mdmat = (squeeze(Cmat_GMX(:,:,3)));
mdmat(:,nj+1) = nan;
mdmat(ni+1,:) = nan;

%pvals
mpsig = flipud(squeeze(Pmat_GMX(:,:,1)));
mpsigy = ymat(mpsig < 0.05);
mpsigf = fmat(mpsig < 0.05);
%psig(0.1 > psig > 0.05) = "."

mbsig = flipud(squeeze(Pmat_GMX(:,:,2)));
mbsigy = ymat(mbsig < 0.05);
mbsigf = fmat(mbsig < 0.05);

mzsig = flipud(squeeze(Pmat_GMX(:,:,4)));
mzsigy = ymat(mzsig < 0.05);
mzsigf = fmat(mzsig < 0.05);

mlsig = flipud(squeeze(Pmat_GMX(:,:,5)));
mlsigy = ymat(mlsig < 0.05);
mlsigf = fmat(mlsig < 0.05);

mdsig = flipud(squeeze(Pmat_GMX(:,:,3)));
mdsigy = ymat(mdsig < 0.05);
mdsigf = fmat(mdsig < 0.05);

%% HI

%coeffs
hpmat = (squeeze(Cmat_HI(:,:,1)));
[ni,nj] = size(hpmat);
hpmat(:,nj+1) = nan;
hpmat(ni+1,:) = nan;

hbmat = (squeeze(Cmat_HI(:,:,2)));
hbmat(:,nj+1) = nan;
hbmat(ni+1,:) = nan;

hzmat = (squeeze(Cmat_HI(:,:,4)));
hzmat(:,nj+1) = nan;
hzmat(ni+1,:) = nan;

hlmat = (squeeze(Cmat_HI(:,:,5)));
hlmat(:,nj+1) = nan;
hlmat(ni+1,:) = nan;

hdmat = (squeeze(Cmat_HI(:,:,3)));
hdmat(:,nj+1) = nan;
hdmat(ni+1,:) = nan;

%pvals
hpsig = flipud(squeeze(Pmat_HI(:,:,1)));
hpsigy = ymat(hpsig < 0.05);
hpsigf = fmat(hpsig < 0.05);
%psig(0.1 > psig > 0.05) = "."

hbsig = flipud(squeeze(Pmat_HI(:,:,2)));
hbsigy = ymat(hbsig < 0.05);
hbsigf = fmat(hbsig < 0.05);

hzsig = flipud(squeeze(Pmat_HI(:,:,4)));
hzsigy = ymat(hzsig < 0.05);
hzsigf = fmat(hzsig < 0.05);

hlsig = flipud(squeeze(Pmat_HI(:,:,5)));
hlsigy = ymat(hlsig < 0.05);
hlsigf = fmat(hlsig < 0.05);

hdsig = flipud(squeeze(Pmat_HI(:,:,3)));
hdsigy = ymat(hdsig < 0.05);
hdsigf = fmat(hdsig < 0.05);


%% NE

%coeffs
npmat = (squeeze(Cmat_NE(:,:,1)));
[ni,nj] = size(npmat);
npmat(:,nj+1) = nan;
npmat(ni+1,:) = nan;

nbmat = (squeeze(Cmat_NE(:,:,2)));
nbmat(:,nj+1) = nan;
nbmat(ni+1,:) = nan;

nzmat = (squeeze(Cmat_NE(:,:,4)));
nzmat(:,nj+1) = nan;
nzmat(ni+1,:) = nan;

nlmat = (squeeze(Cmat_NE(:,:,5)));
nlmat(:,nj+1) = nan;
nlmat(ni+1,:) = nan;

ndmat = (squeeze(Cmat_NE(:,:,3)));
ndmat(:,nj+1) = nan;
ndmat(ni+1,:) = nan;

%pvals
npsig = flipud(squeeze(Pmat_NE(:,:,1)));
npsigy = ymat(npsig < 0.05);
npsigf = fmat(npsig < 0.05);
%psig(0.1 > psig > 0.05) = "."

nbsig = flipud(squeeze(Pmat_NE(:,:,2)));
nbsigy = ymat(nbsig < 0.05);
nbsigf = fmat(nbsig < 0.05);

nzsig = flipud(squeeze(Pmat_NE(:,:,4)));
nzsigy = ymat(nzsig < 0.05);
nzsigf = fmat(nzsig < 0.05);

nlsig = flipud(squeeze(Pmat_NE(:,:,5)));
nlsigy = ymat(nlsig < 0.05);
nlsigf = fmat(nlsig < 0.05);

ndsig = flipud(squeeze(Pmat_NE(:,:,3)));
ndsigy = ymat(ndsig < 0.05);
ndsigf = fmat(ndsig < 0.05);


%% SE

%coeffs
spmat = (squeeze(Cmat_SE(:,:,1)));
[ni,nj] = size(spmat);
spmat(:,nj+1) = nan;
spmat(ni+1,:) = nan;

sbmat = (squeeze(Cmat_SE(:,:,2)));
sbmat(:,nj+1) = nan;
sbmat(ni+1,:) = nan;

szmat = (squeeze(Cmat_SE(:,:,4)));
szmat(:,nj+1) = nan;
szmat(ni+1,:) = nan;

slmat = (squeeze(Cmat_SE(:,:,5)));
slmat(:,nj+1) = nan;
slmat(ni+1,:) = nan;

sdmat = (squeeze(Cmat_SE(:,:,3)));
sdmat(:,nj+1) = nan;
sdmat(ni+1,:) = nan;

%pvals
spsig = flipud(squeeze(Pmat_SE(:,:,1)));
spsigy = ymat(spsig < 0.05);
spsigf = fmat(spsig < 0.05);
%psig(0.1 > psig > 0.05) = "."

sbsig = flipud(squeeze(Pmat_SE(:,:,2)));
sbsigy = ymat(sbsig < 0.05);
sbsigf = fmat(sbsig < 0.05);

szsig = flipud(squeeze(Pmat_SE(:,:,4)));
szsigy = ymat(szsig < 0.05);
szsigf = fmat(szsig < 0.05);

slsig = flipud(squeeze(Pmat_SE(:,:,5)));
slsigy = ymat(slsig < 0.05);
slsigf = fmat(slsig < 0.05);

sdsig = flipud(squeeze(Pmat_SE(:,:,3)));
sdsigy = ymat(sdsig < 0.05);
sdsigf = fmat(sdsig < 0.05);

%% Figures =======================================================
figure(1)
% CCE
subplot(4,5,1)
pcolor(cpmat); shading flat;
cmocean('balance')
caxis([-1 1])
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
title(tanom{1})
text(1,9,'CCE','FontWeight','bold')
%plot(cpsigy,cpsigf,'*');

subplot(4,5,2)
pcolor(czmat); shading flat;
cmocean('balance')
caxis([-1 1])
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
title(tanom{4})
hold on;
%plot(czsigy,czsigf,'*');

subplot(4,5,3)
pcolor(clmat); shading flat;
cmocean('balance')
caxis([-1 1])
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
title(tanom{5})
hold on;
%plot(clsigy,clsigf,'*');

subplot(4,5,4)
pcolor(cbmat); shading flat;
cmocean('balance')
caxis([-1 1])
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
title(tanom{2})
hold on;
%plot(cbsigy,cbsigf,'*');

subplot(4,5,5)
pcolor(cdmat); shading flat;
cmocean('balance')
caxis([-1 1])
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
title(tanom{3})
hold on;
%plot(cdsigy,cdsigf,'*');

% EBS
subplot(4,5,6)
pcolor(epmat); shading flat;
cmocean('balance')
caxis([-1 1])
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
text(1,9,'EBS','FontWeight','bold')
hold on;
%plot(psigy,psigf,'*');

subplot(4,5,7)
pcolor(ezmat); shading flat;
cmocean('balance')
caxis([-1 1])
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
%plot(zsigy,zsigf,'*');

subplot(4,5,8)
pcolor(elmat); shading flat;
cmocean('balance')
caxis([-1 1])
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
%plot(lsigy,lsigf,'*');

subplot(4,5,9)
pcolor(ebmat); shading flat;
cmocean('balance')
caxis([-1 1])
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
%plot(bsigy,bsigf,'*');

subplot(4,5,10)
pcolor(edmat); shading flat;
cmocean('balance')
caxis([-1 1])
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))

%CHK
subplot(4,5,11)
pcolor(kpmat); shading flat;
cmocean('balance')
caxis([-1 1])
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
text(1,9,'CHK','FontWeight','bold')
hold on;
%plot(psigy,psigf,'*');

subplot(4,5,12)
pcolor(kzmat); shading flat;
cmocean('balance')
caxis([-1 1])
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
%plot(zsigy,zsigf,'*');

subplot(4,5,13)
pcolor(klmat); shading flat;
cmocean('balance')
caxis([-1 1])
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
%plot(lsigy,lsigf,'*');

subplot(4,5,14)
pcolor(kbmat); shading flat;
cmocean('balance')
caxis([-1 1])
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
%plot(bsigy,bsigf,'*');

subplot(4,5,15)
pcolor(kdmat); shading flat;
cmocean('balance')
caxis([-1 1])
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))

%NE
subplot(4,5,16)
pcolor(npmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
text(1,9,'NEUS','FontWeight','bold')
hold on;
%plot(psigy,psigf,'*');

subplot(4,5,17)
pcolor(nzmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
hold on;
%plot(zsigy,zsigf,'*');

subplot(4,5,18)
pcolor(nlmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
hold on;
%plot(lsigy,lsigf,'*');

subplot(4,5,19)
pcolor(nbmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))
hold on;
%plot(bsigy,bsigf,'*');

subplot(4,5,20)
pcolor(ndmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',(tfish))

c=colorbar('position',[0.925 0.35 0.025 0.35],'orientation','vertical');
c.Label.String = 'CCE linear regression coefficient';

print('-dpng',[ppath 'pcolor_CCE_EBS_GAK_NE_coef_inputs_fntypes_biom.png'])

