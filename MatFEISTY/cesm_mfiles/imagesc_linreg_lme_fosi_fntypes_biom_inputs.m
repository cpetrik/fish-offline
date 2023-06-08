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
pmat = (squeeze(Cmat_CCE(:,:,1)));
[ni,nj] = size(pmat);

bmat = (squeeze(Cmat_CCE(:,:,2)));

dmat = (squeeze(Cmat_CCE(:,:,3)));

zmat = (squeeze(Cmat_CCE(:,:,4)));

lmat = (squeeze(Cmat_CCE(:,:,5)));



%pvals
psig = flipud(squeeze(Pmat_CCE(:,:,1))); %not sure if these should be flipped or not
[ni,nj] = size(psig);
psigy = ymat(psig < 0.05);
psigf = fmat(psig < 0.05);
%psig(0.1 > psig > 0.05) = "."

bsig = flipud(squeeze(Pmat_CCE(:,:,2)));
bsigy = ymat(bsig < 0.05);
bsigf = fmat(bsig < 0.05);

zsig = flipud(squeeze(Pmat_CCE(:,:,3)));
zsigy = ymat(zsig < 0.05);
zsigf = fmat(zsig < 0.05);

lsig = flipud(squeeze(Pmat_CCE(:,:,5)));
lsigy = ymat(lsig < 0.05);
lsigf = fmat(lsig < 0.05);

dsig = flipud(squeeze(Pmat_CCE(:,:,7)));
dsigy = ymat(dsig < 0.05);
dsigf = fmat(dsig < 0.05);

%
figure(1)
subplot(2,3,1)
imagesc(pmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1:6,'XTickLabel',0:5,...
    'YTick',1:8,'YTickLabel',fliplr(tfish))
title(tanom{1})
hold on;
%plot(psigy,psigf,'*');

subplot(2,3,2)
pcolor(zmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{3})
hold on;
%plot(zsigy,zsigf,'*');

subplot(2,3,3)
pcolor(lmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{5})
hold on;
%plot(lsigy,lsigf,'*');

subplot(2,3,4)
pcolor(bmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{2})
hold on;
%plot(bsigy,bsigf,'*');

subplot(2,3,5)
pcolor(dmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{7})
hold on;
%plot(dsigy,dsigf,'*');
c=colorbar('position',[0.65 0.1 0.025 0.35],'orientation','vertical');
c.Label.String = 'CCE linear regression coefficient';

print('-dpng',[ppath 'pcolor_CCE_FOSI_biom_regress_drivers_div2SD.png'])

%% AI
ymat = repmat([1.5:7.5],9,1);
fmat = repmat([1.5:9.5],7,1);
fmat = fmat';

%coeffs
pmat = flipud(squeeze(Cmat_AI(:,:,1)));
[ni,nj] = size(pmat);
pmat(:,nj+1) = nan;
pmat(ni+1,:) = nan;

bmat = flipud(squeeze(Cmat_AI(:,:,2)));
bmat(:,nj+1) = nan;
bmat(ni+1,:) = nan;

zmat = flipud(squeeze(Cmat_AI(:,:,3)));
zmat(:,nj+1) = nan;
zmat(ni+1,:) = nan;

lmat = flipud(squeeze(Cmat_AI(:,:,5)));
lmat(:,nj+1) = nan;
lmat(ni+1,:) = nan;

dmat = flipud(squeeze(Cmat_AI(:,:,7)));
dmat(:,nj+1) = nan;
dmat(ni+1,:) = nan;

%pvals
psig = flipud(squeeze(Pmat_AI(:,:,1)));
[ni,nj] = size(psig);
psigy = ymat(psig < 0.05);
psigf = fmat(psig < 0.05);
%psig(0.1 > psig > 0.05) = "."

bsig = flipud(squeeze(Pmat_AI(:,:,2)));
bsigy = ymat(bsig < 0.05);
bsigf = fmat(bsig < 0.05);

zsig = flipud(squeeze(Pmat_AI(:,:,3)));
zsigy = ymat(zsig < 0.05);
zsigf = fmat(zsig < 0.05);

lsig = flipud(squeeze(Pmat_AI(:,:,5)));
lsigy = ymat(lsig < 0.05);
lsigf = fmat(lsig < 0.05);

dsig = flipud(squeeze(Pmat_AI(:,:,7)));
dsigy = ymat(dsig < 0.05);
dsigf = fmat(dsig < 0.05);

%
figure(2)
subplot(2,3,1)
pcolor(pmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{1})
hold on;
%plot(psigy,psigf,'*');

subplot(2,3,2)
pcolor(zmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{3})
hold on;
%plot(zsigy,zsigf,'*');

subplot(2,3,3)
pcolor(lmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{5})
hold on;
%plot(lsigy,lsigf,'*');

subplot(2,3,4)
pcolor(bmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{2})
hold on;
%plot(bsigy,bsigf,'*');

subplot(2,3,5)
pcolor(dmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{7})
hold on;
%plot(dsigy,dsigf,'*');
c=colorbar('position',[0.65 0.1 0.025 0.35],'orientation','vertical');
c.Label.String = 'AI linear regression coefficient';

print('-dpng',[ppath 'pcolor_AI_FOSI_biom_regress_drivers_div2SD.png'])

%% CHK
ymat = repmat([1.5:7.5],9,1);
fmat = repmat([1.5:9.5],7,1);
fmat = fmat';

%coeffs
pmat = flipud(squeeze(Cmat_CHK(:,:,1)));
[ni,nj] = size(pmat);
pmat(:,nj+1) = nan;
pmat(ni+1,:) = nan;

bmat = flipud(squeeze(Cmat_CHK(:,:,2)));
bmat(:,nj+1) = nan;
bmat(ni+1,:) = nan;

zmat = flipud(squeeze(Cmat_CHK(:,:,3)));
zmat(:,nj+1) = nan;
zmat(ni+1,:) = nan;

lmat = flipud(squeeze(Cmat_CHK(:,:,5)));
lmat(:,nj+1) = nan;
lmat(ni+1,:) = nan;

dmat = flipud(squeeze(Cmat_CHK(:,:,7)));
dmat(:,nj+1) = nan;
dmat(ni+1,:) = nan;

%pvals
psig = flipud(squeeze(Pmat_CHK(:,:,1)));
[ni,nj] = size(psig);
psigy = ymat(psig < 0.05);
psigf = fmat(psig < 0.05);
%psig(0.1 > psig > 0.05) = "."

bsig = flipud(squeeze(Pmat_CHK(:,:,2)));
bsigy = ymat(bsig < 0.05);
bsigf = fmat(bsig < 0.05);

zsig = flipud(squeeze(Pmat_CHK(:,:,3)));
zsigy = ymat(zsig < 0.05);
zsigf = fmat(zsig < 0.05);

lsig = flipud(squeeze(Pmat_CHK(:,:,5)));
lsigy = ymat(lsig < 0.05);
lsigf = fmat(lsig < 0.05);

dsig = flipud(squeeze(Pmat_CHK(:,:,7)));
dsigy = ymat(dsig < 0.05);
dsigf = fmat(dsig < 0.05);

%
figure(3)
subplot(2,3,1)
pcolor(pmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{1})
hold on;
%plot(psigy,psigf,'*');

subplot(2,3,2)
pcolor(zmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{3})
hold on;
%plot(zsigy,zsigf,'*');

subplot(2,3,3)
pcolor(lmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{5})
hold on;
%plot(lsigy,lsigf,'*');

subplot(2,3,4)
pcolor(bmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{2})
hold on;
%plot(bsigy,bsigf,'*');

subplot(2,3,5)
pcolor(dmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{7})
hold on;
%plot(dsigy,dsigf,'*');
c=colorbar('position',[0.65 0.1 0.025 0.35],'orientation','vertical');
c.Label.String = 'CHK linear regression coefficient';

print('-dpng',[ppath 'pcolor_CHK_FOSI_biom_regress_drivers_div2SD.png'])

%% EBS
ymat = repmat([1.5:7.5],9,1);
fmat = repmat([1.5:9.5],7,1);
fmat = fmat';

%coeffs
pmat = flipud(squeeze(Cmat_EBS(:,:,1)));
[ni,nj] = size(pmat);
pmat(:,nj+1) = nan;
pmat(ni+1,:) = nan;

bmat = flipud(squeeze(Cmat_EBS(:,:,2)));
bmat(:,nj+1) = nan;
bmat(ni+1,:) = nan;

zmat = flipud(squeeze(Cmat_EBS(:,:,3)));
zmat(:,nj+1) = nan;
zmat(ni+1,:) = nan;

lmat = flipud(squeeze(Cmat_EBS(:,:,5)));
lmat(:,nj+1) = nan;
lmat(ni+1,:) = nan;

dmat = flipud(squeeze(Cmat_EBS(:,:,7)));
dmat(:,nj+1) = nan;
dmat(ni+1,:) = nan;

%pvals
psig = flipud(squeeze(Pmat_EBS(:,:,1)));
[ni,nj] = size(psig);
psigy = ymat(psig < 0.05);
psigf = fmat(psig < 0.05);
%psig(0.1 > psig > 0.05) = "."

bsig = flipud(squeeze(Pmat_EBS(:,:,2)));
bsigy = ymat(bsig < 0.05);
bsigf = fmat(bsig < 0.05);

zsig = flipud(squeeze(Pmat_EBS(:,:,3)));
zsigy = ymat(zsig < 0.05);
zsigf = fmat(zsig < 0.05);

lsig = flipud(squeeze(Pmat_EBS(:,:,5)));
lsigy = ymat(lsig < 0.05);
lsigf = fmat(lsig < 0.05);

dsig = flipud(squeeze(Pmat_EBS(:,:,7)));
dsigy = ymat(dsig < 0.05);
dsigf = fmat(dsig < 0.05);

%
figure(4)
subplot(2,3,1)
pcolor(pmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{1})
hold on;
%plot(psigy,psigf,'*');

subplot(2,3,2)
pcolor(zmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{3})
hold on;
%plot(zsigy,zsigf,'*');

subplot(2,3,3)
pcolor(lmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{5})
hold on;
%plot(lsigy,lsigf,'*');

subplot(2,3,4)
pcolor(bmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{2})
hold on;
%plot(bsigy,bsigf,'*');

subplot(2,3,5)
pcolor(dmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{7})
hold on;
%plot(dsigy,dsigf,'*');
c=colorbar('position',[0.65 0.1 0.025 0.35],'orientation','vertical');
c.Label.String = 'EBS linear regression coefficient';

print('-dpng',[ppath 'pcolor_EBS_FOSI_biom_regress_drivers_div2SD.png'])

%% GAK
ymat = repmat([1.5:7.5],9,1);
fmat = repmat([1.5:9.5],7,1);
fmat = fmat';

%coeffs
pmat = flipud(squeeze(Cmat_GAK(:,:,1)));
[ni,nj] = size(pmat);
pmat(:,nj+1) = nan;
pmat(ni+1,:) = nan;

bmat = flipud(squeeze(Cmat_GAK(:,:,2)));
bmat(:,nj+1) = nan;
bmat(ni+1,:) = nan;

zmat = flipud(squeeze(Cmat_GAK(:,:,3)));
zmat(:,nj+1) = nan;
zmat(ni+1,:) = nan;

lmat = flipud(squeeze(Cmat_GAK(:,:,5)));
lmat(:,nj+1) = nan;
lmat(ni+1,:) = nan;

dmat = flipud(squeeze(Cmat_GAK(:,:,7)));
dmat(:,nj+1) = nan;
dmat(ni+1,:) = nan;

%pvals
psig = flipud(squeeze(Pmat_GAK(:,:,1)));
[ni,nj] = size(psig);
psigy = ymat(psig < 0.05);
psigf = fmat(psig < 0.05);
%psig(0.1 > psig > 0.05) = "."

bsig = flipud(squeeze(Pmat_GAK(:,:,2)));
bsigy = ymat(bsig < 0.05);
bsigf = fmat(bsig < 0.05);

zsig = flipud(squeeze(Pmat_GAK(:,:,3)));
zsigy = ymat(zsig < 0.05);
zsigf = fmat(zsig < 0.05);

lsig = flipud(squeeze(Pmat_GAK(:,:,5)));
lsigy = ymat(lsig < 0.05);
lsigf = fmat(lsig < 0.05);

dsig = flipud(squeeze(Pmat_GAK(:,:,7)));
dsigy = ymat(dsig < 0.05);
dsigf = fmat(dsig < 0.05);

%
figure(5)
subplot(2,3,1)
pcolor(pmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{1})
hold on;
%plot(psigy,psigf,'*');

subplot(2,3,2)
pcolor(zmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{3})
hold on;
%plot(zsigy,zsigf,'*');

subplot(2,3,3)
pcolor(lmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{5})
hold on;
%plot(lsigy,lsigf,'*');

subplot(2,3,4)
pcolor(bmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{2})
hold on;
%plot(bsigy,bsigf,'*');

subplot(2,3,5)
pcolor(dmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{7})
hold on;
%plot(dsigy,dsigf,'*');
c=colorbar('position',[0.65 0.1 0.025 0.35],'orientation','vertical');
c.Label.String = 'GAK linear regression coefficient';

print('-dpng',[ppath 'pcolor_GAK_FOSI_biom_regress_drivers_div2SD.png'])

%% GMX
ymat = repmat([1.5:7.5],9,1);
fmat = repmat([1.5:9.5],7,1);
fmat = fmat';

%coeffs
pmat = flipud(squeeze(Cmat_GMX(:,:,1)));
[ni,nj] = size(pmat);
pmat(:,nj+1) = nan;
pmat(ni+1,:) = nan;

bmat = flipud(squeeze(Cmat_GMX(:,:,2)));
bmat(:,nj+1) = nan;
bmat(ni+1,:) = nan;

zmat = flipud(squeeze(Cmat_GMX(:,:,3)));
zmat(:,nj+1) = nan;
zmat(ni+1,:) = nan;

lmat = flipud(squeeze(Cmat_GMX(:,:,5)));
lmat(:,nj+1) = nan;
lmat(ni+1,:) = nan;

dmat = flipud(squeeze(Cmat_GMX(:,:,7)));
dmat(:,nj+1) = nan;
dmat(ni+1,:) = nan;

%pvals
psig = flipud(squeeze(Pmat_GMX(:,:,1)));
[ni,nj] = size(psig);
psigy = ymat(psig < 0.05);
psigf = fmat(psig < 0.05);
%psig(0.1 > psig > 0.05) = "."

bsig = flipud(squeeze(Pmat_GMX(:,:,2)));
bsigy = ymat(bsig < 0.05);
bsigf = fmat(bsig < 0.05);

zsig = flipud(squeeze(Pmat_GMX(:,:,3)));
zsigy = ymat(zsig < 0.05);
zsigf = fmat(zsig < 0.05);

lsig = flipud(squeeze(Pmat_GMX(:,:,5)));
lsigy = ymat(lsig < 0.05);
lsigf = fmat(lsig < 0.05);

dsig = flipud(squeeze(Pmat_GMX(:,:,7)));
dsigy = ymat(dsig < 0.05);
dsigf = fmat(dsig < 0.05);

%
figure(6)
subplot(2,3,1)
pcolor(pmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{1})
hold on;
%plot(psigy,psigf,'*');

subplot(2,3,2)
pcolor(zmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{3})
hold on;
%plot(zsigy,zsigf,'*');

subplot(2,3,3)
pcolor(lmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{5})
hold on;
%plot(lsigy,lsigf,'*');

subplot(2,3,4)
pcolor(bmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{2})
hold on;
%plot(bsigy,bsigf,'*');

subplot(2,3,5)
pcolor(dmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{7})
hold on;
%plot(dsigy,dsigf,'*');
c=colorbar('position',[0.65 0.1 0.025 0.35],'orientation','vertical');
c.Label.String = 'GMX linear regression coefficient';

print('-dpng',[ppath 'pcolor_GMX_FOSI_biom_regress_drivers_div2SD.png'])

%% HI
ymat = repmat([1.5:7.5],9,1);
fmat = repmat([1.5:9.5],7,1);
fmat = fmat';

%coeffs
pmat = flipud(squeeze(Cmat_HI(:,:,1)));
[ni,nj] = size(pmat);
pmat(:,nj+1) = nan;
pmat(ni+1,:) = nan;

bmat = flipud(squeeze(Cmat_HI(:,:,2)));
bmat(:,nj+1) = nan;
bmat(ni+1,:) = nan;

zmat = flipud(squeeze(Cmat_HI(:,:,3)));
zmat(:,nj+1) = nan;
zmat(ni+1,:) = nan;

lmat = flipud(squeeze(Cmat_HI(:,:,5)));
lmat(:,nj+1) = nan;
lmat(ni+1,:) = nan;

dmat = flipud(squeeze(Cmat_HI(:,:,7)));
dmat(:,nj+1) = nan;
dmat(ni+1,:) = nan;

%pvals
psig = flipud(squeeze(Pmat_HI(:,:,1)));
[ni,nj] = size(psig);
psigy = ymat(psig < 0.05);
psigf = fmat(psig < 0.05);
%psig(0.1 > psig > 0.05) = "."

bsig = flipud(squeeze(Pmat_HI(:,:,2)));
bsigy = ymat(bsig < 0.05);
bsigf = fmat(bsig < 0.05);

zsig = flipud(squeeze(Pmat_HI(:,:,3)));
zsigy = ymat(zsig < 0.05);
zsigf = fmat(zsig < 0.05);

lsig = flipud(squeeze(Pmat_HI(:,:,5)));
lsigy = ymat(lsig < 0.05);
lsigf = fmat(lsig < 0.05);

dsig = flipud(squeeze(Pmat_HI(:,:,7)));
dsigy = ymat(dsig < 0.05);
dsigf = fmat(dsig < 0.05);

%
figure(7)
subplot(2,3,1)
pcolor(pmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{1})
hold on;
%plot(psigy,psigf,'*');

subplot(2,3,2)
pcolor(zmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{3})
hold on;
%plot(zsigy,zsigf,'*');

subplot(2,3,3)
pcolor(lmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{5})
hold on;
%plot(lsigy,lsigf,'*');

subplot(2,3,4)
pcolor(bmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{2})
hold on;
%plot(bsigy,bsigf,'*');

subplot(2,3,5)
pcolor(dmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{7})
hold on;
%plot(dsigy,dsigf,'*');
c=colorbar('position',[0.65 0.1 0.025 0.35],'orientation','vertical');
c.Label.String = 'HI linear regression coefficient';

print('-dpng',[ppath 'pcolor_HI_FOSI_biom_regress_drivers_div2SD.png'])

%% NE
ymat = repmat([1.5:7.5],9,1);
fmat = repmat([1.5:9.5],7,1);
fmat = fmat';

%coeffs
pmat = flipud(squeeze(Cmat_NE(:,:,1)));
[ni,nj] = size(pmat);
pmat(:,nj+1) = nan;
pmat(ni+1,:) = nan;

bmat = flipud(squeeze(Cmat_NE(:,:,2)));
bmat(:,nj+1) = nan;
bmat(ni+1,:) = nan;

zmat = flipud(squeeze(Cmat_NE(:,:,3)));
zmat(:,nj+1) = nan;
zmat(ni+1,:) = nan;

lmat = flipud(squeeze(Cmat_NE(:,:,5)));
lmat(:,nj+1) = nan;
lmat(ni+1,:) = nan;

dmat = flipud(squeeze(Cmat_NE(:,:,7)));
dmat(:,nj+1) = nan;
dmat(ni+1,:) = nan;

%pvals
psig = flipud(squeeze(Pmat_NE(:,:,1)));
[ni,nj] = size(psig);
psigy = ymat(psig < 0.05);
psigf = fmat(psig < 0.05);
%psig(0.1 > psig > 0.05) = "."

bsig = flipud(squeeze(Pmat_NE(:,:,2)));
bsigy = ymat(bsig < 0.05);
bsigf = fmat(bsig < 0.05);

zsig = flipud(squeeze(Pmat_NE(:,:,3)));
zsigy = ymat(zsig < 0.05);
zsigf = fmat(zsig < 0.05);

lsig = flipud(squeeze(Pmat_NE(:,:,5)));
lsigy = ymat(lsig < 0.05);
lsigf = fmat(lsig < 0.05);

dsig = flipud(squeeze(Pmat_NE(:,:,7)));
dsigy = ymat(dsig < 0.05);
dsigf = fmat(dsig < 0.05);

%
figure(8)
subplot(2,3,1)
pcolor(pmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{1})
hold on;
%plot(psigy,psigf,'*');

subplot(2,3,2)
pcolor(zmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{3})
hold on;
%plot(zsigy,zsigf,'*');

subplot(2,3,3)
pcolor(lmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{5})
hold on;
%plot(lsigy,lsigf,'*');

subplot(2,3,4)
pcolor(bmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{2})
hold on;
%plot(bsigy,bsigf,'*');

subplot(2,3,5)
pcolor(dmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{7})
hold on;
%plot(dsigy,dsigf,'*');
c=colorbar('position',[0.65 0.1 0.025 0.35],'orientation','vertical');
c.Label.String = 'NE linear regression coefficient';

print('-dpng',[ppath 'pcolor_NE_FOSI_biom_regress_drivers_div2SD.png'])


%% SE
ymat = repmat([1.5:7.5],9,1);
fmat = repmat([1.5:9.5],7,1);
fmat = fmat';

%coeffs
pmat = flipud(squeeze(Cmat_SE(:,:,1)));
[ni,nj] = size(pmat);
pmat(:,nj+1) = nan;
pmat(ni+1,:) = nan;

bmat = flipud(squeeze(Cmat_SE(:,:,2)));
bmat(:,nj+1) = nan;
bmat(ni+1,:) = nan;

zmat = flipud(squeeze(Cmat_SE(:,:,3)));
zmat(:,nj+1) = nan;
zmat(ni+1,:) = nan;

lmat = flipud(squeeze(Cmat_SE(:,:,5)));
lmat(:,nj+1) = nan;
lmat(ni+1,:) = nan;

dmat = flipud(squeeze(Cmat_SE(:,:,7)));
dmat(:,nj+1) = nan;
dmat(ni+1,:) = nan;

%pvals
psig = flipud(squeeze(Pmat_SE(:,:,1)));
[ni,nj] = size(psig);
psigy = ymat(psig < 0.05);
psigf = fmat(psig < 0.05);
%psig(0.1 > psig > 0.05) = "."

bsig = flipud(squeeze(Pmat_SE(:,:,2)));
bsigy = ymat(bsig < 0.05);
bsigf = fmat(bsig < 0.05);

zsig = flipud(squeeze(Pmat_SE(:,:,3)));
zsigy = ymat(zsig < 0.05);
zsigf = fmat(zsig < 0.05);

lsig = flipud(squeeze(Pmat_SE(:,:,5)));
lsigy = ymat(lsig < 0.05);
lsigf = fmat(lsig < 0.05);

dsig = flipud(squeeze(Pmat_SE(:,:,7)));
dsigy = ymat(dsig < 0.05);
dsigf = fmat(dsig < 0.05);

%
figure(9)
subplot(2,3,1)
pcolor(pmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{1})
hold on;
%plot(psigy,psigf,'*');

subplot(2,3,2)
pcolor(zmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{3})
hold on;
%plot(zsigy,zsigf,'*');

subplot(2,3,3)
pcolor(lmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{5})
hold on;
%plot(lsigy,lsigf,'*');

subplot(2,3,4)
pcolor(bmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{2})
hold on;
%plot(bsigy,bsigf,'*');

subplot(2,3,5)
pcolor(dmat); shading flat;
cmocean('balance')
caxis([-1 1])
xlabel('Lag')
ylabel('Type')
set(gca,'XTick',1.5:6.5,'XTickLabel',0:5,...
    'YTick',1.5:8.5,'YTickLabel',fliplr(tfish))
title(tanom{7})
hold on;
%plot(dsigy,dsigf,'*');
c=colorbar('position',[0.65 0.1 0.025 0.35],'orientation','vertical');
c.Label.String = 'SE linear regression coefficient';

print('-dpng',[ppath 'pcolor_SE_FOSI_biom_regress_drivers_div2SD.png'])


