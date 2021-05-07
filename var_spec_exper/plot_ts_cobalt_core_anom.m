% Plot CORE seasonal anomalies

clear all
close all

spath='/Volumes/MIP/GCM_DATA/CORE-forced/';

% anomaly time series
load([spath 'cobalt_core_anom_1950_2007.mat'])

pp = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

%% time 
[nr,nc,nt] =size(det_anom);
time = yid;

%% reshape to save only land-free cells
% TP
nidxmat = reshape(1:nr*nc, nr, nc);
xTP = reshape(tp_anom, nr*nc, nt);
nnan = reshape(sum(isnan(xTP')), nr, nc);
nidx = nidxmat(nnan < 1);
xTP = xTP(nidx, :)';
% TB
xTB = reshape(tb_anom, nr*nc, nt);
xTB = xTB(nidx, :)';
% DET
xDet = reshape(det_anom, nr*nc, nt);
xDet = xDet(nidx, :)';
% ZOOMED_INT100
xMZ = reshape(mz_anom, nr*nc, nt);
xMZ = xMZ(nidx, :)';
% ZOOLRG_INT100
xLZ = reshape(lz_anom, nr*nc, nt);
xLZ = xLZ(nidx, :)';
% HPLOSS MED
xHPM = reshape(hpmz_anom, nr*nc, nt);
xHPM = xHPM(nidx, :)';
% HPLOSS LRG
xHPL = reshape(hplz_anom, nr*nc, nt);
xHPL = xHPL(nidx, :)';

%% time means
mtp = nanmean(xTP,2);
mtb = nanmean(xTB,2);
mdet = nanmean(xDet,2);
mzm = nanmean(xMZ,2);
mzl = nanmean(xLZ,2);
mhpmz = nanmean(xHPM,2);
mhplz = nanmean(xHPL,2);

%%
cm9=[0.5 0.5 0.5 ;...   %grey
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0 0.7 0]; %...    %green
    %0 0 0];...      %black
    
set(groot,'defaultAxesColorOrder',cm9);

%%
figure(1)
%subplot(2,2,1)
plot(time,mtp,'k','LineWidth',2);

figure
%subplot(2,2,3)
plot(time,mtb,'k','LineWidth',2);

figure
%subplot(2,2,4)
plot(time,mdet,'k','LineWidth',2);
%print('-dpng',[pp 'ts_CORE_temp_det_anom.png'])    

figure
%subplot(2,2,1)
plot(time,mzm,'k','LineWidth',2);

figure
%subplot(2,2,2)
plot(time,mhpmz,'k','LineWidth',2);

figure
%subplot(2,2,3)
plot(time,mzl,'k','LineWidth',2);

figure
%subplot(2,2,4)
plot(time,mhplz,'k','LineWidth',2);
%print('-dpng',[pp 'ts_CORE_zoo_anom.png'])    

%%
figure(1)
%subplot(2,2,1)
plot(time,xTP(:,100),'k','LineWidth',2);

figure
%subplot(2,2,3)
plot(time,xTB(:,100),'k','LineWidth',2);

figure
%subplot(2,2,4)
plot(time,xDet(:,100),'k','LineWidth',2);
%print('-dpng',[pp 'ts_CORE_temp_det_anom.png'])    

figure
%subplot(2,2,1)
plot(time,xMZ(:,100),'k','LineWidth',2);

figure
%subplot(2,2,2)
plot(time,xHPM(:,100),'k','LineWidth',2);

figure
%subplot(2,2,3)
plot(time,xLZ(:,100),'k','LineWidth',2);

figure
%subplot(2,2,4)
plot(time,xHPL(:,100),'k','LineWidth',2);