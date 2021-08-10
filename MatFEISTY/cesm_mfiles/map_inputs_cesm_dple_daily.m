% CESM DPLE output 

clear all
close all

%% Paths

fpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
dpath='/Volumes/MIP/GCM_DATA/CESM/DPLE/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/DPLE/';

load([fpath 'gridspec_POP_gx1v6.mat']);
load([fpath 'Data_grid_POP_gx1v6.mat'],'GRD');
load([dpath 'Data_grid_POP_gx1v6_DPLE.mat']);
load([dpath 'Data_cesm_dple_daily_M1_Y1954_10.mat']);

%%
c_Tp = mean(ESM.Tp,2);
c_Tb = mean(ESM.Tb,2);
c_D = mean(ESM.det,2);
c_Z = mean(ESM.Zm,2);
c_Zl = mean(ESM.dZm,2);

[mi,mj]=size(TLONG);
cTp=NaN*ones(mi,mj);
cTb=NaN*ones(mi,mj);
cD=NaN*ones(mi,mj);
cZ=NaN*ones(mi,mj);
cZl=NaN*ones(mi,mj);

cTp(WID)=c_Tp;
cTb(WID)=c_Tb;
cD(WID) =c_D;
cZ(WID) =c_Z;
cZl(WID) =c_Zl;

%%
clatlim=[-90 90];
clonlim=[-180 180];

figure(2)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cTp)
cmocean('thermal')
%caxis([0 35])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('Tp')
% load coast;                     
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cTb)
cmocean('thermal')
%caxis([0 35])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Tb')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZ))
cmocean('tempo')
caxis([0 2])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 Zoo')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cD)
cmocean('tempo')
caxis([0 3])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('Det')
print('-dpng',[pp 'Map_CESM_DPLE_M1_Y1954_10_daily_interp_forcings.png'])

% WORKS NOW!!!

%%
tp = (nanmean(ESM.Tp));
tb = (nanmean(ESM.Tb));
zm = (nanmean(ESM.Zm))* 1e-9 * 1e4 * 12.01 * 9.0;
zl = (nanmean(ESM.dZm))* 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
det = (nanmean(ESM.det))* 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

figure(3)
subplot(2,3,1)
plot(1:365,squeeze(tp))
xlim([1 365])
title('Tp')

subplot(2,3,2)
plot(1:365,squeeze(tb))
xlim([1 365])
title('Tb')

subplot(2,3,3)
plot(1:365,squeeze(det))
xlim([1 365])
title('det')

subplot(2,3,4)
plot(1:365,squeeze(zm))
xlim([1 365])
title('zoo')

subplot(2,3,5)
plot(1:365,squeeze(zl))
xlim([1 365])
title('zoo loss')


