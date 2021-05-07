% CESM 4P-4Z output 

clear all
close all

%% Paths

fpath='/Volumes/MIP/GCM_DATA/CESM/4P4Z/';
gpath='/Volumes/MIP/GCM_DATA/CESM/FOSI/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/4P4Z/';

load([fpath 'gridspec_POP_gx1v6_4p4z.mat']);
load([fpath 'Data_grid_POP_gx1v6_4p4z.mat'],'GRD');
load([fpath 'Data_cesm_4p4z_daily_1.mat']);

%%
c_Tp = mean(ESM.Tp,2);
c_Tb = mean(ESM.Tb,2);
c_D = mean(ESM.det,2);
c_Zm = mean(ESM.Zm,2);
c_Zl = mean(ESM.Zl,2);
c_dZm = mean(ESM.dZm,2);
c_dZl = mean(ESM.dZl,2);

[mi,mj]=size(TLONG);
cTp=NaN*ones(mi,mj);
cTb=NaN*ones(mi,mj);
cD=NaN*ones(mi,mj);
cZm=NaN*ones(mi,mj);
cZl=NaN*ones(mi,mj);
cdZm=NaN*ones(mi,mj);
cdZl=NaN*ones(mi,mj);

cTp(GRD.ID) = c_Tp;
cTb(GRD.ID) = c_Tb;
cD(GRD.ID) = c_D;
cZm(GRD.ID) = c_Zm;
cZl(GRD.ID) = c_Zl;
cdZm(GRD.ID) = c_dZm;
cdZl(GRD.ID) = c_dZl;

%%
clatlim=[-90 90];
clonlim=[-280 80];

figure(2)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cTp)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('Tp')
% load coast;                     
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,cTb)
cmocean('thermal')
caxis([0 35])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('Tb')

% subplot('Position',[0.5 0.51 0.5 0.5])
% axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(TLAT,TLONG,log10(cZ))
% cmocean('tempo')
% caxis([0 2])
% colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
% title('log_1_0 mesoZoo g m^-^2')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cD))
cmocean('tempo')
caxis([-0.5 1.5])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('log_1_0 Det g m^-^2 d^-^1')
print('-dpng',[pp 'Map_CESM_4P4Z_yr1_from_daily_interp_temp_det.png'])

%%
figure(3)
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZm))
cmocean('tempo')
caxis([0 2])
colorbar('Position',[0.05 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 Zoo3 g m^-^2')
% load coast;                     
% h=patchm(lat+0.5,long+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cdZm))
cmocean('tempo')
caxis([-1.5 0.5])
colorbar('Position',[0.05 0.05 0.4 0.03],'orientation','horizontal')
title('log_1_0 Zoo3 loss g m^-^2 d^-^1')

subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cZl))
cmocean('tempo')
caxis([0 2])
colorbar('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 Zoo4 g m^-^2')

subplot('Position',[0.5 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(cdZl))
cmocean('tempo')
caxis([-2 -1])
colorbar('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('log_1_0 Zoo4 loss g m^-^2 d^-^1')
print('-dpng',[pp 'Map_CESM_4P4Z_yr1_from_daily_interp_zoop.png'])



