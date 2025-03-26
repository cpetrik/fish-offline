% CESM DPLE ms figs

clear 
close all

%% Paths
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/CESM_MAPP/DPLE/';

% cpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';
% load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);

%% Fig 1 - TB
fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/DPLE_ms/DAIO_save_nc/';
fold1 = 'Fig.1/';
load([fpath fold1 'temp_bottom_pval.mat'],'predic','pval','diff','diffpval',...
    'lat','lon')
% predic, pval, diff, diffpval

[LAT,LON] = meshgrid(lat,lon);
clatlim=[-90 90];
clonlim=[-180 180];
load coastlines;                     

tb_acc = double(squeeze(predic(:,:,1)));
tb_pval = double(squeeze(pval(:,:,1)));
tb_diff = double(squeeze(diff(:,:,1)));
tb_dpval = double(squeeze(diffpval(:,:,1)));

tb_acc(tb_acc<-1e2)=nan;
tb_pval(tb_pval<-1e2)=nan;
tb_diff(tb_diff<-1e2)=nan;
tb_dpval(tb_dpval<-1e2)=nan;

clear predic pval diff diffpval

figure(1)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tb_acc)
cmocean('balance')
clim([-1 1])
colorbar %('orientation','horizontal')
title('TB LY 1-3')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% Fig 2 - Zloss
fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/DPLE_ms/DAIO_save_nc/';
fold2 = 'Fig.2/';
load([fpath fold2 'LzooC_150m_pval.mat'],'predic','pval','diff','diffpval')
% predic, pval, diff, diffpval

zl_acc = double(squeeze(predic(:,:,1)));
zl_pval = double(squeeze(pval(:,:,1)));
zl_diff = double(squeeze(diff(:,:,1)));
zl_dpval = double(squeeze(diffpval(:,:,1)));

zl_acc(zl_acc<-1e2)=nan;
zl_pval(zl_pval<-1e2)=nan;
zl_diff(zl_diff<-1e2)=nan;
zl_dpval(zl_dpval<-1e2)=nan;

clear predic pval diff diffpval

figure(2)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,zl_acc)
cmocean('balance')
clim([-1 1])
colorbar %('orientation','horizontal')
title('ZmLoss LY 1-3')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% Fig 4 - Fishes
fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/DPLE_ms/DAIO_save_nc/';
fold4 = 'Fig.4/';
load([fpath fold4 'tvb_pval.mat'],'predic','pval','diff','diffpval')
% predic, pval, diff, diffpval

tot_acc = double(squeeze(predic(:,:,1)));
tot_pval = double(squeeze(pval(:,:,1)));
tot_diff = double(squeeze(diff(:,:,1)));
tot_dpval = double(squeeze(diffpval(:,:,1)));

tot_acc(tot_acc<-1e2)=nan;
tot_pval(tot_pval<-1e2)=nan;
tot_diff(tot_diff<-1e2)=nan;
tot_dpval(tot_dpval<-1e2)=nan;

clear predic pval diff diffpval

figure(3)
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tot_acc)
cmocean('balance')
clim([-1 1])
colorbar %('orientation','horizontal')
title('Total fish biomass LY 1-3')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

%% only sig and greater than persist?
tb_acc_sig = tb_acc;
tb_acc_sig(tb_pval>0.05)=nan;
tb_acc_sig(tb_diff<=0)=nan;

zl_acc_sig = zl_acc;
zl_acc_sig(zl_pval>0.05)=nan;
zl_acc_sig(zl_diff<=0)=nan;

tot_acc_sig = tot_acc;
tot_acc_sig(tot_pval>0.05)=nan;
tot_acc_sig(tot_diff<=0)=nan;

%% NEUS maps
latlim=[20 60];
lonlim=[-100 -30];

figure(4)
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tb_acc_sig)
cmocean('balance')
clim([-1 1])
colorbar %('orientation','horizontal')
title('TB LY 1-3')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[pp 'NEUS_TB_DPLE_ACC.png'])

figure(5)
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,zl_acc_sig)
cmocean('balance')
clim([-1 1])
colorbar %('orientation','horizontal')
title('ZmLoss LY 1-3')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[pp 'NEUS_ZmLoss_DPLE_ACC.png'])

figure(6)
axesm ('Mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LON,tot_acc_sig)
cmocean('balance')
clim([-1 1])
colorbar %('orientation','horizontal')
title('Total fish biomass LY 1-3')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
print('-dpng',[pp 'NEUS_AllFish_DPLE_ACC.png'])


