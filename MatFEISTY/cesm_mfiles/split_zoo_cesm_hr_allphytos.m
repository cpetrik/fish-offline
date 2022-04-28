% CESM High Res JRA forcing, Cal Curr only
% Use phyto size classes to split zoop into meso and micro

% there are explicit coccolithophores in this simulation 
% coccos as in the transition between small phytoplankton and diatoms - 
% some are quite large, some are small. They can be eaten by copepods 
% and also by microzooplankton. Perhaps you can count around 50-75% of
% cocco biomass going towards mesozooplankton partitioning?

clear all
close all

fpath='/Volumes/MIP/GCM_DATA/CESM_HR/';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/HR/';

%% Plankton inputs
load([fpath 'g.e22.G1850ECO_JRA_HR.TL319_t13.004.FIESTY-forcing_hiresJRA_CAcurr.mat'])

%%
% Size:       320x384x816
% Dimensions: nlon,nlat,time
%FillValue   = NaN;
%missing_value = 9.969209968386869e+36;

% units = 'mmol/m^3 cm';
% units = 'mmol/m^3/s cm';

%% doubles
zooC_100m = double(zooC_100m); 
zoo_loss_100m = double(zoo_loss_100m); 
diatC_100m = double(diatC_100m); 
diazC_100m = double(diazC_100m); 
coccoC_100m = double(coccoC_100m); 
spC_100m = double(spC_100m); 

%% fraction large phyto -> frac large zoo
fracL = (diatC_100m + 0.5*coccoC_100m) ./ ...
    (diatC_100m + spC_100m + diazC_100m + coccoC_100m + eps);

edges = [0:0.1:1];
histogram(fracL(:),edges)

%%
LzooC_100m = fracL .* zooC_100m;
Lzoo_loss_100m = fracL .* zoo_loss_100m;

SzooC_100m = (1.0-fracL) .* zooC_100m;

save([fpath 'g.e22.G1850ECO_JRA_HR.TL319_t13.004.FIESTY-forcing_hiresJRA_CAcurr_meszoo_totloss_allphytoC.mat'],...
    'fracL','SzooC_100m','LzooC_100m','Lzoo_loss_100m','zooC_100m_units',...
    'zoo_loss_100m_units');

%% Units
% meso zoo: nmolC cm-2 to g(WW) m-2
zooC_100m = zooC_100m * 1e-9 * 1e4 * 12.01 * 9.0;
LzooC_100m = LzooC_100m * 1e-9 * 1e4 * 12.01 * 9.0;
SzooC_100m = SzooC_100m * 1e-9 * 1e4 * 12.01 * 9.0;

% meso zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
zoo_loss_100m = zoo_loss_100m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
Lzoo_loss_100m = Lzoo_loss_100m * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

%%
cFL = mean(fracL,3);
cZ  = mean(LzooC_100m,3);
cZt = mean(Lzoo_loss_100m,3);
cT  = mean(zooC_100m,3);
cTt = mean(zoo_loss_100m,3);
cZS  = mean(SzooC_100m,3);

%% 
cZ(cZ<=0) = eps;
cT(cT<=0) = eps;
cTt(cTt<=0) = eps;
cZt(cZt<=0) = eps;
cZS(cZS<=0) = eps;

cZ(cZ>= 9.9e+36) = nan;
cT(cT>= 9.9e+36) = nan;
cTt(cTt>= 9.9e+36) = nan;
cZt(cZt>= 9.9e+36) = nan;
cZS(cZS>= 9.9e+36) = nan;

%%
% clatlim=[-90 90];
% clonlim=[-280 80];
plotminlat=28; %Set these bounds for your data
plotmaxlat=53;
plotminlon=220;
plotmaxlon=255;
clatlim=[plotminlat plotmaxlat];
clonlim=[plotminlon plotmaxlon];

load coastlines

LAT = squeeze(TLAT(:,:,1));
LONG = squeeze(TLONG(:,:,1));

%%
figure(1)
subplot('Position',[0.01 0.64 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LONG,log10(cT))
cmocean('tempo')
caxis([0.5 1.25])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 TotZoo')
%text(0.2,1.65,'log_1_0 TotZoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.64 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LONG,log10(cTt))
cmocean('tempo')
caxis([-0.5 0.5])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 TotZoo loss')
text(0.2,1.65,'log_1_0 TotZoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.325 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LONG,log10(cZS))
cmocean('tempo')
caxis([0 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 SmZoo')
%text(0.2,1.65,'log_1_0 SmZoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.325 0.4 0.3])
%subplot('Position',[0.21 0.31 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LONG,(cFL))
cmocean('balance')
caxis([0 1])
colorbar%('Position',[0.55 0.05 0.4 0.03],'orientation','horizontal')
title('Frac Lg')
%text(0.2,1.65,'Frac Lg','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.01 0.01 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LONG,log10(cZ))
cmocean('tempo')
caxis([0 1])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 LgZoo')
%text(0.2,1.65,'log_1_0 LgZoo','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

subplot('Position',[0.41 0.01 0.4 0.3])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(LAT,LONG,log10(cZt))
cmocean('tempo')
caxis([-0.5 0.5])
colorbar%('Position',[0.55 0.56 0.4 0.03],'orientation','horizontal')
title('log_1_0 LgZoo loss')
%text(0.2,1.65,'log_1_0 LgZoo loss','HorizontalAlignment','center','FontWeight','bold')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);

print('-dpng',[pp 'Map_CESM_JRA_HR_mean_zoo_loss_Bsplit_allphytos_50cocco.png'])




