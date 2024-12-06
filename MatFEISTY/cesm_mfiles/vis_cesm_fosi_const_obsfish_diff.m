% Visualize output of FEISTY
% CESM FOSI
% Time series plots and maps
% Difference between constant and observed fishing

clear 
close all

%% Map data
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);

[ni,nj]=size(TLONG);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

% colors
cm10=[0.5 0.5 0;... %tan/army
    0 0.7 0;...   %g
    1 0 1;...     %m
    1 0 0;...     %r
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    0 0.5 0.75;...   %med blue
    0 0 0.75;...    %b
    0.5 0.5 0.5; ...    %med grey
    0 0 0];...      %black

set(groot,'defaultAxesColorOrder',cm10);
cmBP50=cbrewer('seq','BuPu',50,'PCHIP');

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod1 = 'v15_All_fish03';
mod2 = 'v15_obsfish';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%% Const fishing
load([fpath 'Time_Means_FOSI_' mod1 '_' cfile '.mat'], 'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean',...
    'mf_ttc','mp_ttc','md_ttc','lp_ttc','ld_ttc');
load([fpath 'Space_Means_FOSI_' mod1 '_' cfile '.mat'], 'sf_sbio','sp_sbio','sd_sbio',...
    'mf_sbio','mp_sbio','md_sbio',...
    'lp_sbio','ld_sbio',...
    'mf_stc','mp_stc','md_stc','lp_stc','ld_stc');

% Types together
CtFb = sf_tmean + mf_tmean;
CtPb = sp_tmean + mp_tmean + lp_tmean;
CtDb = sd_tmean + md_tmean + ld_tmean;

CtFc = mf_ttc;
CtPc = mp_ttc + lp_ttc;
CtDc = md_ttc + ld_ttc;

CsFb = sf_sbio + mf_sbio;
CsPb = sp_sbio + mp_sbio + lp_sbio;
CsDb = sd_sbio + md_sbio + ld_sbio;

CsFc = mf_stc;
CsPc = mp_stc + lp_stc;
CsDc = md_stc + ld_stc;

%%
clear sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean lp_tmean 
clear ld_tmean mf_ttc mp_ttc md_ttc lp_ttc ld_ttc sf_sbio sp_sbio sd_sbio
clear mf_sbio mp_sbio md_sbio lp_sbio ld_sbio mf_stc mp_stc md_stc lp_stc ld_stc

%% Obs fishing
load([fpath 'Time_Means_FOSI_' mod2 '_' cfile '.mat'], 'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean',...
    'mf_ttc','mp_ttc','md_ttc','lp_ttc','ld_ttc');
load([fpath 'Space_Means_FOSI_' mod2 '_' cfile '.mat'], 'sf_sbio','sp_sbio','sd_sbio',...
    'mf_sbio','mp_sbio','md_sbio',...
    'lp_sbio','ld_sbio',...
    'mf_stc','mp_stc','md_stc','lp_stc','ld_stc');

% Types together
OtFb = sf_tmean + mf_tmean;
OtPb = sp_tmean + mp_tmean + lp_tmean;
OtDb = sd_tmean + md_tmean + ld_tmean;

OtFc = mf_ttc;
OtPc = mp_ttc + lp_ttc;
OtDc = md_ttc + ld_ttc;

OsFb = sf_sbio + mf_sbio;
OsPb = sp_sbio + mp_sbio + lp_sbio;
OsDb = sd_sbio + md_sbio + ld_sbio;

OsFc = mf_stc;
OsPc = mp_stc + lp_stc;
OsDc = md_stc + ld_stc;

%%
clear sf_tmean sp_tmean sd_tmean mf_tmean mp_tmean md_tmean lp_tmean 
clear ld_tmean mf_ttc mp_ttc md_ttc lp_ttc ld_ttc sf_sbio sp_sbio sd_sbio
clear mf_sbio mp_sbio md_sbio lp_sbio ld_sbio mf_stc mp_stc md_stc lp_stc ld_stc

%% Plots in time
t1 = 1:length(CtFb); %time; 1948 to 2015
t2 = 1:length(OtFb);
y1 = t1/12;
y2 = t2/12;

% All types
figure(1)
plot(y1,log10(CtFb),'r','Linewidth',2); hold on;
plot(y1,log10(CtPb),'b','Linewidth',2); hold on;
plot(y1,log10(CtDb),'k','Linewidth',2); hold on;
plot(y2,log10(OtFb),'--r','Linewidth',2); hold on;
plot(y2,log10(OtPb),'--b','Linewidth',2); hold on;
plot(y2,log10(OtDb),'--k','Linewidth',2); hold on;
legend('Fc','Pc','Dc','Fo','Po','Do')
legend('location','west')
xlim([y(1) y(end)])
%ylim([0.05 0.5])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('FOSI')
stamp('')
print('-dpng',[ppath 'FOSI_comp_biomass_all_types.png'])

figure(2)
plot(y1,log10(CtFc),'r','Linewidth',2); hold on;
plot(y1,log10(CtPc),'b','Linewidth',2); hold on;
plot(y1,log10(CtDc),'k','Linewidth',2); hold on;
plot(y2,log10(OtFc),'--r','Linewidth',2); hold on;
plot(y2,log10(OtPc),'--b','Linewidth',2); hold on;
plot(y2,log10(OtDc),'--k','Linewidth',2); hold on;
legend('Fc','Pc','Dc','Fo','Po','Do')
legend('location','west')
xlim([y(1) y(end)])
%ylim([0.05 0.5])
xlabel('Time (y)')
ylabel('log_1_0 Catch (g)')
title('FOSI')
stamp('')
print('-dpng',[ppath 'FOSI_comp_catch_all_types.png'])

%% Plots in space

ZCfb=NaN*ones(ni,nj);
ZCpb=NaN*ones(ni,nj);
ZCdb=NaN*ones(ni,nj);

ZCfb(GRD.ID)=CsFb;
ZCpb(GRD.ID)=CsPb;
ZCdb(GRD.ID)=CsDb;

ZCfc=NaN*ones(ni,nj);
ZCpc=NaN*ones(ni,nj);
ZCdc=NaN*ones(ni,nj);

ZCfc(GRD.ID)=CsFc;
ZCpc(GRD.ID)=CsPc;
ZCdc(GRD.ID)=CsDc;

ZOfb=NaN*ones(ni,nj);
ZOpb=NaN*ones(ni,nj);
ZOdb=NaN*ones(ni,nj);

ZOfb(GRD.ID)=OsFb;
ZOpb(GRD.ID)=OsPb;
ZOdb(GRD.ID)=OsDb;

ZOfc=NaN*ones(ni,nj);
ZOpc=NaN*ones(ni,nj);
ZOdc=NaN*ones(ni,nj);

ZOfc(GRD.ID)=OsFc;
ZOpc(GRD.ID)=OsPc;
ZOdc(GRD.ID)=OsDc;

%% save outputs for comparison
% save([fpath 'Plot_Means_FOSI_' mod '_' cfile '.mat'],'F','P','D','B',...
%     'AllF','AllP','AllD','AllS','AllM','AllL');

%% All 4 on subplots
figure(3)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(ZOfb - ZCfb))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-5 5]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('mean All F (g m^-^2)')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(ZOdb - ZCdb))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-5 5]);
set(gcf,'renderer','painters')
title('mean All D (g m^-^2)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,(ZOpb - ZCpb))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-5 5]);
set(gcf,'renderer','painters')
title('mean All P (g m^-^2)')

% All
% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(TLAT,TLONG,(ZOb - ZCb))
% cmocean('balance')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-5 5]);
% set(gcf,'renderer','painters')
% title(' mean All fishes (g m^-^2)')

stamp('')
print('-dpng',[ppath 'FOSI_comp_global_biomass_subplot.png'])

%% All 4 on subplots - Catch
figure(4)
% all F
subplot('Position',[0 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(ZOfc) - log10(ZCfc))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-5 5]);
colorbar('Position',[0.25 0.5 0.5 0.05],'orientation','horizontal')
set(gcf,'renderer','painters')
title('All F (g)')

% all D
subplot('Position',[0 0 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(ZOdc) - log10(ZCdc))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-5 5]);
set(gcf,'renderer','painters')
title('All D (g)')

% All P
subplot('Position',[0.5 0.51 0.5 0.5])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,log10(ZOpc) - log10(ZCpc))
cmocean('balance')
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-5 5]);
set(gcf,'renderer','painters')
title('All P (g)')

% All
% subplot('Position',[0.5 0 0.5 0.5])
% axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
%     'Grid','off','FLineWidth',1)
% surfm(TLAT,TLONG,log10(ZOc) - log10(ZCc))
% cmocean('balance')
% h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
% clim([-5 5]);
% set(gcf,'renderer','painters')
% title(' mean All fishes (g)')

stamp('')
print('-dpng',[ppath 'FOSI_comp_global_catch_subplot.png'])





