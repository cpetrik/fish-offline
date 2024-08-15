% Visualize output of FEISTY
% CESM FOSI
% Time series plots and maps

clear 
close all

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod = 'v15_All_fish03';
%mod = 'All_fish03';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

% load means
load([fpath 'Plot_Means_FOSI_' mod '_' cfile '.mat'],...
    'AllF','AllP','AllD');
All = AllF + AllP + AllD;

load([fpath 'Plot_Nu_Means_FOSI_' mod '_' cfile '.mat'],...
    'Pmf','Plp','Pld','PA');

% annual
Pmf = Pmf * 365;
Plp = Plp * 365;
Pld = Pld * 365;
PA = PA * 365;

%% CVs
load([fpath 'FEISTY_FOSI_',mod,'_interann_var.mat'],...
    'cvf','cvp','cvd','cva');
bcvf = cvf;
bcvp = cvp;
bcvd = cvd;
bcva = cva;

clear cvf cvp cvd cva

load([fpath 'FEISTY_FOSI_',mod,'_nu_interann_var.mat'],...
    'cvf','cvp','cvd','cva');
pcvf = cvf;
pcvp = cvp;
pcvd = cvd;
pcva = cva;

clear cvf cvp cvd cva

%% Map data
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
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

cmYOR=cbrewer('seq','YlOrRd',66,'pchip');
cmYOB=cbrewer('seq','YlOrBr',66,'pchip');
cmOR=cbrewer('seq','OrRd',66,'pchip');
cmBP50=cbrewer('seq','BuPu',50,'PCHIP');

%% fix seam
[TLAT2,TLON2,AllF2] = cyclic_map_seam_cesm(TLAT,TLONG,AllF);
[~,~,AllP2] = cyclic_map_seam_cesm(TLAT,TLONG,AllP);
[~,~,AllD2] = cyclic_map_seam_cesm(TLAT,TLONG,AllD);
[~,~,All2]  = cyclic_map_seam_cesm(TLAT,TLONG,All);

[~,~,pF2] = cyclic_map_seam_cesm(TLAT,TLONG,Pmf);
[~,~,pP2] = cyclic_map_seam_cesm(TLAT,TLONG,Plp);
[~,~,pD2] = cyclic_map_seam_cesm(TLAT,TLONG,Pld);
[~,~,pA2]  = cyclic_map_seam_cesm(TLAT,TLONG,PA);

[~,~,cvFB] = cyclic_map_seam_cesm(TLAT,TLONG,bcvf);
[~,~,cvPB] = cyclic_map_seam_cesm(TLAT,TLONG,bcvp);
[~,~,cvDB] = cyclic_map_seam_cesm(TLAT,TLONG,bcvd);
[~,~,cvAB] = cyclic_map_seam_cesm(TLAT,TLONG,bcva);

[~,~,cvFP] = cyclic_map_seam_cesm(TLAT,TLONG,pcvf);
[~,~,cvPP] = cyclic_map_seam_cesm(TLAT,TLONG,pcvp);
[~,~,cvDP] = cyclic_map_seam_cesm(TLAT,TLONG,pcvd);
[~,~,cvAP] = cyclic_map_seam_cesm(TLAT,TLONG,pcva);

%% Remove low biom areas
pF2(AllF2<1e-9) = nan;
pP2(AllP2<1e-6) = nan;
pD2(AllD2<1e-6) = nan;
pA2(All2<1e-6) = nan;

cvFB(AllF2<1e-9) = nan;
cvPB(AllP2<1e-6) = nan;
cvDB(AllD2<1e-6) = nan;
cvAB(All2<1e-6) = nan;

cvFP(AllF2<1e-9) = nan;
cvPP(AllP2<1e-6) = nan;
cvDP(AllD2<1e-6) = nan;
cvAP(All2<1e-6) = nan;

AllF2(AllF2<1e-9) = nan;
AllP2(AllP2<1e-6) = nan;
AllD2(AllD2<1e-6) = nan;
All2(All2<1e-6) = nan;

% pF2(pF2<0) = nan;
% pP2(pP2<0) = nan;
% pD2(pD2<0) = nan;
% pA2(pA2<0) = nan;

%% 8plot by type
f1 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - F biom
subplot('Position',[0.015 0.75 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllF2))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1.5 2.0])
set(gcf,'renderer','painters')
text(0,1.75,'F biom','HorizontalAlignment','center')

%B - P biom
subplot('Position',[0.015 0.5 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllP2))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1.5 2.0])
set(gcf,'renderer','painters')
text(0,1.75,'P biom','HorizontalAlignment','center')

%C - D biom
subplot('Position',[0.015 0.25 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(AllD2))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1.5 2.0])
set(gcf,'renderer','painters')
text(0,1.75,'D biom','HorizontalAlignment','center')
colorbar('Position',[0.445 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out',...
    'Ticks',-1:2,'TickLabels',-1:2)

%D - A biom
subplot('Position',[0.015 0.0 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(All2))
colormap(cmBP50)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1.5 2.0])
set(gcf,'renderer','painters')
text(0,1.75,'All biom','HorizontalAlignment','center')


%E - F prod
subplot('Position',[0.51 0.75 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(pF2+eps))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1])
%colorbar('Position',[0.95 0.765 0.025 0.225],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'F prod','HorizontalAlignment','center')

%F - P prod
subplot('Position',[0.51 0.5 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(pP2+eps))
colormap(cmBP50)
colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1])
colorbar('Position',[0.945 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out',...
    'Ticks',-1:1,'TickLabels',-1:1)
set(gcf,'renderer','painters')
text(0,1.75,'P prod','HorizontalAlignment','center')

%G - D prod
subplot('Position',[0.51 0.25 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(pD2+eps))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1])
set(gcf,'renderer','painters')
text(0,1.75,'D','HorizontalAlignment','center')

%H - All prod
subplot('Position',[0.51 0.0 0.43 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,log10(pA2+eps))
colormap(cmBP50)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([-1 1])
set(gcf,'renderer','painters')
text(0,1.75,'All prod','HorizontalAlignment','center')

print('-dpng',[ppath 'Map_FEISTY_FOSI_',mod,'_mean_biomass_prod_types.png'])


%% 8plot by type CoeffVar
f2 = figure('Units','inches','Position',[1 3 6.5 8]);
%f1.Units = 'inches';

%A - F
subplot('Position',[0.015 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,(cvFB))
colormap(cmYOR)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'F biom','HorizontalAlignment','center')

%B -P
subplot('Position',[0.015 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,(cvPB))
colormap(cmYOR)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'P biom','HorizontalAlignment','center')

%C - D
subplot('Position',[0.015 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,(cvDB))
colormap(cmYOR)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'D biom','HorizontalAlignment','center')

%D - A
subplot('Position',[0.015 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,(cvAB))
colormap(cmYOR)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'All biom','HorizontalAlignment','center')

%E - F
subplot('Position',[0.47 0.75 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,(cvFP))
colormap(cmYOR)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'F prod','HorizontalAlignment','center')

%F - P
subplot('Position',[0.47 0.5 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,(cvPP))
colormap(cmYOR)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1])
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
set(gcf,'renderer','painters')
text(0,1.75,'P prod','HorizontalAlignment','center')

%G - D
subplot('Position',[0.47 0.25 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,(cvDP))
colormap(cmYOR)
%colorbar
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'D prod','HorizontalAlignment','center')

%H - all
subplot('Position',[0.47 0.0 0.44 0.25])
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1,'origin',[0 -100 0])
surfm(TLAT2,TLON2,(cvAP))
colormap(cmYOR)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
clim([0 1])
set(gcf,'renderer','painters')
text(0,1.75,'All prod','HorizontalAlignment','center')

print('-dpng',[ppath 'Map_FEISTY_FOSI_',mod,'_CV_biom_prod_types.png'])
