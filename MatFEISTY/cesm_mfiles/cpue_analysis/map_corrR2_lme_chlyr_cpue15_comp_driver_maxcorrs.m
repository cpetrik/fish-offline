% Map corr R2s of ind driver combos with cpue
% Subplots together for comparison
% chl yrs only 1997-2015

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regress_cpue/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',...
    cfile,'/corrs_cpue/'];

mod = 'v15_All_fish03';
mod2 = 'v15_obsfish2015';

%%
load([spath,'LMEs_corr_cpue_chlyrs15_inputs_feisty_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab','lid')

AtabC = LAtab;
FtabC = LFtab;
PtabC = LPtab;
DtabC = LDtab;

clear LAtab LFtab LPtab LDtab

%%
load([spath,'LMEs_corr_cpue_chlyrs15_driver_obsfish2015_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab') 

AtabO = LAtab;
FtabO = LFtab;
PtabO = LPtab;
DtabO = LDtab;

clear LAtab LFtab LPtab LDtab

%%
load([spath,'LMEs_corr_cpue_chlyrs15_driver_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

AtabD = LAtab(lid,:);
FtabD = LFtab(lid,:);
PtabD = LPtab(lid,:);
DtabD = LDtab(lid,:);

clear LAtab LFtab LPtab LDtab

%%
load([spath,'LMEs_corr_cpue_chlyrs15_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

AtabS = LAtab(lid,:);
FtabS = LFtab(lid,:);
PtabS = LPtab(lid,:);
DtabS = LDtab(lid,:);

sst = (AtabS(:,4)==1);
chl = (AtabS(:,4)==2);
AtabS(sst,4) = 5;
AtabS(chl,4) = 6;

clear LAtab LFtab LPtab LDtab

%%  ---------------------------------------------------------
cnam = {'coef','p','lag','idriver','driver'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

%% Map
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

tlme = double(lme_mask);

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
clatlim=[plotminlat plotmaxlat];
clonlim=[plotminlon plotmaxlon];

load coastlines;

cmR = cbrewer('seq','Reds',10,'PCHIP');

%% Put R2 on grid

R2FS = nan(ni,nj);
R2PS = nan(ni,nj);
R2DS = nan(ni,nj);
R2AS = nan(ni,nj);

R2FD = nan(ni,nj);
R2PD = nan(ni,nj);
R2DD = nan(ni,nj);
R2AD = nan(ni,nj);

R2FC = nan(ni,nj);
R2PC = nan(ni,nj);
R2DC = nan(ni,nj);
R2AC = nan(ni,nj);

R2FO = nan(ni,nj);
R2PO = nan(ni,nj);
R2DO = nan(ni,nj);
R2AO = nan(ni,nj);

for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);
   
    if (AtabS(i,2) <= 0.05)
        R2AS(id) = AtabS(i,1).^2;
    end
    if (AtabD(i,2) <= 0.05)
        R2AD(id) = AtabD(i,1).^2;
    end
    if (AtabC(i,2) <= 0.05)
        R2AC(id) = AtabC(i,1).^2;
    end
    if (AtabO(i,2) <= 0.05)
        R2AO(id) = AtabO(i,1).^2;
    end

    if (FtabS(i,2) <= 0.05)
        R2FS(id) = FtabS(i,1).^2;
    end
    if (FtabD(i,2) <= 0.05)
        R2FD(id) = FtabD(i,1).^2;
    end
    if (FtabC(i,2) <= 0.05)
        R2FC(id) = FtabC(i,1).^2;
    end
    if (FtabO(i,2) <= 0.05)
        R2FO(id) = FtabO(i,1).^2;
    end

    if (PtabS(i,2) <= 0.05)
        R2PS(id) = PtabS(i,1).^2;
    end
    if (PtabD(i,2) <= 0.05)
        R2PD(id) = PtabD(i,1).^2;
    end
    if (PtabC(i,2) <= 0.05)
        R2PC(id) = PtabC(i,1).^2;
    end
    if (PtabO(i,2) <= 0.05)
        R2PO(id) = PtabO(i,1).^2;
    end

    if (DtabS(i,2) <= 0.05)
        R2DS(id) = DtabS(i,1).^2;
    end
    if (DtabD(i,2) <= 0.05)
        R2DD(id) = DtabD(i,1).^2;
    end
    if (DtabC(i,2) <= 0.05)
        R2DC(id) = DtabC(i,1).^2;
    end
    if (DtabO(i,2) <= 0.05)
        R2DO(id) = DtabO(i,1).^2;
    end
end

%% R2 of Dominant driver map - All
f1 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.05 0.555 0.45 0.4]) %Top L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2AS)
colormap(cmR);
clim([0 1])
title('CPUE Sat')

subplot('Position',[0.5 0.555 0.45 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2AD)
colormap(cmR);
clim([0 1])
title('CPUE Sat+OBGC')

subplot('Position',[0.05 0.105 0.45 0.4]) %Bottom L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2AC)
colormap(cmR);
clim([0 1])
title('CPUE Const Effort')

subplot('Position',[0.5 0.105 0.45 0.4]) %Bottom R
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2AO)
colormap(cmR);
clim([0 1])
colorbar('Position',[0.25 0.055 0.5 0.03],'orientation','horizontal')
title('CPUE Obs Effort')
print('-dpng',[ppath 'Map_LMEs_chlyr_cpue15_driver_comp_maxcorr_R2_All.png'])

%% R2 of Dominant driver map - Forage
f2 = figure('Units','inches','Position',[1 4 7.5 5]);
subplot('Position',[0.05 0.555 0.45 0.4]) %Top L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2FS)
colormap(cmR);
clim([0 1])
title('CPUE Sat')

subplot('Position',[0.5 0.555 0.45 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2FD)
colormap(cmR);
clim([0 1])
title('CPUE Sat+BGC')

subplot('Position',[0.05 0.105 0.45 0.4]) %Bottom L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2FC)
colormap(cmR);
clim([0 1])
title('CPUE Const Effort')

subplot('Position',[0.5 0.105 0.45 0.4]) %Bottom R
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2FO)
colormap(cmR);
clim([0 1])
colorbar('Position',[0.25 0.055 0.5 0.03],'orientation','horizontal')
title('CPUE Obs Effort')
print('-dpng',[ppath 'Map_LMEs_chlyr_cpue15_driver_comp_maxcorr_R2_Forage.png'])

%% R2 of Dominant driver map - Lg Pel
f3 = figure('Units','inches','Position',[1 5 7.5 5]);
subplot('Position',[0.05 0.555 0.45 0.4]) %Top L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2PS)
colormap(cmR);
clim([0 1])
title('CPUE Sat')

subplot('Position',[0.5 0.555 0.45 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2PD)
colormap(cmR);
clim([0 1])
title('CPUE Sat+BGC')

subplot('Position',[0.05 0.105 0.45 0.4]) %Bottom L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2PC)
colormap(cmR);
clim([0 1])
title('CPUE Const Effort')

subplot('Position',[0.5 0.105 0.45 0.4]) %Bottom R
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2PO)
colormap(cmR);
clim([0 1])
colorbar('Position',[0.25 0.055 0.5 0.03],'orientation','horizontal')
title('CPUE Obs Effort')
print('-dpng',[ppath 'Map_LMEs_chlyr_cpue15_driver_comp_maxcorr_R2_LgPel.png'])

%% R2 of Dominant driver map - Dem
f4 = figure('Units','inches','Position',[1 6 7.5 5]);

subplot('Position',[0.05 0.555 0.45 0.4]) %Top L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2DS)
colormap(cmR);
clim([0 1])
title('CPUE Sat')

subplot('Position',[0.5 0.555 0.45 0.4])
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2DD)
colormap(cmR);
clim([0 1])
title('CPUE Sat+BGC')

subplot('Position',[0.05 0.105 0.45 0.4]) %Bottom L
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2DC)
colormap(cmR);
clim([0 1])
title('CPUE Const Effort')

subplot('Position',[0.5 0.105 0.45 0.4]) %Bottom R
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2DO)
colormap(cmR);
clim([0 1])
colorbar('Position',[0.25 0.055 0.5 0.03],'orientation','horizontal')
title('CPUE Obs Effort')
print('-dpng',[ppath 'Map_LMEs_chlyr_cpue15_driver_comp_maxcorr_R2_Dem.png'])

%% Create 8 plot of 4 fn types w/ CPUE and Catch

subpos = [0.015 0.75 0.43 0.25;...
    0.015 0.5 0.43 0.25;...
    0.015 0.25 0.43 0.25;...
    0.015 0.0 0.43 0.25;...
    0.48 0.75 0.43 0.25;...
    0.48 0.5 0.43 0.25;...
    0.48 0.25 0.43 0.25;...
    0.48 0.0 0.43 0.25];

f10 = figure('Units','inches','Position',[1 3 6.5 8]);

% F Const
subplot('Position',subpos(1,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2FC)
colormap(cmR);
clim([0 1])
text(0,1.75,'F CPUE Const Effort','HorizontalAlignment','center')

% F Obs
subplot('Position',subpos(5,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2FO)
colormap(cmR);
clim([0 1])
text(0,1.75,'F CPUE Obs Effort','HorizontalAlignment','center')

% P Const
subplot('Position',subpos(2,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2PC)
colormap(cmR);
clim([0 1])
text(0,1.75,'P CPUE Const Effort','HorizontalAlignment','center')

%P Obs
subplot('Position',subpos(6,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2PO)
colormap(cmR);
clim([0 1])
text(0,1.75,'P CPUE Obs Effort','HorizontalAlignment','center')

% D Const
subplot('Position',subpos(3,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2DC)
colormap(cmR);
clim([0 1])
text(0,1.75,'D CPUE Const Effort','HorizontalAlignment','center')

%D Obs
subplot('Position',subpos(7,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2DO)
colormap(cmR);
clim([0 1])
text(0,1.75,'D CPUE Obs Effort','HorizontalAlignment','center')

% All CPUE Const
subplot('Position',subpos(4,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2AC)
colormap(cmR);
clim([0 1])
text(0,1.75,'All CPUE Const Effort','HorizontalAlignment','center')

% All CPUE Obs
subplot('Position',subpos(8,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,R2AO)
colormap(cmR);
clim([0 1])
colorbar('Position',[0.92 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
text(0,1.75,'All CPUE Obs Effort','HorizontalAlignment','center')

print('-dpng',[ppath 'Map_LMEs_chlyr_cpue15_feisty_obsfish2015_maxcorr_R2_types.png'])

