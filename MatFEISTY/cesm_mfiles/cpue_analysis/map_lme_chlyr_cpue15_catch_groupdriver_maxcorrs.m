% Maps of driver with max corr with cpue & catch
% Group drivers by temp, prey, fish
% Subplots together for comparison
% Restricted analysis to chl yrs 1997-2015

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regress_cpue/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',...
    cfile,'/corrs_cpue/'];

mod = 'v15_All_fish03';
mod2 = 'v15_obsfish';

%% cpue const
load([spath,'LMEs_corr_cpue_chlyrs15_inputs_feisty_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab','lid')

EAtabC = LAtab;
EFtabC = LFtab;
EPtabC = LPtab;
EDtabC = LDtab;

clear LAtab LFtab LPtab LDtab

%% cpue obs
load([spath,'LMEs_corr_cpue_chlyrs15_driver_obsfish2015_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab') 

EAtabO = LAtab;
EFtabO = LFtab;
EPtabO = LPtab;
EDtabO = LDtab;

clear LAtab LFtab LPtab LDtab

%% catch const
load([spath,'LMEs_corr_catch_chlyrs15_inputs_feisty_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

CAtabC = LAtab;
CFtabC = LFtab;
CPtabC = LPtab;
CDtabC = LDtab;

clear LAtab LFtab LPtab LDtab

%% catch obs
load([spath,'LMEs_corr_catch_chlyrs15_inputs_obsfish2015_maxcorr_posfood.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

CAtabO = LAtab;
CFtabO = LFtab;
CPtabO = LPtab;
CDtabO = LDtab;

clear LAtab LFtab LPtab LDtab

%%  ---------------------------------------------------------
cnam = {'coef','p','lag','idriver','driver','group'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
ctex = {'TP','TB','Det','ZmLoss','SST','Chl','Biom','Prod','Yield'};

tid = [1,2,5];
pid = [3,4,6];
fid = 7:9;

%create new column with driver group
EAtabC(:,5) = nan;
eait = logical(sum((EAtabC(:,4)==tid),2));
eaip = logical(sum((EAtabC(:,4)==pid),2));
eaif = logical(sum((EAtabC(:,4)==fid),2));
EAtabC(eait,5) = 1;
EAtabC(eaip,5) = 2;
EAtabC(eaif,5) = 3;

EFtabC(:,5) = nan;
efit = logical(sum((EFtabC(:,4)==tid),2));
efip = logical(sum((EFtabC(:,4)==pid),2));
efif = logical(sum((EFtabC(:,4)==fid),2));
EFtabC(efit,5) = 1;
EFtabC(efip,5) = 2;
EFtabC(efif,5) = 3;

EPtabC(:,5) = nan;
epit = logical(sum((EPtabC(:,4)==tid),2));
epip = logical(sum((EPtabC(:,4)==pid),2));
epif = logical(sum((EPtabC(:,4)==fid),2));
EPtabC(epit,5) = 1;
EPtabC(epip,5) = 2;
EPtabC(epif,5) = 3;

EDtabC(:,5) = nan;
edit = logical(sum((EDtabC(:,4)==tid),2));
edip = logical(sum((EDtabC(:,4)==pid),2));
edif = logical(sum((EDtabC(:,4)==fid),2));
EDtabC(edit,5) = 1;
EDtabC(edip,5) = 2;
EDtabC(edif,5) = 3;

clear eait eaip eaif efit efip efif epit epip epif edit edip edif

%
EAtabO(:,5) = nan;
eait = logical(sum((EAtabO(:,4)==tid),2));
eaip = logical(sum((EAtabO(:,4)==pid),2));
eaif = logical(sum((EAtabO(:,4)==fid),2));
EAtabO(eait,5) = 1;
EAtabO(eaip,5) = 2;
EAtabO(eaif,5) = 3;

EFtabO(:,5) = nan;
efit = logical(sum((EFtabO(:,4)==tid),2));
efip = logical(sum((EFtabO(:,4)==pid),2));
efif = logical(sum((EFtabO(:,4)==fid),2));
EFtabO(efit,5) = 1;
EFtabO(efip,5) = 2;
EFtabO(efif,5) = 3;

EPtabO(:,5) = nan;
epit = logical(sum((EPtabO(:,4)==tid),2));
epip = logical(sum((EPtabO(:,4)==pid),2));
epif = logical(sum((EPtabO(:,4)==fid),2));
EPtabO(epit,5) = 1;
EPtabO(epip,5) = 2;
EPtabO(epif,5) = 3;

EDtabO(:,5) = nan;
edit = logical(sum((EDtabO(:,4)==tid),2));
edip = logical(sum((EDtabO(:,4)==pid),2));
edif = logical(sum((EDtabO(:,4)==fid),2));
EDtabO(edit,5) = 1;
EDtabO(edip,5) = 2;
EDtabO(edif,5) = 3;

clear eait eaip eaif efit efip efif epit epip epif edit edip edif

% Catch const
CAtabC(:,5) = nan;
eait = logical(sum((CAtabC(:,4)==tid),2));
eaip = logical(sum((CAtabC(:,4)==pid),2));
eaif = logical(sum((CAtabC(:,4)==fid),2));
CAtabC(eait,5) = 1;
CAtabC(eaip,5) = 2;
CAtabC(eaif,5) = 3;

CFtabC(:,5) = nan;
efit = logical(sum((CFtabC(:,4)==tid),2));
efip = logical(sum((CFtabC(:,4)==pid),2));
efif = logical(sum((CFtabC(:,4)==fid),2));
CFtabC(efit,5) = 1;
CFtabC(efip,5) = 2;
CFtabC(efif,5) = 3;

CPtabC(:,5) = nan;
epit = logical(sum((CPtabC(:,4)==tid),2));
epip = logical(sum((CPtabC(:,4)==pid),2));
epif = logical(sum((CPtabC(:,4)==fid),2));
CPtabC(epit,5) = 1;
CPtabC(epip,5) = 2;
CPtabC(epif,5) = 3;

CDtabC(:,5) = nan;
edit = logical(sum((CDtabC(:,4)==tid),2));
edip = logical(sum((CDtabC(:,4)==pid),2));
edif = logical(sum((CDtabC(:,4)==fid),2));
CDtabC(edit,5) = 1;
CDtabC(edip,5) = 2;
CDtabC(edif,5) = 3;

clear eait eaip eaif efit efip efif epit epip epif edit edip edif

%
CAtabO(:,5) = nan;
eait = logical(sum((CAtabO(:,4)==tid),2));
eaip = logical(sum((CAtabO(:,4)==pid),2));
eaif = logical(sum((CAtabO(:,4)==fid),2));
CAtabO(eait,5) = 1;
CAtabO(eaip,5) = 2;
CAtabO(eaif,5) = 3;

CFtabO(:,5) = nan;
efit = logical(sum((CFtabO(:,4)==tid),2));
efip = logical(sum((CFtabO(:,4)==pid),2));
efif = logical(sum((CFtabO(:,4)==fid),2));
CFtabO(efit,5) = 1;
CFtabO(efip,5) = 2;
CFtabO(efif,5) = 3;

CPtabO(:,5) = nan;
epit = logical(sum((CPtabO(:,4)==tid),2));
epip = logical(sum((CPtabO(:,4)==pid),2));
epif = logical(sum((CPtabO(:,4)==fid),2));
CPtabO(epit,5) = 1;
CPtabO(epip,5) = 2;
CPtabO(epif,5) = 3;

CDtabO(:,5) = nan;
edit = logical(sum((CDtabO(:,4)==tid),2));
edip = logical(sum((CDtabO(:,4)==pid),2));
edif = logical(sum((CDtabO(:,4)==fid),2));
CDtabO(edit,5) = 1;
CDtabO(edip,5) = 2;
CDtabO(edif,5) = 3;

clear eait eaip eaif efit efip efif epit epip epif edit edip edif

%% colorblind friendly
load('paul_tol_cmaps.mat')

%colorblind friendly - subselect & re-order drainbow

%group by temp, resources, fish
gtex = {'Temp','Res','Fish'};
gcol(1,:) = drainbow(12,:)/255; % orange
gcol(2,:) = drainbow(7,:)/255; %green 
gcol(3,:) = drainbow(3,:)/255; %dk purp

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

%% Put on grid

EACcorr  = nan(ni,nj);
EFCcorr  = nan(ni,nj);
EPCcorr  = nan(ni,nj);
EDCcorr  = nan(ni,nj);
EAOcorr  = nan(ni,nj);
EFOcorr  = nan(ni,nj);
EPOcorr  = nan(ni,nj);
EDOcorr  = nan(ni,nj);
CACcorr  = nan(ni,nj);
CFCcorr  = nan(ni,nj);
CPCcorr  = nan(ni,nj);
CDCcorr  = nan(ni,nj);
CAOcorr  = nan(ni,nj);
CFOcorr  = nan(ni,nj);
CPOcorr  = nan(ni,nj);
CDOcorr  = nan(ni,nj);

for i=1:length(lid)
    L=lid(i);
    id = find(tlme==L);

    if (EAtabC(i,2) <= 0.05)
        EACcorr(id) = EAtabC(i,5);
    end
    if (EFtabC(i,2) <= 0.05)
        EFCcorr(id) = EFtabC(i,5);
    end
    if (EPtabC(i,2) <= 0.05)
        EPCcorr(id) = EPtabC(i,5);
    end
    if (EDtabC(i,2) <= 0.05)
        EDCcorr(id) = EDtabC(i,5);
    end
    if (EAtabO(i,2) <= 0.05)
        EAOcorr(id) = EAtabO(i,5);
    end
    if (EFtabO(i,2) <= 0.05)
        EFOcorr(id) = EFtabO(i,5);
    end
    if (EPtabO(i,2) <= 0.05)
        EPOcorr(id) = EPtabO(i,5);
    end
    if (EDtabO(i,2) <= 0.05)
        EDOcorr(id) = EDtabO(i,5);
    end

    
    if (CAtabC(i,2) <= 0.05)
        CACcorr(id) = CAtabC(i,5);
    end
    if (CFtabC(i,2) <= 0.05)
        CFCcorr(id) = CFtabC(i,5);
    end
    if (CPtabC(i,2) <= 0.05)
        CPCcorr(id) = CPtabC(i,5);
    end
    if (CDtabC(i,2) <= 0.05)
        CDCcorr(id) = CDtabC(i,5);
    end
    if (CAtabO(i,2) <= 0.05)
        CAOcorr(id) = CAtabO(i,5);
    end
    if (CFtabO(i,2) <= 0.05)
        CFOcorr(id) = CFtabO(i,5);
    end
    if (CPtabO(i,2) <= 0.05)
        CPOcorr(id) = CPtabO(i,5);
    end
    if (CDtabO(i,2) <= 0.05)
        CDOcorr(id) = CDtabO(i,5);
    end

    
end

%% CPUE Dominant driver map
subpos = [0.015 0.75 0.43 0.25;...
    0.015 0.5 0.43 0.25;...
    0.015 0.25 0.43 0.25;...
    0.015 0.0 0.43 0.25;...
    0.46 0.75 0.43 0.25;...
    0.46 0.5 0.43 0.25;...
    0.46 0.25 0.43 0.25;...
    0.46 0.0 0.43 0.25];

f1 = figure('Units','inches','Position',[1 3 6.5 8]);

% F Const
subplot('Position',subpos(1,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,EFCcorr)
colormap(gcol);
text(0,1.75,'Forage Const Effort','HorizontalAlignment','center')

% F Obs
subplot('Position',subpos(5,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,EFOcorr)
colormap(gcol);
text(0,1.75,'Forage Obs Effort','HorizontalAlignment','center')

% P Const
subplot('Position',subpos(2,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,EPCcorr)
colormap(gcol);
text(0,1.75,'LgPel Const Effort','HorizontalAlignment','center')

%P Obs
subplot('Position',subpos(6,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,EPOcorr)
colormap(gcol);
text(0,1.75,'LgPel Obs Effort','HorizontalAlignment','center')
text(3.5,2.75,'CPUE','HorizontalAlignment','center')

% D Const
subplot('Position',subpos(3,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,EDCcorr)
colormap(gcol);
text(0,1.75,'Demers Const Effort','HorizontalAlignment','center')

%D Obs
subplot('Position',subpos(7,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,EDOcorr)
colormap(gcol);
text(0,1.75,'Demers Obs Effort','HorizontalAlignment','center')

% All CPUE Const
subplot('Position',subpos(4,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,EACcorr)
colormap(gcol);
text(0,1.75,'All Const Effort','HorizontalAlignment','center')

% All CPUE Obs
subplot('Position',subpos(8,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,EAOcorr)
colormap(gcol);
%clim([0 1])
colorbar('Ticks',1:3,'TickLabels',gtex, 'Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out');
text(0,1.75,'All Obs Effort','HorizontalAlignment','center')

print('-dpng',[ppath 'Map_LMEs_chlyr15_cpue15_feisty_obsfish2015_groupdriver_types.png'])


%% Catch Dominant driver map

f2 = figure('Units','inches','Position',[1 3 6.5 8]);

% F Const
subplot('Position',subpos(1,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,CFCcorr)
colormap(gcol);
text(0,1.75,'Forage Const Effort','HorizontalAlignment','center')

% F Obs
subplot('Position',subpos(5,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,CFOcorr)
colormap(gcol);
text(0,1.75,'Forage Obs Effort','HorizontalAlignment','center')

% P Const
subplot('Position',subpos(2,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,CPCcorr)
colormap(gcol);
text(0,1.75,'LgPel Const Effort','HorizontalAlignment','center')

%P Obs
subplot('Position',subpos(6,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,CPOcorr)
colormap(gcol);
text(0,1.75,'LgPel Obs Effort','HorizontalAlignment','center')
text(3.5,2.75,'Catch','HorizontalAlignment','center')

% D Const
subplot('Position',subpos(3,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,CDCcorr)
colormap(gcol);
text(0,1.75,'Demers Const Effort','HorizontalAlignment','center')

%D Obs
subplot('Position',subpos(7,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,CDOcorr)
colormap(gcol);
text(0,1.75,'Demers Obs Effort','HorizontalAlignment','center')

% All CPUE Const
subplot('Position',subpos(4,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,CACcorr)
colormap(gcol);
text(0,1.75,'All Const Effort','HorizontalAlignment','center')

% All CPUE Obs
subplot('Position',subpos(8,:))
axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hold on
surfm(TLAT,TLONG,CAOcorr)
colormap(gcol);
%clim([0 1])
colorbar('Ticks',1:3,'TickLabels',gtex, 'Position',[0.9 0.25 0.025 0.5],'orientation','vertical','AxisLocation','out')
text(0,1.75,'All Obs Effort','HorizontalAlignment','center')

print('-dpng',[ppath 'Map_LMEs_chlyr15_catch_feisty_obsfish2015_groupdriver_types.png'])



