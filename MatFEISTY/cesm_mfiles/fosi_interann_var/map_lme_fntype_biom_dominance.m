% Use fish rel biomass to define 3-4 ecosys/foodweb structures
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/' cfile '/corrs/'];

mod = 'v15_All_fish03_';

%% inputs
load([cpath 'lme_means_g.e11_LENS.GECOIAF.T62_g16.009.mat'],...
    'lme_tp_fosi','lme_tb_fosi','lme_det_fosi','lme_mz_fosi','lme_mzloss_fosi')

% fish
load([fpath 'LME_fosi_fished_',mod,cfile '.mat'],'lme_area','lme_mtype',...
    'lme_mnu');

lme_nu = lme_mnu;

%%  ---------------------------------------------------------
ftex = {'F','P','D','A','B'};
ctex = {'TP','TB','Det','Zmeso','ZmLoss'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

% US LMEs
lid = [54,1:2,10,3,5:7,65]; %ADD 65 = Aleutian Islands
lname = {'CHK','EBS','GAK','HI','CCE','GMX','SE','NE','AI'};

load('paul_tol_cmaps.mat')

mcol = [238/255 102/255 119/255;... %red
    0/255 68/255 136/255;...    %blue
    34/255 136/255 51/255;...   %green
    51/255 187/255 238/255;...  %cyan
    170/255 51/255 119/255;...  %purple
    ];

%colorblind friendly
% cb=[34/255 136/255 51/255;...   %green
%     153/255 153/255 51/255;...  %olive
%     51/255 187/255 238/255;...  %cyan
%     0/255 68/255 136/255;...    %blue
%     238/255 102/255 119/255;... %red
%     170/255 51/255 119/255;...  %purple
%     0 0 0;...                   %black
%     0.25 0.25 0.25;...             %dk grey
%     0.50 0.50 0.50;...             % grey
%     0.75 0.75 0.75];               %lt grey

set(groot,'defaultAxesColorOrder',mcol);

cmB=cbrewer('seq','Blues',50,'PCHIP');
cmR=cbrewer('seq','Reds',50,'PCHIP');
cmG=cbrewer('seq','Greens',50,'PCHIP');

%% Relative amounts
lme_type = lme_mtype(:,1:3);
lme_btot = sum(lme_type,2);

lbio = lme_type ./ repmat(lme_btot,1,3);

%% Define foodweb types
% F dom: F>=0.4     %1
% D&F dom: P<0.2    %2
% D dom: D>=0.5     %3
% P&F dom: D<0.2    %4
% even              %5

etype = 5*ones(66,1);
etype(lbio(:,1)>=0.4) = 1;
etype(lbio(:,2)<0.2) = 2;
etype(lbio(:,3)>=0.5) = 3;
etype(lbio(:,3)<0.2) = 4;

%NaNs are 23, 33, 62 (inland seas)
etype([23, 33, 62],1) = nan;

alltex = {'Forage',...    % = 1
    'F & D',...          % = 2
    'Demersal',...           % = 3
    'F & P',...           % = 4
    'Even'};             % = 5

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

%%
% on grid
Ebiom  = nan(ni,nj);

for i=1:66
    L=i;
    id = find(tlme==L);

    Ebiom(id) = etype(i,1);
    
end

%%
f1 = figure('Units','inches','Position',[1 3 7.5 5]);
ax1=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Ebiom)
colormap(ax1,mcol)
caxis([1 5]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('Food web structure')
colorbar('Ticks',1.5:0.75:5,'TickLabels',alltex)

print('-dpng',[ppath 'Map_LMEs_rel_biom_ecosys_types.png'])

%% bright
bri = bright./255;
f2 = figure('Units','inches','Position',[1 3 7.5 5]);
ax1=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Ebiom)
colormap(ax1,bri(1:5,:))
caxis([1 5]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('Food web structure')
colorbar('Ticks',1.5:0.75:5,'TickLabels',alltex)
print('-dpng',[ppath 'Map_LMEs_rel_biom_ecosys_types_bright.png'])

%% night
nit = flipud(night(1:3:end,:)./255);
f3 = figure('Units','inches','Position',[1 3 7.5 5]);
ax1=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Ebiom)
colormap(ax1,nit(2:6,:))
caxis([1 5]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('Food web structure')
colorbar('Ticks',1.5:0.75:5,'TickLabels',alltex)
print('-dpng',[ppath 'Map_LMEs_rel_biom_ecosys_types_night.png'])

%% save etypes for mapping with other results
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

etex = alltex;

save([fpath 'LME_fosi_fished_',mod,cfile '.mat'],'etype','etex',...
    'Ebiom','-append');

