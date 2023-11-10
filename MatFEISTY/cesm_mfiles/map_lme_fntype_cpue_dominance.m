% Use fish rel cpue & cme to define 3-4 ecosys/foodweb structures
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/' cfile '/corrs/'];

mod = 'v15_All_fish03_';

%% catch/effort
%fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/Fish-MIP/Phase3/fishing/';
fpath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

load([fpath 'FishMIP_Phase3a_LME_CPUE_1961-2010_interann_var.mat'],...
    'cpue_f_mean','cpue_p_mean','cpue_d_mean','cpue_a_mean');

load([fpath 'FishMIP_Phase3a_LME_catch_minus_effort_1961-2010_interann_var.mat'],...
    'cme_f_mean','cme_p_mean','cme_d_mean','cme_a_mean');

%%  ---------------------------------------------------------
ftex = {'F','P','D','A','B'};
% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)

% mcol = [238/255 102/255 119/255;... %red
%     0/255 68/255 136/255;...    %blue
%     34/255 136/255 51/255;...   %green
%     51/255 187/255 238/255;...  %cyan
%     170/255 51/255 119/255;...  %purple
%     ];

% colorblind friendly
load('paul_tol_cmaps.mat')

%ecosystem type
mcol = drainbow(2:3:end,:) ./ 255;

set(groot,'defaultAxesColorOrder',mcol);

%% Relative amounts
cpue_type(:,1) = cpue_f_mean;
cpue_type(:,2) = cpue_p_mean;
cpue_type(:,3) = cpue_d_mean;
cpue_btot = sum(cpue_type,2);
fcpue = cpue_type ./ repmat(cpue_btot,1,3);

cme_type(:,1) = cme_f_mean;
cme_type(:,2) = cme_p_mean;
cme_type(:,3) = cme_d_mean;
cme_btot = sum(cme_type,2);
fcme = cme_type ./ repmat(cme_btot,1,3);

%% Define foodweb types
% F dom: F>=0.4     %1
% D&F dom: P<0.2    %2
% D dom: D>=0.5     %3
% P&F dom: D<0.2    %4
% even              %5

Cetype = 5*ones(66,1);
Cetype(fcpue(:,1)>=0.4) = 1;
Cetype(fcpue(:,2)<0.2) = 2;
Cetype(fcpue(:,3)>=0.5) = 3;
Cetype(fcpue(:,3)<0.2) = 4;

Metype = 5*ones(66,1);
Metype(fcme(:,1)>=0.4) = 1;
Metype(fcme(:,2)<0.2) = 2;
Metype(fcme(:,3)>=0.5) = 3;
Metype(fcme(:,3)<0.2) = 4;

%NaNs are 23, 33, 62 (inland seas)
Cetype([23, 33, 62],1) = nan;
Metype([23, 33, 62],1) = nan;

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
Cbiom  = nan(ni,nj);
Mbiom  = nan(ni,nj);

for i=1:66
    L=i;
    id = find(tlme==L);

    Cbiom(id) = Cetype(i,1);
    Mbiom(id) = Metype(i,1);
    
end

%%
f1 = figure('Units','inches','Position',[1 3 7.5 5]);
subplot('Position',[0.05 0.5 0.4 0.45]) 
ax1=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Cbiom)
colormap(ax1,mcol)
caxis([1 5]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('CPUE structure')
%colorbar('Ticks',1.5:0.75:5,'TickLabels',alltex)

subplot('Position',[0.5 0.5 0.4 0.45]) 
ax1=axesm ('Robinson','MapLatLimit',clatlim,'MapLonLimit',clonlim,'frame','on',...
    'Grid','off','FLineWidth',1);
surfm(TLAT,TLONG,Mbiom)
colormap(ax1,mcol)
caxis([1 5]);
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
title('CME structure')
colorbar('Ticks',1.5:0.75:5,'TickLabels',alltex)

print('-dpng',[ppath 'Map_LMEs_rel_biom_ecosys_types.png'])

%% save etypes for mapping with other results
ftex = alltex;

fpath1 ='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/Fish-MIP/Phase3/fishing/';
fpath2 = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

save([fpath1 'FishMIP_Phase3a_LME_catch_effort_1961-2010_fntype_dom.mat'],...
    'Cetype','Metype','ftex','Cbiom','Mbiom');

save([fpath2 'FishMIP_Phase3a_LME_catch_effort_1961-2010_fntype_dom.mat'],...
    'Cetype','Metype','ftex','Cbiom','Mbiom');

