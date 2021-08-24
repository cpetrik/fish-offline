% Visualize output of FEISTY
% CESM Spinup first year of FOSI
% Time series plots and maps

clear all
close all

%% Zoo data
% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6.mat']);
load([cpath 'Data_grid_POP_gx1v6.mat']);

load([cpath 'Data_cesm_fosi_quad_v3.3_daily_1.mat']);

[ni,nj]=size(TLONG);

%% Fish ids
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
mod = 'quad_v11';
% fish v6,8,10 = zoo v2, v2.2, v2.3
% fish v7,9,11 = zoo v3, v3.2, v3.3

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];

%orig loactions of low biomass
mod0 = 'quad_v7_All_fish03';
load([fpath 'Means_Spinup_',mod0,'_lowbiom.mat']);

%%
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;     

%% colors
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

%cmBP50=cbrewer('seq','BuPu',50,'PCHIP');
cmBP50=cbrewer('seq','YlGnBu',10,'PCHIP');

%% Locations of low biom
lowF = zeros(size(TLONG));
lowF(fid) = ones;

lowP = zeros(size(TLONG));
lowP(pid) = ones;

Fid = zeros(size(TLONG));
Fid(fid) = fid;

Pid = zeros(size(TLONG));
Pid(pid) = pid;

%% Time series of those cells
% grid cell id > 10e4
fid2 = fid(fid>10e4);
pid2 = pid(pid>10e4);
[fx,fid3] = intersect(GRD.ID,fid2);
[px,pid3] = intersect(GRD.ID,pid2);

%%
Fzb_tmean=mean(ESM.Zm(fid3,:),1);
Pzb_tmean=mean(ESM.Zm(pid3,:),1);
Fzq_tmean=mean(ESM.dZm(fid3,:),1);
Pzq_tmean=mean(ESM.dZm(pid3,:),1);

%% Plots in time
t = 1:length(Fzb_tmean); %time;
y = t;

%% Types together

figure(1)
subplot(2,1,1)
plot(y,log10(Fzb_tmean),'r','Linewidth',2); hold on;
plot(y,log10(Pzb_tmean),'b','Linewidth',2); hold on;
legend('low F','low P')
legend('location','south')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (y)')
ylabel('log_1_0 Biomass (g m^-^2)')
title('FOSI day 1')

subplot(2,1,2)
plot(y,log10(Fzq_tmean),'r','Linewidth',2); hold on;
plot(y,log10(Pzq_tmean),'b','Linewidth',2); hold on;
xlim([y(1) y(end)])
%ylim([-0.05 0.8])
xlabel('Time (y)')
ylabel('log_1_0 Quad loss (g m^-^2 d^-^1)')
stamp(mod)
print('-dpng',[pp 'ts_mean_zoo_FOSI_day1_',mod,'_lowbiom_cells.png'])

%%
rnum = round(length(fid3)*rand(40,1));

figure(3)
subplot(2,2,1)
plot(y,(ESM.dZm(fid3(rnum(1:10)),:))); hold on;
xlim([y(1) y(end)])
%ylim([-8 0])
xlabel('Time (y)')
title('Quad loss (g m^-^2 d^-^1)')

subplot(2,2,2)
plot(y,(ESM.dZm(fid3(rnum(11:20)),:))); hold on;
xlim([y(1) y(end)])
%ylim([-8 0])
xlabel('Time (y)')
title('low F locations')

subplot(2,2,3)
plot(y,(ESM.dZm(fid3(rnum(21:30)),:))); hold on;
xlim([y(1) y(end)])
%ylim([-8 0])
xlabel('Time (y)')

subplot(2,2,4)
plot(y,(ESM.dZm(fid3(rnum(31:40)),:))); hold on;
xlim([y(1) y(end)])
%ylim([-8 0])
xlabel('Time (y)')
stamp(mod)
print('-dpng',[pp 'ts_zoo_quad_FOSI_day1_',mod,'_rand_lowFbiom_cells.png'])
 

