% Visualize output of FEISTY
% CESM Spinup first year of FOSI
% Time series plots and maps

clear all
close all

%% Zoo data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6.mat']);
load([cpath 'Data_grid_POP_gx1v6.mat']);
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';

load([cpath 'g.e11_LENS.GECOIAF.T62_g16.009.FIESTY-forcing.mat']);
load([cpath 'g.e11_LENS.GECOIAF.T62_g16.009.meszoo_v3.mat']);

%% Fish ids
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
mod = 'quad_v7_All_fish03';
% fish v6 = zoo v2
% fish v7 = zoo v3

fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
load([fpath 'Means_Spinup_',mod,'_lowbiom.mat']);   

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

%% Locations of low biom
% grid cell id > 10e4
fid2 = fid(fid>10e4);
pid2 = pid(pid>10e4);

%% 1st year
[ni,nj,nt] = size(diatC_150m);
diatC = diatC_150m(:,:,1:12);
LzooC = LzooC_150m(:,:,1:12);
Lzoo_loss = Lzoo_loss_150m(:,:,1:12);
Lzoo_quad = Lzoo_quad_150m(:,:,1:12);
zoo_loss = zoo_loss_150m(:,:,1:12);
fracDiat = fracL(:,:,1:12);

diatC = reshape(diatC,ni*nj,12);
LzooC = reshape(LzooC,ni*nj,12);
Lzoo_loss = reshape(Lzoo_loss,ni*nj,12);
Lzoo_quad = reshape(Lzoo_quad,ni*nj,12);
zoo_loss = reshape(zoo_loss,ni*nj,12);
fracDiat = reshape(fracDiat,ni*nj,12);

%%
diat = double(mean(diatC(fid2,:),1));
Lzoob = double(mean(LzooC(fid2,:),1));
Lzool = double(mean(Lzoo_loss(fid2,:),1));
Lzooq = double(mean(Lzoo_quad(fid2,:),1));
zool = double(mean(zoo_loss(fid2,:),1));
frac = double(mean(fracDiat(fid2,:),1));

%% Plots in time
y = 1:12;

figure(1)
subplot(2,1,1)
plot(y,log10(diat),'r','Linewidth',2); hold on;
plot(y,log10(Lzoob),'b','Linewidth',2); hold on;
legend('diatC','zooC')
legend('location','south')
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (mo)')
ylabel('log_1_0 Biomass/Loss')
title('FOSI year 1')

subplot(2,1,2)
plot(y,(zool),'r','Linewidth',2); hold on;
plot(y,(Lzool),'b','Linewidth',2); hold on;
plot(y,(Lzooq),'k','Linewidth',2); hold on;
legend('Zloss','LZloss','LZquad')
xlim([y(1) y(end)])
ylim([-1e-3 0.012])
xlabel('Time (mo)')
ylabel('Zoo Loss')
stamp(mod)
print('-dpng',[pp 'ts_mean_diat_zoo_FOSI_yr1_',mod,'_lowFbiom_cells.png'])

%%
rnum = round(length(fid2)*rand(40,1));

figure(3)
subplot(2,2,1)
plot(y,(Lzoo_quad(fid2(rnum(1:10)),:))); hold on;
xlim([y(1) y(end)])
ylim([-1e-3 0.03])
title('Quad loss (g m^-^2 d^-^1)')

subplot(2,2,2)
plot(y,(Lzoo_quad(fid2(rnum(11:20)),:))); hold on;
xlim([y(1) y(end)])
ylim([-1e-3 0.03])
title('low F locations')

subplot(2,2,3)
plot(y,(Lzoo_quad(fid2(rnum(21:30)),:))); hold on;
xlim([y(1) y(end)])
ylim([-1e-3 0.03])
xlabel('Time (mo)')

subplot(2,2,4)
plot(y,(Lzoo_quad(fid2(rnum(31:40)),:))); hold on;
xlim([y(1) y(end)])
ylim([-1e-3 0.03])
xlabel('Time (mo)')
stamp(mod)
print('-dpng',[pp 'ts_Lzoo_quad_FOSI_year1_',mod,'_rand_lowFbiom_cells.png'])

%%
figure(4)
subplot(2,2,1)
plot(y,(Lzoo_loss(fid2(rnum(1:10)),:))); hold on;
xlim([y(1) y(end)])
ylim([-1e-3 0.03])
title('LgZ loss (g m^-^2 d^-^1)')

subplot(2,2,2)
plot(y,(Lzoo_loss(fid2(rnum(11:20)),:))); hold on;
xlim([y(1) y(end)])
ylim([-1e-3 0.03])
title('low F locations')

subplot(2,2,3)
plot(y,(Lzoo_loss(fid2(rnum(21:30)),:))); hold on;
xlim([y(1) y(end)])
ylim([-1e-3 0.03])
xlabel('Time (mo)')

subplot(2,2,4)
plot(y,(Lzoo_loss(fid2(rnum(31:40)),:))); hold on;
xlim([y(1) y(end)])
ylim([-1e-3 0.03])
xlabel('Time (mo)')
stamp(mod)
print('-dpng',[pp 'ts_Lzoo_loss_FOSI_year1_',mod,'_rand_lowFbiom_cells.png'])

%%
figure(5)
subplot(2,2,1)
plot(y,(zoo_loss(fid2(rnum(1:10)),:))); hold on;
xlim([y(1) y(end)])
ylim([-1e-3 0.03])
title('Total loss (g m^-^2 d^-^1)')

subplot(2,2,2)
plot(y,(zoo_loss(fid2(rnum(11:20)),:))); hold on;
xlim([y(1) y(end)])
ylim([-1e-3 0.03])
title('low F locations')

subplot(2,2,3)
plot(y,(zoo_loss(fid2(rnum(21:30)),:))); hold on;
xlim([y(1) y(end)])
ylim([-1e-3 0.03])
xlabel('Time (mo)')

subplot(2,2,4)
plot(y,(zoo_loss(fid2(rnum(31:40)),:))); hold on;
xlim([y(1) y(end)])
ylim([-1e-3 0.03])
xlabel('Time (mo)')
stamp(mod)
print('-dpng',[pp 'ts_tzoo_loss_FOSI_year1_',mod,'_rand_lowFbiom_cells.png'])

%% Frac Lg
figure(6)
plot(y,frac,'color',[0 0.75 0.5],'Linewidth',2); hold on;
xlim([y(1) y(end)])
%ylim([-5 2])
xlabel('Time (mo)')
ylabel('Frac diat')
title('FOSI year 1')
stamp(mod)
print('-dpng',[pp 'ts_mean_diatFrac_FOSI_yr1_',mod,'_lowFbiom_cells.png'])

figure(7)
subplot(2,2,1)
plot(y,(fracDiat(fid2(rnum(1:10)),:))); hold on;
xlim([y(1) y(end)])
ylim([-0.1 1])
title('Diat Frac')

subplot(2,2,2)
plot(y,(fracDiat(fid2(rnum(11:20)),:))); hold on;
xlim([y(1) y(end)])
ylim([-0.1 1])
title('low F locations')

subplot(2,2,3)
plot(y,(fracDiat(fid2(rnum(21:30)),:))); hold on;
xlim([y(1) y(end)])
ylim([-0.1 1])
xlabel('Time (mo)')

subplot(2,2,4)
plot(y,(fracDiat(fid2(rnum(31:40)),:))); hold on;
xlim([y(1) y(end)])
xlabel('Time (mo)')
ylim([-0.1 1])
stamp(mod)
print('-dpng',[pp 'ts_fracL_FOSI_year1_',mod,'_rand_lowFbiom_cells.png'])

%%
zid = find(Lzoo_quad_150m==0);
length(zid) / (ni*nj*nt)

zid2 = find(Lzoo_quad_150m(:,:,1:12)==0);
length(zid2) / (ni*nj*12)



