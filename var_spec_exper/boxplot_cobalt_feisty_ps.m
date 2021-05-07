% Distribution plots of power sepctra from CORE-forced
% Both forcing (phys-BGC) and fish
% Look for patterns of incr phys -> BGC -> Sfish -> Mfish -> Lfish
% Summarize global and by biomes?

clear all
close all

%% COBALT --------------------------------------------------------
spath='/Volumes/MIP/GCM_DATA/CORE-forced/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([cpath 'COBALT_hist_biomes_1950_2005.mat']);

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
ppath = [pp cfile '/'];

%% powerspectra
load([spath 'powerspec_cobalt_core_1950_2007.mat']);
load([fpath 'powerspec_feisty_core_1950_2007.mat']);

biome = {'LC','ECCS','ECSS','Coast'};

lid = find(biome4_hist==1);
uid = find(biome4_hist==2);
sid = find(biome4_hist==3);
cid = find(biome4_hist==4);

%% boxplots - Global
figure(1)
subplot(2,2,1)
boxplot([tp_ps(GRD.ID),tb_ps(GRD.ID),det_ps(GRD.ID),zm_ps(GRD.ID),hpmz_ps(GRD.ID),...
    zl_ps(GRD.ID),hplz_ps(GRD.ID)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TP','TB','Det','MZ','MZhp','LZ','LZhp'})
ylabel('PS Slope')
title('Forcing')

subplot(2,2,3)
boxplot([sf_ps(GRD.ID),sp_ps(GRD.ID),sd_ps(GRD.ID),mf_ps(GRD.ID),mp_ps(GRD.ID),...
    md_ps(GRD.ID),lp_ps(GRD.ID),ld_ps(GRD.ID)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'SF','SP','SD','MF','MP','MD','LP','LD'})
ylabel('PS Slope')
title('Stages')

subplot(2,2,4)
boxplot([F_ps(GRD.ID),P_ps(GRD.ID),D_ps(GRD.ID),S_ps(GRD.ID),M_ps(GRD.ID),...
    L_ps(GRD.ID),All_ps(GRD.ID)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'F','P','D','S','M','L','All'})
ylabel('PS Slope')
title('Groups')
print('-dpng',[ppath 'CORE_boxplot_ps_global.png'])

%%
figure(2)
subplot(2,1,1)
boxplot([tp_ps(GRD.ID),zm_ps(GRD.ID),hpmz_ps(GRD.ID),zl_ps(GRD.ID),...
    hplz_ps(GRD.ID),sf_ps(GRD.ID),sp_ps(GRD.ID),sd_ps(GRD.ID),mf_ps(GRD.ID),...
    mp_ps(GRD.ID),lp_ps(GRD.ID)],...
    'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TP','MZ','MZhp','LZ','LZhp','SF','SP','SD','MF','MP','LP'})
ylabel('PS Slope')
title('Pelagic')

%
subplot(2,2,3)
boxplot([tb_ps(GRD.ID),det_ps(GRD.ID),md_ps(GRD.ID),ld_ps(GRD.ID)],...
    'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TB','Det','MD','LD'})
ylabel('PS Slope')
title('Benthic')
print('-dpng',[ppath 'CORE_boxplot_ps_global_habitat.png'])

%% LC
figure(3)
subplot(2,2,1)
boxplot([tp_ps(lid),tb_ps(lid),det_ps(lid),zm_ps(lid),hpmz_ps(lid),...
    zl_ps(lid),hplz_ps(lid)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TP','TB','Det','MZ','MZhp','LZ','LZhp'})
ylabel('PS Slope')
title('Forcing')

subplot(2,2,3)
boxplot([sf_ps(lid),sp_ps(lid),sd_ps(lid),mf_ps(lid),mp_ps(lid),...
    md_ps(lid),lp_ps(lid),ld_ps(lid)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'SF','SP','SD','MF','MP','MD','LP','LD'})
ylabel('PS Slope')
title('Stages')

subplot(2,2,4)
boxplot([F_ps(lid),P_ps(lid),D_ps(lid),S_ps(lid),M_ps(lid),...
    L_ps(lid),All_ps(lid)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'F','P','D','S','M','L','All'})
ylabel('PS Slope')
title('Groups')
print('-dpng',[ppath 'CORE_boxplot_ps_LC.png'])

%
figure(4)
subplot(2,1,1)
boxplot([tp_ps(lid),zm_ps(lid),hpmz_ps(lid),zl_ps(lid),...
    hplz_ps(lid),sf_ps(lid),sp_ps(lid),sd_ps(lid),mf_ps(lid),...
    mp_ps(lid),lp_ps(lid)],...
    'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TP','MZ','MZhp','LZ','LZhp','SF','SP','SD','MF','MP','LP'})
ylabel('PS Slope')
title('Pelagic')

%
subplot(2,2,3)
boxplot([tb_ps(lid),det_ps(lid),md_ps(lid),ld_ps(lid)],...
    'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TB','Det','MD','LD'})
ylabel('PS Slope')
title('Benthic')
print('-dpng',[ppath 'CORE_boxplot_ps_LC_habitat.png'])

%% ECCS
figure(5)
subplot(2,2,1)
boxplot([tp_ps(uid),tb_ps(uid),det_ps(uid),zm_ps(uid),hpmz_ps(uid),...
    zl_ps(uid),hplz_ps(uid)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TP','TB','Det','MZ','MZhp','LZ','LZhp'})
ylabel('PS Slope')
title('Forcing')

subplot(2,2,3)
boxplot([sf_ps(uid),sp_ps(uid),sd_ps(uid),mf_ps(uid),mp_ps(uid),...
    md_ps(uid),lp_ps(uid),ld_ps(uid)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'SF','SP','SD','MF','MP','MD','LP','LD'})
ylabel('PS Slope')
title('Stages')

subplot(2,2,4)
boxplot([F_ps(uid),P_ps(uid),D_ps(uid),S_ps(uid),M_ps(uid),...
    L_ps(uid),All_ps(uid)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'F','P','D','S','M','L','All'})
ylabel('PS Slope')
title('Groups')
print('-dpng',[ppath 'CORE_boxplot_ps_eccs.png'])

%
figure(6)
subplot(2,1,1)
boxplot([tp_ps(uid),zm_ps(uid),hpmz_ps(uid),zl_ps(uid),...
    hplz_ps(uid),sf_ps(uid),sp_ps(uid),sd_ps(uid),mf_ps(uid),...
    mp_ps(uid),lp_ps(uid)],...
    'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TP','MZ','MZhp','LZ','LZhp','SF','SP','SD','MF','MP','LP'})
ylabel('PS Slope')
title('Pelagic')

%
subplot(2,2,3)
boxplot([tb_ps(uid),det_ps(uid),md_ps(uid),ld_ps(uid)],...
    'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TB','Det','MD','LD'})
ylabel('PS Slope')
title('Benthic')
print('-dpng',[ppath 'CORE_boxplot_ps_eccs_habitat.png'])

%% ECSS
figure(7)
subplot(2,2,1)
boxplot([tp_ps(sid),tb_ps(sid),det_ps(sid),zm_ps(sid),hpmz_ps(sid),...
    zl_ps(sid),hplz_ps(sid)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TP','TB','Det','MZ','MZhp','LZ','LZhp'})
ylabel('PS Slope')
title('Forcing')

subplot(2,2,3)
boxplot([sf_ps(sid),sp_ps(sid),sd_ps(sid),mf_ps(sid),mp_ps(sid),...
    md_ps(sid),lp_ps(sid),ld_ps(sid)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'SF','SP','SD','MF','MP','MD','LP','LD'})
ylabel('PS Slope')
title('Stages')

subplot(2,2,4)
boxplot([F_ps(sid),P_ps(sid),D_ps(sid),S_ps(sid),M_ps(sid),...
    L_ps(sid),All_ps(sid)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'F','P','D','S','M','L','All'})
ylabel('PS Slope')
title('Groups')
print('-dpng',[ppath 'CORE_boxplot_ps_ecss.png'])

%
figure(8)
subplot(2,1,1)
boxplot([tp_ps(sid),zm_ps(sid),hpmz_ps(sid),zl_ps(sid),...
    hplz_ps(sid),sf_ps(sid),sp_ps(sid),sd_ps(sid),mf_ps(sid),...
    mp_ps(sid),lp_ps(sid)],...
    'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TP','MZ','MZhp','LZ','LZhp','SF','SP','SD','MF','MP','LP'})
ylabel('PS Slope')
title('Pelagic')

%
subplot(2,2,3)
boxplot([tb_ps(sid),det_ps(sid),md_ps(sid),ld_ps(sid)],...
    'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TB','Det','MD','LD'})
ylabel('PS Slope')
title('Benthic')
print('-dpng',[ppath 'CORE_boxplot_ps_ecss_habitat.png'])

%% Coastal
figure(9)
subplot(2,2,1)
boxplot([tp_ps(cid),tb_ps(cid),det_ps(cid),zm_ps(cid),hpmz_ps(cid),...
    zl_ps(cid),hplz_ps(cid)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TP','TB','Det','MZ','MZhp','LZ','LZhp'})
ylabel('PS Slope')
title('Forcing')

subplot(2,2,3)
boxplot([sf_ps(cid),sp_ps(cid),sd_ps(cid),mf_ps(cid),mp_ps(cid),...
    md_ps(cid),lp_ps(cid),ld_ps(cid)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'SF','SP','SD','MF','MP','MD','LP','LD'})
ylabel('PS Slope')
title('Stages')

subplot(2,2,4)
boxplot([F_ps(cid),P_ps(cid),D_ps(cid),S_ps(cid),M_ps(cid),...
    L_ps(cid),All_ps(cid)],'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'F','P','D','S','M','L','All'})
ylabel('PS Slope')
title('Groups')
print('-dpng',[ppath 'CORE_boxplot_ps_coast.png'])

%
figure(10)
subplot(2,1,1)
boxplot([tp_ps(cid),zm_ps(cid),hpmz_ps(cid),zl_ps(cid),...
    hplz_ps(cid),sf_ps(cid),sp_ps(cid),sd_ps(cid),mf_ps(cid),...
    mp_ps(cid),lp_ps(cid)],...
    'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TP','MZ','MZhp','LZ','LZhp','SF','SP','SD','MF','MP','LP'})
ylabel('PS Slope')
title('Pelagic')

%
subplot(2,2,3)
boxplot([tb_ps(cid),det_ps(cid),md_ps(cid),ld_ps(cid)],...
    'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
set(gca,'XTickLabel',{'TB','Det','MD','LD'})
ylabel('PS Slope')
title('Benthic')
print('-dpng',[ppath 'CORE_boxplot_ps_coast_habitat.png'])

%%
lab1 = repmat({'Glb'},length(GRD.ID),1);
lab2 = repmat({'LC'},length(lid),1);
lab3 = repmat({'CS'},length(uid),1);
lab4 = repmat({'SS'},length(sid),1);
lab5 = repmat({'Co'},length(cid),1);
blab = [lab1;lab2;lab3;lab4;lab5];
tp_all = [tp_ps(GRD.ID);tp_ps(lid);tp_ps(uid);tp_ps(sid);tp_ps(cid)];
tb_all = [tb_ps(GRD.ID);tb_ps(lid);tb_ps(uid);tb_ps(sid);tb_ps(cid)];
det_all = [det_ps(GRD.ID);det_ps(lid);det_ps(uid);det_ps(sid);det_ps(cid)];
mz_all = [zm_ps(GRD.ID);zm_ps(lid);zm_ps(uid);zm_ps(sid);zm_ps(cid)];
mhp_all = [hpmz_ps(GRD.ID);hpmz_ps(lid);hpmz_ps(uid);hpmz_ps(sid);hpmz_ps(cid)];
lz_all = [zl_ps(GRD.ID);zl_ps(lid);zl_ps(uid);zl_ps(sid);zl_ps(cid)];
lhp_all = [hplz_ps(GRD.ID);hplz_ps(lid);hplz_ps(uid);hplz_ps(sid);hplz_ps(cid)];

figure(11)
subplot(3,3,1)
boxplot(tp_all,blab,...
    'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('TP')

subplot(3,3,2)
boxplot(tb_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('TB')

subplot(3,3,3)
boxplot(det_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('Det')

subplot(3,3,4)
boxplot(mz_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('MZ')

subplot(3,3,5)
boxplot(mhp_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('MZhp')

subplot(3,3,7)
boxplot(lz_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('LZ')

subplot(3,3,8)
boxplot(lhp_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('LZhp')
print('-dpng',[ppath 'CORE_boxplot_ps_forcing_biome.png'])

%%
sf_all = [sf_ps(GRD.ID);sf_ps(lid);sf_ps(uid);sf_ps(sid);sf_ps(cid)];
sp_all = [sp_ps(GRD.ID);sp_ps(lid);sp_ps(uid);sp_ps(sid);sp_ps(cid)];
sd_all = [sd_ps(GRD.ID);sd_ps(lid);sd_ps(uid);sd_ps(sid);sd_ps(cid)];
mf_all = [mf_ps(GRD.ID);mf_ps(lid);mf_ps(uid);mf_ps(sid);mf_ps(cid)];
mp_all = [mp_ps(GRD.ID);mp_ps(lid);mp_ps(uid);mp_ps(sid);mp_ps(cid)];
md_all = [md_ps(GRD.ID);md_ps(lid);md_ps(uid);md_ps(sid);md_ps(cid)];
lp_all = [lp_ps(GRD.ID);lp_ps(lid);lp_ps(uid);lp_ps(sid);lp_ps(cid)];
ld_all = [ld_ps(GRD.ID);ld_ps(lid);ld_ps(uid);ld_ps(sid);ld_ps(cid)];
F_all = [F_ps(GRD.ID);F_ps(lid);F_ps(uid);F_ps(sid);F_ps(cid)];
P_all = [P_ps(GRD.ID);P_ps(lid);P_ps(uid);P_ps(sid);P_ps(cid)];
D_all = [D_ps(GRD.ID);D_ps(lid);D_ps(uid);D_ps(sid);D_ps(cid)];
S_all = [S_ps(GRD.ID);S_ps(lid);S_ps(uid);S_ps(sid);S_ps(cid)];
M_all = [M_ps(GRD.ID);M_ps(lid);M_ps(uid);M_ps(sid);M_ps(cid)];
L_all = [L_ps(GRD.ID);L_ps(lid);L_ps(uid);L_ps(sid);L_ps(cid)];
All_all = [All_ps(GRD.ID);All_ps(lid);All_ps(uid);All_ps(sid);All_ps(cid)];

figure(12)
subplot(3,3,1)
boxplot(sf_all,blab,...
    'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('SF')

subplot(3,3,2)
boxplot(sp_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('SP')

subplot(3,3,3)
boxplot(sd_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('SD')

subplot(3,3,4)
boxplot(mf_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('MF')

subplot(3,3,5)
boxplot(mp_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('MP')

subplot(3,3,6)
boxplot(md_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('MD')

subplot(3,3,8)
boxplot(lp_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('LP')

subplot(3,3,9)
boxplot(ld_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('LD')
print('-dpng',[ppath 'CORE_boxplot_ps_stages_biome.png'])

%%
figure(13)
subplot(3,3,1)
boxplot(F_all,blab,...
    'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('F')

subplot(3,3,2)
boxplot(P_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('P')

subplot(3,3,3)
boxplot(D_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('D')

subplot(3,3,4)
boxplot(S_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('S')

subplot(3,3,5)
boxplot(M_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('M')

subplot(3,3,6)
boxplot(L_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('L')

subplot(3,3,8)
boxplot(All_all,blab,'BoxStyle','outline','Colors','k',...
    'OutlierSize',4,'Symbol','k.')
title('All')

print('-dpng',[ppath 'CORE_boxplot_ps_groups_biome.png'])
