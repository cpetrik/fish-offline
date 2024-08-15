% Plot ts of LME CPUE against most significant driver
% that is not satellite SST or chl
% CESM FOSI
% 1997-2010 

clear
close all

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',...
    cfile,'/corrs_cpue/'];
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

mod = 'v15_All_fish03';
mod2 = 'v15_obsfish';

%% Sat
load([fpath 'lme_satellite_sst_chl_ann_mean_anoms_1997_2010_2015.mat'],...
    'achl10','asst10','eyr');

%% FOSI input forcing
% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom_1997_2010_2015.mat'],...
    'adet10','atb10','atp10','azlos10');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% FEISTY outputs
% Constant effort
%Biomass
load([fpath 'FEISTY_FOSI_',mod,'_lme_biom_ann_mean_anoms_1997_2010_2015.mat'],...
    'aba10','abd10','abf10','abp10');

Caba = aba10;
Cabd = abd10;
Cabf = abf10;
Cabp = abp10;

clear aba10 abd10 abf10 abp10

%Nu
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms_1997_2010_2015.mat'],...
    'ana10','and10','anf10','anp10');

Cana = ana10;
Cand = and10;
Canf = anf10;
Canp = anp10;

clear ana10 and10 anf10 anp10

%% Obs effort
%Biomass
load([fpath 'FEISTY_FOSI_',mod2,'_lme_biom_ann_mean_anoms_1997_2010.mat'],...
    'aba10','abd10','abf10','abp10');

Oaba = aba10;
Oabd = abd10;
Oabf = abf10;
Oabp = abp10;

clear aba10 abd10 abf10 abp10

%Nu
load([fpath 'FEISTY_FOSI_',mod2,'_lme_nu_ann_mean_anoms_1997_2010.mat'],...
    'ana10','and10','anf10','anp10');

Oana = ana10;
Oand = and10;
Oanf = anf10;
Oanp = anp10;

clear ana10 and10 anf10 anp10

%% Fishing data
% Anoms with linear trend removed
load([ypath 'FishMIP_Phase3a_LME_CPUE_1997-2010_ann_mean_anoms.mat'])

%% colors
load('paul_tol_cmaps.mat')

%colorblind friendly - subselect & re-order drainbow
ctex = {'TP','TB','Det','ZmLoss','SST','Chl','Biom','Prod'};
% orange, dk blue, grey, lt blue, dk purp, lt purp, red, green
mcol(1,:) = drainbow(12,:)/255; % orange
mcol(2,:) = drainbow(4,:)/255; %dk blue
mcol(3,:) = drainbow(15,:)/255; %grey
mcol(4,:) = drainbow(6,:)/255; %lt blue
mcol(5,:) = drainbow(14,:)/255; %red
mcol(6,:) = drainbow(7,:)/255; %green
mcol(7,:) = drainbow(3,:)/255; %dk purp
mcol(8,:) = drainbow(1,:)/255; %lt purp

colororder(mcol)
close all

%% TS by LME
nrows=11;
ncols=6;
pos = subfigrid(nrows,ncols,[0.05 0.025 0.025 0.05],[0.72 0.75]);
yr = eyr;

%% SST
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,asst10(i,:),'color',mcol(5,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_1997-2010_LMEs_SST.png'])

%% Chl
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,achl10(i,:),'color',mcol(6,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_1997-2010_LMEs_Chl.png'])

%% TP
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,atp10(i,:),'color',mcol(1,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_1997-2010_LMEs_TP.png'])

%% TB
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,atb10(i,:),'color',mcol(2,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_1997-2010_LMEs_TB.png'])

%% Det
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,adet10(i,:),'color',mcol(3,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_1997-2010_LMEs_Det.png'])

%% ZmL
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,azlos10(i,:),'color',mcol(4,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_1997-2010_LMEs_ZmL.png'])

%% Biom
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,Caba(i,:),'color',mcol(7,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_1997-2010_LMEs_Abio_const.png'])

%% Biom
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,Oaba(i,:),'color',mcol(7,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_1997-2010_LMEs_Abio_obsfish.png'])

%% Prod
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,Cana(i,:),'color',mcol(8,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_1997-2010_LMEs_Aprod_const.png'])

%% Prod
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,Oana(i,:),'color',mcol(8,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_1997-2010_LMEs_Aprod_obsfish.png'])

%% CPUE all
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,aa_cpue97(i,:),'k','LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_1997-2010_LMEs_CPUE_all.png'])


