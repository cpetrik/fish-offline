% Plot ts of LME CPUE against most significant driver
% that is not satellite SST or chl
% CESM FOSI

clear
close all

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',cfile,'/corrs/'];
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

mod = 'v15_All_fish03';
mod2 = 'v15_obsfish';

%% Sat
load([fpath 'lme_satellite_sst_chl_ann_mean_anoms.mat']);

%% FOSI input forcing
% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'adety','atb','atp','azlosy');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% FEISTY outputs
% Constant effort
%Biomass
load([fpath 'FEISTY_FOSI_',mod,'_lme_ann_mean_anoms.mat'],...
    'aa','ad','af','ap');

Caba = aa;
Cabd = ad;
Cabf = af;
Cabp = ap;

clear aa ad af ap

%Nu
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms.mat'],...
    'aa','ad','af','ap');

Cana = aa;
Cand = ad;
Canf = af;
Canp = ap;

clear aa ad af ap

%% Obs effort
%Biomass
load([fpath 'FEISTY_FOSI_',mod2,'_lme_ann_mean_anoms.mat'],...
    'aa','ad','af','ap');

Oaba = aa;
Oabd = ad;
Oabf = af;
Oabp = ap;

clear aa ad af ap

%Nu
load([fpath 'FEISTY_FOSI_',mod2,'_lme_nu_ann_mean_anoms.mat'],...
    'aa','ad','af','ap');

Oana = aa;
Oand = ad;
Oanf = af;
Oanp = ap;

clear aa ad af ap

%% Fishing data
% Anoms with linear trend removed
load([ypath 'FishMIP_Phase3a_LME_CPUE_1961-2010_ann_mean_anoms.mat'])

%% subset effort years
fyr = 1948:2015;
eyr = 1961:2010;
yr = 1997:2010;
nt = length(yr);

lid=1:66;

[~,cid] = intersect(cyr,yr);
[~,sid] = intersect(tyr,yr);
[~,eid] = intersect(eyr,yr);
[~,fid] = intersect(fyr,yr);

achl   = achl(lid,cid);
asst   = asst(lid,sid);
adety  = adety(lid,fid);
atb    = atb(lid,fid);
atp    = atp(lid,fid);
azlosy = azlosy(lid,fid);
Caba   = Caba(lid,fid);
Cana   = Cana(lid,fid);
Oaba   = Oaba(lid,fid);
Oana   = Oana(lid,fid);
aall   = aall(lid,eid);

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

%% SST
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,asst(i,:),'color',mcol(5,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_LMEs_SST.png'])

%% Chl
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,achl(i,:),'color',mcol(6,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_LMEs_Chl.png'])

%% TP
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,atp(i,:),'color',mcol(1,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_LMEs_TP.png'])

%% TB
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,atb(i,:),'color',mcol(2,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_LMEs_TB.png'])

%% Det
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,adety(i,:),'color',mcol(3,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_LMEs_Det.png'])

%% ZmL
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,azlosy(i,:),'color',mcol(4,:),'LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_LMEs_ZmL.png'])

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
print('-dpng',[ppath 'ts_LMEs_Abio_const.png'])

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
print('-dpng',[ppath 'ts_LMEs_Abio_obsfish.png'])

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
print('-dpng',[ppath 'ts_LMEs_Aprod_const.png'])

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
print('-dpng',[ppath 'ts_LMEs_Aprod_obsfish.png'])

%% CPUE all
figure('Units','inches','Position',[1 3 6.5 8.5]);
i=0;
for m = 1:nrows
    for n = 1:ncols
        i=i+1;
        subplot('position',pos(m,:,n))
        plot(yr,aall(i,:),'k','LineWidth',1); hold on
        title(num2str(i))
        if (i<61)
            set(gca,'XTickLabel',[])
        end
    end
end
print('-dpng',[ppath 'ts_LMEs_CPUE_all.png'])


