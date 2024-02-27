% Plot ts of LME CPUE against most significant driver
% that is not satellite SST or chl
% CESM FOSI

clear
close all

%% FOSI input forcing

% cpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';
% spath = cpath;
% ypath = cpath;

cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat'],...
    'adety','atb','atp','azlosy','azoo');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% FEISTY outputs
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/',cfile,'/corrs/'];

mod = 'v15_All_fish03';

% Anoms with linear trend removed
%Biomass
load([fpath 'FEISTY_FOSI_',mod,'_lme_ann_mean_anoms.mat'],...
    'aa','ad','af','ap');

aba = aa;
abd = ad;
abf = af;
abp = ap;

clear aa ad af ap

%Nu
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms.mat'],...
    'aa','ad','af','ap');

ana = aa;
and = ad;
anf = af;
anp = ap;

clear aa ad af ap

cnam = {'corr','p','lag','idriver','driver'};

%% subset effort years
fyr = 1948:2015;

yst = 1;
yen = length(fyr);

%% divide by 2 std
lme_tp_stda = std(atp,0,2);
lme_tb_stda = std(atb,0,2);
lme_dety_stda = std(adety,0,2);
lme_mzly_stda = std(azlosy,0,2);
Bstd = std(aba,0,2);
Pstd = std(ana,0,2);

%%
adety  = adety ./ repmat(2*lme_dety_stda,1,68);
atb    = atb ./ repmat(2*lme_tb_stda,1,68);
atp    = atp ./ repmat(2*lme_tp_stda,1,68);
azlosy = azlosy ./ repmat(2*lme_mzly_stda,1,68);
aba   = aba ./ repmat(2*Bstd,1,68);
ana   = ana ./ repmat(2*Pstd,1,68);

%% US LMEs
lid = [54,1:2,65,10,3,5:7]; %ADD 65 = Aleutian Islands
lname = {'CHK','EBS','GAK','AI','HI','CCE','GMX','SE','NE'};

%% colors
load('paul_tol_cmaps.mat')

mcol(1,:) = drainbow(12,:)/255; % orange
mcol(2,:) = drainbow(4,:)/255; %dk blue
mcol(3,:) = drainbow(15,:)/255; %grey
mcol(4,:) = drainbow(6,:)/255; %lt blue
mcol(5,:) = drainbow(3,:)/255; %dk purp
mcol(6,:) = drainbow(1,:)/255; %lt purp

colororder(mcol)

%% plot time series  =================================

f2 = figure('Units','inches','Position',[1 3 14 8]);
for i=1:length(lid)
    clf
    tiledlayout(3,2, 'TileSpacing', 'compact')

    lme = lid(i);
    
    % t = LAtab(lme,3);
    % drive = ana(lme,yst:yen-t);
    % cpue = aall(lme,yst+t:yen);

    nexttile % TP
    plot(fyr,atp(lme,:),'color',mcol(1,:),'LineWidth',2); hold on;
    ylabel('TP')
    xlim([fyr(1) fyr(end)])
    text(2018,1.7,['LME ' num2str(lme) ' - ' lname{i}],...
        'HorizontalAlignment','Center','FontSize',14)
    
    nexttile % Zm
    plot(fyr,azlosy(lme,:),'color',mcol(4,:),'LineWidth',2); hold on;
    ylabel('ZmLoss')
    xlim([fyr(1) fyr(end)])
    
    nexttile % TB
    plot(fyr,atb(lme,:),'color',mcol(2,:),'LineWidth',2); hold on;
    ylabel('TB')
    xlim([fyr(1) fyr(end)])
    
    nexttile % Det
    plot(fyr,adety(lme,:),'color',mcol(3,:),'LineWidth',2); hold on;
    ylabel('Det')
    xlim([fyr(1) fyr(end)])
    
    nexttile % Biom
    plot(fyr,aba(lme,:),'color',mcol(5,:),'LineWidth',2); hold on;
    ylabel('Fish Biom')
    xlim([fyr(1) fyr(end)])
    
    nexttile % Prod
    plot(fyr,ana(lme,:),'color',mcol(6,:),'LineWidth',2); hold on;
    ylabel('Fish Prod')
    xlim([fyr(1) fyr(end)])
    
    print('-dpng',[ppath 'ts_',lname{i},'_drivers_fish_subplot.png'])

end


%% plot time series  =================================

figure(2);
for i=1:length(lid)
    clf
    lme = lid(i);
    
    plot(fyr,atp(lme,:),'color',mcol(1,:),'LineWidth',2); hold on;
    plot(fyr,azlosy(lme,:),'color',mcol(4,:),'LineWidth',2); hold on;
    plot(fyr,atb(lme,:),'color',mcol(2,:),'LineWidth',2); hold on;
    plot(fyr,adety(lme,:),'color',mcol(3,:),'LineWidth',2); hold on;
    plot(fyr,aba(lme,:),'color',mcol(5,:),'LineWidth',2); hold on;
    plot(fyr,ana(lme,:),'color',mcol(6,:),'LineWidth',2); hold on;
    xlim([fyr(1) fyr(end)])
    title(['LME ' num2str(lme) ' - ' lname{i}],'FontSize',14)
   
    print('-dpng',[ppath 'ts_',lname{i},'_drivers_fish_alltogether.png'])

end


