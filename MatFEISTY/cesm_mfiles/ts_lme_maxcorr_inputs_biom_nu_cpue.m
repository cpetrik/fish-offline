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
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/',cfile,'/corrs'];

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

%% Fish data
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

% Anoms with linear trend removed
load([ypath 'FishMIP_Phase3a_LME_CPUE_1961-2010_ann_mean_anoms.mat'])

%Corr of forcing 
tanom = {'TP','TB','Det','ZmLoss','Biom','Prod','Rec','SST','chl'};
cnam = {'corr','p','lag','idriver','driver'};
load([spath,'LMEs_corr_cpue_sat_driver_feisty_maxcorrs.mat']);

%% subset effort years
fyr = 1948:2015;
eyr = 1961:2010;
[yr,fid] = intersect(fyr,eyr);

adety  = adety(lid,fid);
atb    = atb(lid,fid);
atp    = atp(lid,fid);
azlosy = azlosy(lid,fid);
aba    = aba(lid,fid);
ana    = ana(lid,fid);
aall    = aall(lid,:);

yst = 1;
yen = length(eyr);

%% divide by 2 std
% adety  = adety ./ repmat(2*lme_dety_stda,1,50);
% atb    = atb ./ repmat(2*lme_tb_stda,1,50);
% atp    = atp ./ repmat(2*lme_tp_stda,1,50);
% azlosy = azlosy ./ repmat(2*lme_mzly_stda,1,50);
% Ball   = Ball ./ repmat(2*Bstd,1,50);
% Pall   = Pall ./ repmat(2*Pstd,1,50);
% Call   = Call ./ repmat(2*Cstd,1,50);
% Mall   = Mall ./ repmat(2*Mstd,1,50); 

%%
% LMEs - find for each driver
tpid = find(LAtab(:,4)==1);
tbid = find(LAtab(:,4)==2);
did = find(LAtab(:,4)==3);
zid = find(LAtab(:,4)==4);
bid = find(LAtab(:,4)==5);
pid = find(LAtab(:,4)==6);

%lname = {'CHK','EBS','GAK','HI','CCE','GMX','SE','NE','AI'};

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

% Prod
figure(1)
clf
tiledlayout(2,2, 'TileSpacing', 'compact')

for i=1:length(pid)

    lme = pid(i);
    t = LAtab(lme,3);

    drive = ana(lme,yst:yen-t);
    cpue = aall(lme,yst+t:yen);


    nexttile % TP-biom
    %driver
    yyaxis left
    plot(eyr(yst:yen-t),drive,'color',mcol(6,:)); hold on;
    xlim([eyr(1) eyr(end)])
    ylabel('Prod')
    %cpue
    yyaxis right
    plot(eyr(yst:yen-t),cpue,'-k'); hold on;
    xlim([eyr(1) eyr(end)])
    ylabel('CPUE')
    title(['LME ' num2str(lme)])

end
print('-dpng',[ppath 'ts_LMEs_maxcorr_nu_cpue.png'])


%% Biom
figure(2)
tiledlayout(2,2, 'TileSpacing', 'compact')

for i=1:length(bid)

    lme = bid(i);
    t = LAtab(lme,3);

    drive = aba(lme,yst:yen-t);
    cpue = aall(lme,yst+t:yen);


    nexttile % TP-biom
    %driver
    yyaxis left
    plot(eyr(yst:yen-t),drive,'color',mcol(5,:)); hold on;
    xlim([eyr(1) eyr(end)])
    ylabel('Biom')
    %cpue
    yyaxis right
    plot(eyr(yst:yen-t),cpue,'-k'); hold on;
    xlim([eyr(1) eyr(end)])
    ylabel('CPUE')
    title(['LME ' num2str(lme)])

end
print('-dpng',[ppath 'ts_LMEs_maxcorr_biom_cpue.png'])

%% TB
figure(3)
tiledlayout(3,3, 'TileSpacing', 'compact')

for i=1:length(tbid)

    lme = tbid(i);
    t = LAtab(lme,3);

    drive = ana(lme,yst:yen-t);
    cpue = aall(lme,yst+t:yen);


    nexttile % 
    %driver
    yyaxis left
    plot(eyr(yst:yen-t),drive,'color',mcol(2,:)); hold on;
    xlim([eyr(1) eyr(end)])
    ylabel('TB')
    %cpue
    yyaxis right
    plot(eyr(yst:yen-t),cpue,'-k'); hold on;
    xlim([eyr(1) eyr(end)])
    ylabel('CPUE')
    title(['LME ' num2str(lme)])

end
print('-dpng',[ppath 'ts_LMEs_maxcorr_TB_cpue.png'])

%% TP, Det, ZL
figure(4)
tiledlayout(2,2, 'TileSpacing', 'compact')
% TP
lme = tpid;
t = LAtab(lme,3);
drive = ana(lme,yst:yen-t);
cpue = aall(lme,yst+t:yen);

nexttile % TP
%driver
yyaxis left
plot(eyr(yst:yen-t),drive,'color',mcol(1,:)); hold on;
xlim([eyr(1) eyr(end)])
ylabel('Prod')
%cpue
yyaxis right
plot(eyr(yst:yen-t),cpue,'-k'); hold on;
xlim([eyr(1) eyr(end)])
ylabel('CPUE')
title(['LME ' num2str(lme)])

% Det
lme = did;
t = LAtab(lme,3);
drive = ana(lme,yst:yen-t);
cpue = aall(lme,yst+t:yen);

nexttile % TP-biom
%driver
yyaxis left
plot(eyr(yst:yen-t),drive,'color',mcol(3,:)); hold on;
xlim([eyr(1) eyr(end)])
ylabel('Prod')
%cpue
yyaxis right
plot(eyr(yst:yen-t),cpue,'-k'); hold on;
xlim([eyr(1) eyr(end)])
ylabel('CPUE')
title(['LME ' num2str(lme)])


% ZL
lme = zid;
t = LAtab(lme,3);
drive = ana(lme,yst:yen-t);
cpue = aall(lme,yst+t:yen);

nexttile % 
%driver
yyaxis left
plot(eyr(yst:yen-t),drive,'color',mcol(4,:)); hold on;
xlim([eyr(1) eyr(end)])
ylabel('Prod')
%cpue
yyaxis right
plot(eyr(yst:yen-t),cpue,'-k'); hold on;
xlim([eyr(1) eyr(end)])
ylabel('CPUE')
title(['LME ' num2str(lme)])

print('-dpng',[ppath 'ts_LMEs_maxcorr_TP_Det_ZL_cpue.png'])


