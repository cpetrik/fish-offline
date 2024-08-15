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

%% divide by 2 std to see together
lme_sst_stda = std(asst,0,2);
lme_chl_stda = std(achl,0,2);
lme_tp_stda = std(atp,0,2);
lme_tb_stda = std(atb,0,2);
lme_dety_stda = std(adety,0,2);
lme_mzly_stda = std(azlosy,0,2);
CBstd = std(Caba,0,2);
CPstd = std(Cana,0,2);
OBstd = std(Oaba,0,2);
OPstd = std(Oana,0,2);
Estd = std(aall,0,2);

%%
nt = length(yr);
asst    = asst ./ repmat(2*lme_sst_stda,1,nt);
achl    = achl ./ repmat(2*lme_chl_stda,1,nt);
adety  = adety ./ repmat(2*lme_dety_stda,1,nt);
atb    = atb ./ repmat(2*lme_tb_stda,1,nt);
atp    = atp ./ repmat(2*lme_tp_stda,1,nt);
azlosy = azlosy ./ repmat(2*lme_mzly_stda,1,nt);
Caba   = Caba ./ repmat(2*CBstd,1,nt);
Cana   = Cana ./ repmat(2*CPstd,1,nt);
Oaba   = Oaba ./ repmat(2*OBstd,1,nt);
Oana   = Oana ./ repmat(2*OPstd,1,nt);
aall   = aall ./ repmat(2*Estd,1,nt);

%% plot time series  =================================
yst = 1;
yen = nt;

%% 1  - Biom
lme = 1; %1 = sst(4), det(0), Caba(0), Oaba(0)
sy = yst:(yen-4);
sye = (yst+4):(yen);
cy = yst:yen;

f1 = figure('Units','inches','Position',[1 3 6.5 8]);
subplot(4,2,1)
plot(sye,asst(lme,sy),'color',mcol(5,:),'LineWidth',2); hold on;
plot(cy,adety(lme,cy),'color',mcol(3,:),'LineWidth',2); hold on;
plot(cy,Caba(lme,cy),'color',mcol(7,:),'LineWidth',2); hold on;
plot(cy,Oaba(lme,cy),'-.','color',mcol(7,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('SST','Det','C-Bio','O-Bio','CPUE')
%legend('location','westoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',14)
   
%%  2 - Biom
lme = 2; %2 = chl(4), chl(4), Caba(3), Oaba(0)
sy = yst:(yen-4);
sye = (yst+4):yen;
cy = yst:(yen-3);
cye = (yst+3):yen;
oy = yst:(yen);

subplot(4,2,2)
plot(sye,achl(lme,sy),'color',mcol(6,:),'LineWidth',2); hold on;
plot(cye,Caba(lme,cy),'color',mcol(7,:),'LineWidth',2); hold on;
plot(oy,Oaba(lme,oy),'-.','color',mcol(7,:),'LineWidth',2); hold on;
plot(oy,aall(lme,oy),'k','LineWidth',2); hold on;
%legend('Chl','C-Bio','O-Bio','CPUE')
%legend('location','eastoutside')
xlim([oy(1) oy(end)])
set(gca,'XTick',oy(1):4:oy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%%  14 - Biom
lme = 14; %14 = sst(3), tp(4), Caba(0), Oaba(1)
sy = yst:(yen-3);
sye = (yst+3):yen;
ty = yst:(yen-4);
tye = (yst+4):yen;
cy = (yst):yen;
oy = yst:(yen-1);
oye = (yst+1):yen;

subplot(4,2,3)
plot(sye,-1*asst(lme,sy),'color',mcol(5,:),'LineWidth',2); hold on;
plot(tye,-1*atp(lme,ty),'color',mcol(1,:),'LineWidth',2); hold on;
plot(cy,Caba(lme,cy),'color',mcol(7,:),'LineWidth',2); hold on;
plot(oye,Oaba(lme,oy),'-.','color',mcol(7,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('SST','TP','C-Bio','O-Bio','CPUE')
%legend('location','westoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%% 40  - Biom
lme = 40; %40 = chl(4), ZmL(1), ZmL(1), Oaba(3)
sy = yst:(yen-4);
sye = (yst+4):yen;
zy = yst:(yen-1);
zye = (yst+1):yen;
cy = (yst):yen;
oy = yst:(yen-3);
oye = (yst+3):yen;

subplot(4,2,4)
plot(sye,achl(lme,sy),'color',mcol(6,:),'LineWidth',2); hold on;
plot(tye,azlosy(lme,ty),'color',mcol(4,:),'LineWidth',2); hold on;
plot(oye,Oaba(lme,oy),'-.','color',mcol(7,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('Chl','ZmL','O-Bio','CPUE')
%legend('location','eastoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%%  44 - Biom
lme = 44; %44 = chl(3), ZmL(4), ZmL(4), Oaba(3)
sy = yst:(yen-3);
sye = (yst+3):yen;
zy = yst:(yen-4);
zye = (yst+4):yen;
cy = (yst):yen;
oy = yst:(yen-3);
oye = (yst+3):yen;

subplot(4,2,5)
plot(sye,achl(lme,sy),'color',mcol(6,:),'LineWidth',2); hold on;
plot(tye,azlosy(lme,ty),'color',mcol(4,:),'LineWidth',2); hold on;
plot(oye,Oaba(lme,oy),'-.','color',mcol(7,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('Chl','ZmL','O-Bio','CPUE')
%legend('location','westoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%% 50  - Biom
lme = 50; %50 = chl(0), chl(0), chl(0), Oaba(1)
cy = (yst):yen;
oy = yst:(yen-1);
oye = (yst+1):yen;

subplot(4,2,6)
plot(cy,achl(lme,cy),'color',mcol(6,:),'LineWidth',2); hold on;
plot(oye,Oaba(lme,oy),'-.','color',mcol(7,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('Chl','O-Bio','CPUE')
%legend('location','eastoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)


%% 51  - Biom
lme = 51; %51 = sst(3), det(4), Caba(2), Oaba(2)
sy = yst:(yen-3);
sye = (yst+3):yen;
zy = yst:(yen-4);
zye = (yst+4):yen;
cy = (yst):yen;
oy = yst:(yen-2);
oye = (yst+2):yen;

subplot(4,2,7)
plot(sye,asst(lme,sy),'color',mcol(5,:),'LineWidth',2); hold on;
plot(zye,adety(lme,zy),'color',mcol(3,:),'LineWidth',2); hold on;
plot(oye,Caba(lme,oy),'-','color',mcol(7,:),'LineWidth',2); hold on;
plot(oye,Oaba(lme,oy),'-.','color',mcol(7,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('SST','Det','C-Bio','O-Bio','CPUE')
%legend('location','westoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%% 61  - Biom
lme = 61; %61 = chl(4), ZmL(0), ZmL(0), Oaba(2)
sy = yst:(yen-4);
sye = (yst+4):yen;
cy = (yst):yen;
oy = yst:(yen-2);
oye = (yst+2):yen;

subplot(4,2,8)
plot(sye,achl(lme,sy),'color',mcol(6,:),'LineWidth',2); hold on;
plot(cy,azlosy(lme,cy),'color',mcol(4,:),'LineWidth',2); hold on;
plot(oye,Oaba(lme,oy),'-.','color',mcol(7,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('Chl','ZmL','O-Bio','CPUE')
%legend('location','eastoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

print('-dpng',[ppath 'ts_LMEs_maxcorr_Biom_cpue.png'])

%%  Last biom, prod, ZmL, Det
% 64 - Biom
lme = 64; %64 = sst(2), tp(3), tp(3), Oaba(0)
sy = yst:(yen-2);
sye = (yst+2):yen;
ty = yst:(yen-3);
tye = (yst+3):yen;
cy = (yst):yen;

f2 = figure('Units','inches','Position',[1 4 6.5 8]);
subplot(3,2,1)
plot(sye,asst(lme,sy),'color',mcol(5,:),'LineWidth',2); hold on;
plot(tye,-1*atp(lme,ty),'color',mcol(1,:),'LineWidth',2); hold on;
plot(cy,Oaba(lme,cy),'-.','color',mcol(7,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('SST','TP','O-Bio','CPUE')
%legend('location','westoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%% 48 - Prod
lme = 48; %48 = chl(4), ZmL(3), Cana(2), Oana(2)
sy = yst:(yen-4);
sye = (yst+4):yen;
ty = yst:(yen-3);
tye = (yst+3):yen;
cy = yst:(yen-2);
cye = (yst+2):yen;
ey = (yst):yen;

subplot(3,2,2)
plot(sye,achl(lme,sy),'color',mcol(6,:),'LineWidth',2); hold on;
plot(tye,azlosy(lme,ty),'color',mcol(4,:),'LineWidth',2); hold on;
plot(cye,Cana(lme,cy),'-','color',mcol(8,:),'LineWidth',2); hold on;
plot(cye,Oana(lme,cy),'-.','color',mcol(8,:),'LineWidth',2); hold on;
plot(ey,aall(lme,ey),'k','LineWidth',2); hold on;
%legend('Chl','ZmL','C-Prod','O-Prod','CPUE')
%legend('location','eastoutside')
xlim([ey(1) ey(end)])
set(gca,'XTick',ey(1):4:ey(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%% 8 - ZmL
lme = 8; %8 = chl(2), ZmL(3)
sy = yst:(yen-2);
sye = (yst+2):yen;
ty = yst:(yen-3);
tye = (yst+3):yen;
cy = (yst):yen;

subplot(3,2,3)
plot(sye,achl(lme,sy),'color',mcol(6,:),'LineWidth',2); hold on;
plot(tye,azlosy(lme,ty),'color',mcol(4,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('Chl','ZmL','CPUE')
%legend('location','westoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%% 26 - ZmL
lme = 26; %26 = chl(3), ZmL(1)
sy = yst:(yen-3);
sye = (yst+3):yen;
ty = yst:(yen-1);
tye = (yst+1):yen;
cy = (yst):yen;

subplot(3,2,4)
plot(sye,achl(lme,sy),'color',mcol(6,:),'LineWidth',2); hold on;
plot(tye,azlosy(lme,ty),'color',mcol(4,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('Chl','ZmL','CPUE')
%legend('location','eastoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%% 57 - ZmL
lme = 57; %57 = sst(1), ZmL(1)
sy = yst:(yen-1);
sye = (yst+1):yen;
cy = (yst):yen;

subplot(3,2,5)
plot(sye,asst(lme,sy),'color',mcol(5,:),'LineWidth',2); hold on;
plot(sye,azlosy(lme,sy),'color',mcol(4,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('SST','ZmL','CPUE')
%legend('location','westoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%% 11 - Det
lme = 11; %11 = sst(0), Det(0)
cy = (yst):yen;

subplot(3,2,6)
plot(cy,-1*asst(lme,cy),'color',mcol(5,:),'LineWidth',2); hold on;
plot(cy,adety(lme,cy),'color',mcol(3,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('SST','Det','CPUE')
%legend('location','eastoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

print('-dpng',[ppath 'ts_LMEs_maxcorr_Prod_Det_ZmL_cpue.png'])


%%  TP & TB
% 7 - TP
lme = 7; %7 = sst(3), tp(0)
sy = yst:(yen-3);
sye = (yst+3):yen;
cy = (yst):yen;

f3 = figure('Units','inches','Position',[1 5 6.5 8]);
subplot(3,2,1)
plot(sye,-1*asst(lme,sy),'color',mcol(5,:),'LineWidth',2); hold on;
plot(cy,-1*atp(lme,cy),'color',mcol(1,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('SST','TP','CPUE')
%legend('location','westoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%% 3 - TB
lme = 3; %3 = sst(3), TB(3)
sy = yst:(yen-3);
sye = (yst+3):yen;
cy = (yst):yen;

subplot(3,2,2)
plot(sye,asst(lme,sy),'color',mcol(5,:),'LineWidth',2); hold on;
plot(sye,atb(lme,sy),'color',mcol(2,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('SST','TB','CPUE')
%legend('location','eastoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%% 25 - TB
lme = 25; %25 = sst(2), TB(4)
sy = yst:(yen-2);
sye = (yst+2):yen;
ty = yst:(yen-4);
tye = (yst+4):yen;
cy = (yst):yen;

subplot(3,2,3)
plot(sye,asst(lme,sy),'color',mcol(5,:),'LineWidth',2); hold on;
plot(tye,-1*atb(lme,ty),'color',mcol(2,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('SST','TB','CPUE')
%legend('location','westoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%% 38 - TB
lme = 38; %38 = chl(4), TB(2)
sy = yst:(yen-4);
sye = (yst+4):yen;
ty = yst:(yen-2);
tye = (yst+2):yen;
cy = (yst):yen;

subplot(3,2,4)
plot(sye,achl(lme,sy),'color',mcol(6,:),'LineWidth',2); hold on;
plot(tye,-1*atb(lme,ty),'color',mcol(2,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('Chl','TB','CPUE')
%legend('location','eastoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

%% 65 - TB
lme = 65; %65 = chl(2), TB(1)
sy = yst:(yen-2);
sye = (yst+2):yen;
ty = yst:(yen-1);
tye = (yst+1):yen;
cy = (yst):yen;

subplot(3,2,5)
plot(sye,achl(lme,sy),'color',mcol(6,:),'LineWidth',2); hold on;
plot(tye,-1*atb(lme,ty),'color',mcol(2,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
%legend('Chl','TB','CPUE')
%legend('location','westoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

print('-dpng',[ppath 'ts_LMEs_maxcorr_TP_TB_cpue.png'])

%% 3 - all
lme = 3;
cy = (yst):yen;

figure(4)
plot(cy,asst(lme,cy),'color',mcol(5,:),'LineWidth',2); hold on;
plot(cy,achl(lme,cy),'color',mcol(6,:),'LineWidth',2); hold on;
plot(cy,atp(lme,cy),'color',mcol(1,:),'LineWidth',2); hold on;
plot(cy,atb(lme,cy),'color',mcol(2,:),'LineWidth',2); hold on;
plot(cy,adety(lme,cy),'color',mcol(3,:),'LineWidth',2); hold on;
plot(cy,azlosy(lme,cy),'color',mcol(4,:),'LineWidth',2); hold on;
plot(cy,Caba(lme,cy),'-','color',mcol(7,:),'LineWidth',2); hold on;
plot(cy,Oaba(lme,cy),'-.','color',mcol(7,:),'LineWidth',2); hold on;
plot(cy,Cana(lme,cy),'-','color',mcol(8,:),'LineWidth',2); hold on;
plot(cy,Oana(lme,cy),'-.','color',mcol(8,:),'LineWidth',2); hold on;
plot(cy,aall(lme,cy),'k','LineWidth',2); hold on;
legend('TP','TB','Det','ZmLoss','SST','Chl','C-Bio','O-Bio','C-Prod','O-Prod','CPUE')
legend('location','eastoutside')
xlim([cy(1) cy(end)])
set(gca,'XTick',cy(1):4:cy(end),'XTickLabel',yr(1):4:yr(end))
title(['LME ' num2str(lme)],'FontSize',12)

print('-dpng',[ppath 'ts_LMEs_maxcorr_cpue_legend.png'])


