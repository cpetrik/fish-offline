% Plot ts of LME means of FEISTY inputs & outputs w/ climate anoms
% CESM FOSI

clear 
close all

%% Climate indices
apath = '/Users/cpetrik/Dropbox/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
load([apath 'Climate_anomalies_annual_means.mat']);

tanom = canom;
clear canom

%% FOSI input forcing
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_interann_mean_forcings_anom.mat']);
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat']);

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
ppath = [pp cfile '/corrs/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

sims = {'v15_All_fish03_';'v15_climatol_';'v15_varFood_';'v15_varTemp_'};
mod = sims{1};

load([fpath 'FEISTY_FOSI_',mod,'lme_ann_mean_anoms.mat']) % Anoms with linear trend removed

%% 
% LMEs
lid = [54,1:2,10,3,5:7,65]; %ADD 65 = Aleutian Islands
lname = {'CHK','EBS','GAK','HI','CCE','GMX','SE','NE','AI'};
% forcing inputs
iname = {'TP','TB','Det','Zmeso','ZmLoss'};
% FEISTY outputs grouped
gname = {'S','M','L','F','P','D','A','B'};

%%
BSanom(:,1) = [1948:2015]';
AKanom(:,1) = [1948:2015]';
AIanom(:,1) = [1948:2015]';
HIanom(:,1) = [1948:2015]';
CCanom(:,1) = [1948:2015]';
NEanom(:,1) = [1948:2015]';
SEanom(:,1) = [1948:2015]';

BSanom(:,2) = atp(lid(2),:);
BSanom(:,3) = atb(lid(2),:);
BSanom(:,4) = adet(lid(2),:);
BSanom(:,5) = azoo(lid(2),:);
BSanom(:,6) = azlos(lid(2),:);

AKanom(:,2) = atp(lid(3),:);
AKanom(:,3) = atb(lid(3),:);
AKanom(:,4) = adet(lid(3),:);
AKanom(:,5) = azoo(lid(3),:);
AKanom(:,6) = azlos(lid(3),:);

HIanom(:,2) = atp(lid(4),:);
HIanom(:,3) = atb(lid(4),:);
HIanom(:,4) = adet(lid(4),:);
HIanom(:,5) = azoo(lid(4),:);
HIanom(:,6) = azlos(lid(4),:);

CCanom(:,2) = atp(lid(5),:);
CCanom(:,3) = atb(lid(5),:);
CCanom(:,4) = adet(lid(5),:);
CCanom(:,5) = azoo(lid(5),:);
CCanom(:,6) = azlos(lid(5),:);

BSanom(:,7)  = as(lid(2),:);
BSanom(:,8)  = am(lid(2),:);
BSanom(:,9)  = al(lid(2),:);
BSanom(:,10) = af(lid(2),:);
BSanom(:,11) = ap(lid(2),:);
BSanom(:,12) = ad(lid(2),:);
BSanom(:,13) = aa(lid(2),:);
BSanom(:,14) = ab(lid(2),:);

AKanom(:,7)  = as(lid(3),:);
AKanom(:,8)  = am(lid(3),:);
AKanom(:,9)  = al(lid(3),:);
AKanom(:,10) = af(lid(3),:);
AKanom(:,11) = ap(lid(3),:);
AKanom(:,12) = ad(lid(3),:);
AKanom(:,13) = aa(lid(3),:);
AKanom(:,14) = ab(lid(3),:);

HIanom(:,7)  = as(lid(4),:);
HIanom(:,8)  = am(lid(4),:);
HIanom(:,9)  = al(lid(4),:);
HIanom(:,10) = af(lid(4),:);
HIanom(:,11) = ap(lid(4),:);
HIanom(:,12) = ad(lid(4),:);
HIanom(:,13) = aa(lid(4),:);
HIanom(:,14) = ab(lid(4),:);

CCanom(:,7)  = as(lid(5),:);
CCanom(:,8)  = am(lid(5),:);
CCanom(:,9)  = al(lid(5),:);
CCanom(:,10) = af(lid(5),:);
CCanom(:,11) = ap(lid(5),:);
CCanom(:,12) = ad(lid(5),:);
CCanom(:,13) = aa(lid(5),:);
CCanom(:,14) = ab(lid(5),:);

SEanom(:,2) = atp(lid(7),:);
SEanom(:,3) = atb(lid(7),:);
SEanom(:,4) = adet(lid(7),:);
SEanom(:,5) = azoo(lid(7),:);
SEanom(:,6) = azlos(lid(7),:);
SEanom(:,7)  = as(lid(7),:);
SEanom(:,8)  = am(lid(7),:);
SEanom(:,9)  = al(lid(7),:);
SEanom(:,10) = af(lid(7),:);
SEanom(:,11) = ap(lid(7),:);
SEanom(:,12) = ad(lid(7),:);
SEanom(:,13) = aa(lid(7),:);
SEanom(:,14) = ab(lid(7),:);

NEanom(:,2) = atp(lid(8),:);
NEanom(:,3) = atb(lid(8),:);
NEanom(:,4) = adet(lid(8),:);
NEanom(:,5) = azoo(lid(8),:);
NEanom(:,6) = azlos(lid(8),:);
NEanom(:,7)  = as(lid(8),:);
NEanom(:,8)  = am(lid(8),:);
NEanom(:,9)  = al(lid(8),:);
NEanom(:,10) = af(lid(8),:);
NEanom(:,11) = ap(lid(8),:);
NEanom(:,12) = ad(lid(8),:);
NEanom(:,13) = aa(lid(8),:);
NEanom(:,14) = ab(lid(8),:);

AIanom(:,2) = atp(lid(9),:);
AIanom(:,3) = atb(lid(9),:);
AIanom(:,4) = adet(lid(9),:);
AIanom(:,5) = azoo(lid(9),:);
AIanom(:,6) = azlos(lid(9),:);
AIanom(:,7)  = as(lid(9),:);
AIanom(:,8)  = am(lid(9),:);
AIanom(:,9)  = al(lid(9),:);
AIanom(:,10) = af(lid(9),:);
AIanom(:,11) = ap(lid(9),:);
AIanom(:,12) = ad(lid(9),:);
AIanom(:,13) = aa(lid(9),:);
AIanom(:,14) = ab(lid(9),:);

%% plot time series like LeMezo 2016  =================================
%         PUT CLIMATE ON A DIFF Y-AXIS
% PDO
figure(1)
subplot(4,1,1)
%climate = PDO
plot(BSanom(:,1),manom(10,:),'color',0.5*[1 1 1]); hold on;
%TP
plot(BSanom(:,1),BSanom(:,2),'k'); hold on;
%MZ
plot(BSanom(:,1),BSanom(:,5),'b'); hold on;
%All fish
plot(BSanom(:,1),BSanom(:,13),'r'); hold on;
xlim([BSanom(1,1) BSanom(end,1)])
ylabel('EBS')
title('PDO pelagic')

subplot(4,1,2)
%climate = PDO
plot(AKanom(:,1),manom(10,:),'color',0.5*[1 1 1]); hold on;
%TP
plot(AKanom(:,1),AKanom(:,2),'k'); hold on;
%MZ
plot(AKanom(:,1),AKanom(:,5),'b'); hold on;
%All fish
plot(AKanom(:,1),AKanom(:,13),'r'); hold on;
xlim([AKanom(1,1) AKanom(end,1)])
ylabel('GoAK')

subplot(4,1,3)
%climate = PDO
plot(CCanom(:,1),manom(10,:),'color',0.5*[1 1 1]); hold on;
%TP
plot(CCanom(:,1),CCanom(:,2),'k'); hold on;
%MZ
plot(CCanom(:,1),CCanom(:,5),'b'); hold on;
%All fish
plot(CCanom(:,1),CCanom(:,13),'r'); hold on;
xlim([CCanom(1,1) CCanom(end,1)])
ylabel('CCE')

subplot(4,1,4)
%climate = PDO
plot(HIanom(:,1),manom(10,:),'color',0.5*[1 1 1]); hold on;
%TP
plot(HIanom(:,1),HIanom(:,2),'k'); hold on;
%MZ
plot(HIanom(:,1),HIanom(:,5),'b'); hold on;
%All fish
plot(HIanom(:,1),HIanom(:,13),'r'); hold on;
xlim([HIanom(1,1) HIanom(end,1)])
ylabel('HI')
print('-dpng',[ppath 'ts_climate_FOSI_NPac_PDO_',mod,'pel.png'])
    
%
figure(2)
%climate = PDO
plot(BSanom(:,1),manom(10,:),'color',0.5*[1 1 1]); hold on;
%TP
plot(BSanom(:,1),BSanom(:,2),'k'); hold on;
%MZ
plot(BSanom(:,1),BSanom(:,5),'b'); hold on;
%All fish
plot(BSanom(:,1),BSanom(:,13),'r'); hold on;
legend('PDO','Tpel','LTL','HTL')
legend('location','northwest')
print('-dpng',[ppath 'ts_climate_FOSI_NPac_PDO_',mod,'pel_legend.png'])
    
%% bottom
% PDO
figure(3)
subplot(4,1,1)
%climate = PDO
plot(BSanom(:,1),manom(10,:),'color',0.5*[1 1 1]); hold on;
%TP
plot(BSanom(:,1),BSanom(:,3),'k'); hold on;
%MZ
plot(BSanom(:,1),BSanom(:,4),'b'); hold on;
%All fish
plot(BSanom(:,1),BSanom(:,12),'r'); hold on;
xlim([BSanom(1,1) BSanom(end,1)])
ylabel('EBS')
title('PDO benthic')

subplot(4,1,2)
%climate = PDO
plot(AKanom(:,1),manom(10,:),'color',0.5*[1 1 1]); hold on;
%TP
plot(AKanom(:,1),AKanom(:,3),'k'); hold on;
%MZ
plot(AKanom(:,1),AKanom(:,4),'b'); hold on;
%All fish
plot(AKanom(:,1),AKanom(:,12),'r'); hold on;
xlim([AKanom(1,1) AKanom(end,1)])
ylabel('GoAK')

subplot(4,1,3)
%climate = PDO
plot(CCanom(:,1),manom(10,:),'color',0.5*[1 1 1]); hold on;
%TP
plot(CCanom(:,1),CCanom(:,3),'k'); hold on;
%MZ
plot(CCanom(:,1),CCanom(:,4),'b'); hold on;
%All fish
plot(CCanom(:,1),CCanom(:,12),'r'); hold on;
xlim([CCanom(1,1) CCanom(end,1)])
ylabel('CCE')

subplot(4,1,4)
%climate = PDO
plot(HIanom(:,1),manom(10,:),'color',0.5*[1 1 1]); hold on;
%TP
plot(HIanom(:,1),HIanom(:,3),'k'); hold on;
%MZ
plot(HIanom(:,1),HIanom(:,4),'b'); hold on;
%All fish
plot(HIanom(:,1),HIanom(:,12),'r'); hold on;
xlim([HIanom(1,1) HIanom(end,1)])
ylabel('HI')
print('-dpng',[ppath 'ts_climate_FOSI_NPac_PDO_',mod,'btm.png'])
    
%
figure(4)
%climate = PDO
plot(BSanom(:,1),manom(10,:),'color',0.5*[1 1 1]); hold on;
%TP
plot(BSanom(:,1),BSanom(:,3),'k'); hold on;
%MZ
plot(BSanom(:,1),BSanom(:,4),'b'); hold on;
%All fish
plot(BSanom(:,1),BSanom(:,12),'r'); hold on;
legend('PDO','Tbtm','LTL','HTL')
legend('location','northwest')
print('-dpng',[ppath 'ts_climate_FOSI_NPac_PDO_',mod,'btm_legend.png'])
 
%% no climate  ================================================
% NPac & NAtl
figure(5)
subplot(5,1,1)
%TP
plot(BSanom(:,1),BSanom(:,2),'k'); hold on;
%MZ
plot(BSanom(:,1),BSanom(:,5),'b'); hold on;
%All fish
plot(BSanom(:,1),BSanom(:,13),'r'); hold on;
xlim([BSanom(1,1) BSanom(end,1)])
ylabel('EBS')
title('Pelagic')

subplot(5,1,2)
%TP
plot(AKanom(:,1),AKanom(:,2),'k'); hold on;
%MZ
plot(AKanom(:,1),AKanom(:,5),'b'); hold on;
%All fish
plot(AKanom(:,1),AKanom(:,13),'r'); hold on;
xlim([AKanom(1,1) AKanom(end,1)])
ylabel('GoAK')

subplot(5,1,3)
%TP
plot(CCanom(:,1),CCanom(:,2),'k'); hold on;
%MZ
plot(CCanom(:,1),CCanom(:,5),'b'); hold on;
%All fish
plot(CCanom(:,1),CCanom(:,13),'r'); hold on;
xlim([CCanom(1,1) CCanom(end,1)])
ylabel('CCE')

subplot(5,1,4)
%TP
plot(HIanom(:,1),HIanom(:,2),'k'); hold on;
%MZ
plot(HIanom(:,1),HIanom(:,5),'b'); hold on;
%All fish
plot(HIanom(:,1),HIanom(:,13),'r'); hold on;
xlim([HIanom(1,1) HIanom(end,1)])
ylabel('HI')

subplot(5,1,5)
%TP
plot(NEanom(:,1),NEanom(:,2),'k'); hold on;
%MZ
plot(NEanom(:,1),NEanom(:,5),'b'); hold on;
%All fish
plot(NEanom(:,1),NEanom(:,13),'r'); hold on;
xlim([NEanom(1,1) NEanom(end,1)])
ylabel('NE')
print('-dpng',[ppath 'ts_climate_FOSI_',mod,'pel.png'])
    
%
figure(6)
%TP
plot(BSanom(:,1),BSanom(:,2),'k'); hold on;
%MZ
plot(BSanom(:,1),BSanom(:,5),'b'); hold on;
%All fish
plot(BSanom(:,1),BSanom(:,13),'r'); hold on;
legend('Tpel','LTL','HTL')
legend('location','northwest')
print('-dpng',[ppath 'ts_climate_FOSI_NPac_',mod,'pel_legend.png'])
    
%% bottom
figure(7)
subplot(5,1,1)
%TP
plot(BSanom(:,1),BSanom(:,3),'k'); hold on;
%MZ
plot(BSanom(:,1),BSanom(:,4),'b'); hold on;
%All fish
plot(BSanom(:,1),BSanom(:,12),'r'); hold on;
xlim([BSanom(1,1) BSanom(end,1)])
ylabel('EBS')
title('Benthic')

subplot(5,1,2)
%TP
plot(AKanom(:,1),AKanom(:,3),'k'); hold on;
%MZ
plot(AKanom(:,1),AKanom(:,4),'b'); hold on;
%All fish
plot(AKanom(:,1),AKanom(:,12),'r'); hold on;
xlim([AKanom(1,1) AKanom(end,1)])
ylabel('GoAK')

subplot(5,1,3)
%TP
plot(CCanom(:,1),CCanom(:,3),'k'); hold on;
%MZ
plot(CCanom(:,1),CCanom(:,4),'b'); hold on;
%All fish
plot(CCanom(:,1),CCanom(:,12),'r'); hold on;
xlim([CCanom(1,1) CCanom(end,1)])
ylabel('CCE')

subplot(5,1,4)
%TP
plot(HIanom(:,1),HIanom(:,3),'k'); hold on;
%MZ
plot(HIanom(:,1),HIanom(:,4),'b'); hold on;
%All fish
plot(HIanom(:,1),HIanom(:,12),'r'); hold on;
xlim([HIanom(1,1) HIanom(end,1)])
ylabel('HI')

subplot(5,1,5)
%TP
plot(NEanom(:,1),NEanom(:,3),'k'); hold on;
%MZ
plot(NEanom(:,1),NEanom(:,4),'b'); hold on;
%All fish
plot(NEanom(:,1),NEanom(:,12),'r'); hold on;
xlim([NEanom(1,1) NEanom(end,1)])
ylabel('NE')
print('-dpng',[ppath 'ts_climate_FOSI_',mod,'btm.png'])
    
%
figure(8)
%TP
plot(BSanom(:,1),BSanom(:,3),'k'); hold on;
%MZ
plot(BSanom(:,1),BSanom(:,4),'b'); hold on;
%All fish
plot(BSanom(:,1),BSanom(:,12),'r'); hold on;
legend('Tbtm','LTL','HTL')
legend('location','northwest')
print('-dpng',[ppath 'ts_climate_FOSI_',mod,'btm_legend.png'])


