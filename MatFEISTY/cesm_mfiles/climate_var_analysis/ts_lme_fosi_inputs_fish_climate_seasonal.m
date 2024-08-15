% Plot ts of LME means of FEISTY inputs & outputs w/ climate anoms
% CESM FOSI

clear 
close all

%% Climate indices
apath = '/Users/cpetrik/Dropbox/NCAR/MAPP-METF/NCAR3/DPLE_offline/results_dple/climate_indices/';
load([apath 'Climate_anomalies_seasonal_means.mat']);

%% FOSI input forcing
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

% lme means, trend removed, anomaly calc
%load([cpath 'CESM_FOSI_v15_interann_mean_forcings_anom.mat']);
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom.mat']);

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
ppath = [pp cfile '/'];
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
CKanom(:,1) = [1948:2015]';
MXanom(:,1) = [1948:2015]';

CKanom(:,2) = atp(lid(1),:);
CKanom(:,3) = atb(lid(1),:);
CKanom(:,4) = adet(lid(1),:);
CKanom(:,5) = azoo(lid(1),:);
CKanom(:,6) = azlos(lid(1),:);
CKanom(:,7)  = as(lid(1),:);
CKanom(:,8)  = am(lid(1),:);
CKanom(:,9)  = al(lid(1),:);
CKanom(:,10) = af(lid(1),:);
CKanom(:,11) = ap(lid(1),:);
CKanom(:,12) = ad(lid(1),:);
CKanom(:,13) = aa(lid(1),:);
CKanom(:,14) = ab(lid(1),:);

BSanom(:,2) = atp(lid(2),:);
BSanom(:,3) = atb(lid(2),:);
BSanom(:,4) = adet(lid(2),:);
BSanom(:,5) = azoo(lid(2),:);
BSanom(:,6) = azlos(lid(2),:);
BSanom(:,7)  = as(lid(2),:);
BSanom(:,8)  = am(lid(2),:);
BSanom(:,9)  = al(lid(2),:);
BSanom(:,10) = af(lid(2),:);
BSanom(:,11) = ap(lid(2),:);
BSanom(:,12) = ad(lid(2),:);
BSanom(:,13) = aa(lid(2),:);
BSanom(:,14) = ab(lid(2),:);

AKanom(:,2) = atp(lid(3),:);
AKanom(:,3) = atb(lid(3),:);
AKanom(:,4) = adet(lid(3),:);
AKanom(:,5) = azoo(lid(3),:);
AKanom(:,6) = azlos(lid(3),:);
AKanom(:,7)  = as(lid(3),:);
AKanom(:,8)  = am(lid(3),:);
AKanom(:,9)  = al(lid(3),:);
AKanom(:,10) = af(lid(3),:);
AKanom(:,11) = ap(lid(3),:);
AKanom(:,12) = ad(lid(3),:);
AKanom(:,13) = aa(lid(3),:);
AKanom(:,14) = ab(lid(3),:);

HIanom(:,2) = atp(lid(4),:);
HIanom(:,3) = atb(lid(4),:);
HIanom(:,4) = adet(lid(4),:);
HIanom(:,5) = azoo(lid(4),:);
HIanom(:,6) = azlos(lid(4),:);
HIanom(:,7)  = as(lid(4),:);
HIanom(:,8)  = am(lid(4),:);
HIanom(:,9)  = al(lid(4),:);
HIanom(:,10) = af(lid(4),:);
HIanom(:,11) = ap(lid(4),:);
HIanom(:,12) = ad(lid(4),:);
HIanom(:,13) = aa(lid(4),:);
HIanom(:,14) = ab(lid(4),:);

CCanom(:,2) = atp(lid(5),:);
CCanom(:,3) = atb(lid(5),:);
CCanom(:,4) = adet(lid(5),:);
CCanom(:,5) = azoo(lid(5),:);
CCanom(:,6) = azlos(lid(5),:);
CCanom(:,7)  = as(lid(5),:);
CCanom(:,8)  = am(lid(5),:);
CCanom(:,9)  = al(lid(5),:);
CCanom(:,10) = af(lid(5),:);
CCanom(:,11) = ap(lid(5),:);
CCanom(:,12) = ad(lid(5),:);
CCanom(:,13) = aa(lid(5),:);
CCanom(:,14) = ab(lid(5),:);

MXanom(:,2) = atp(lid(6),:);
MXanom(:,3) = atb(lid(6),:);
MXanom(:,4) = adet(lid(6),:);
MXanom(:,5) = azoo(lid(6),:);
MXanom(:,6) = azlos(lid(6),:);
MXanom(:,7)  = as(lid(6),:);
MXanom(:,8)  = am(lid(6),:);
MXanom(:,9)  = al(lid(6),:);
MXanom(:,10) = af(lid(6),:);
MXanom(:,11) = ap(lid(6),:);
MXanom(:,12) = ad(lid(6),:);
MXanom(:,13) = aa(lid(6),:);
MXanom(:,14) = ab(lid(6),:);

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

%% colors
% cm10=[0.5 0.5 0;... %tan/army
%     0 0.7 0;...   %g
%     1 0 1;...     %m
%     1 0 0;...     %r
%     0.5 0 0;...   %maroon
%     0/255 206/255 209/255;... %turq
%     0 0.5 0.75;...   %med blue
%     0 0 0.75;...    %b
%     0.5 0.5 0.5; ...    %med grey
%     0 0 0];...      %black
    
cm10=[0.5 0.5 0.5; ... %med grey
    0 0 0; ...      %black
    0 0 0.75;...    %b
    1 0 0;...       %r
    0 0.5 0.75;...   %med blue
    0.5 0 0;...   %maroon
    0/255 206/255 209/255;... %turq
    1 0 1;...     %m
    0 0.7 0;...   %g
    ];

set(groot,'defaultAxesColorOrder',cm10);

%% plot time series by size  =================================
%WHY ARE FISH AMOUNTS DIFF THAN ANN MEAN FIG?!
% PDO
figure(1)
%CHK
%subplot(4,1,1)
tiledlayout(4,1, 'TileSpacing', 'compact')
nexttile
%climate = PDO
yyaxis left
plot(CKanom(:,1),manom(5,:),'color',0.5*[1 1 1]); hold on;
ylim([-2 2])
ylabel('PDO')
yyaxis right
%Sm
plot(CKanom(:,1),100*CKanom(:,7),'-k'); hold on;
%Md
plot(CKanom(:,1),CKanom(:,8),'-b'); hold on;
%Lg
plot(CKanom(:,1),CKanom(:,9),'-r'); hold on;
ylim([-1.5 1.5])
xlim([CKanom(1,1) CKanom(end,1)])
title('CHK')
ylabel('fish')

nexttile %subplot(4,1,2)
%climate = PDO
yyaxis left
plot(BSanom(:,1),manom(5,:),'color',0.5*[1 1 1]); hold on;
ylim([-2 2])
ylabel('PDO')
yyaxis right
%S
plot(BSanom(:,1),50*BSanom(:,7),'-k'); hold on;
%M
plot(BSanom(:,1),BSanom(:,8),'-b'); hold on;
%L
plot(BSanom(:,1),BSanom(:,9),'-r'); hold on;
ylim([-0.8 0.8])
xlim([BSanom(1,1) BSanom(end,1)])
title('EBS')
ylabel('fish')

% AK
nexttile %subplot(4,1,3)
%climate = PDO
yyaxis left
plot(AKanom(:,1),manom(5,:),'color',0.5*[1 1 1]); hold on;
ylim([-2 2])
ylabel('PDO')
yyaxis right
%S
plot(AKanom(:,1),10*AKanom(:,7),'-k'); hold on;
%M
plot(AKanom(:,1),AKanom(:,8),'-b'); hold on;
%L
plot(AKanom(:,1),AKanom(:,9),'-r'); hold on;
ylim([-0.25 0.25])
xlim([AKanom(1,1) AKanom(end,1)])
title('GAK')
ylabel('fish')

%CCE
nexttile %subplot(4,1,4)
%climate = PDO
yyaxis left
plot(CCanom(:,1),manom(5,:),'color',0.5*[1 1 1]); hold on;
ylim([-2 2])
ylabel('PDO')
yyaxis right
%S
plot(CCanom(:,1),10*CCanom(:,7),'-k'); hold on;
%M
plot(CCanom(:,1),CCanom(:,8),'-b'); hold on;
%L
plot(CCanom(:,1),CCanom(:,9),'-r'); hold on;
ylim([-0.05 0.05])
xlim([CCanom(1,1) CCanom(end,1)])
title('CCE')
ylabel('fish')
%legend({'PDO','Sm','Md','Lg'},'location','southoutside','orientation','horizontal')
lg = legend(nexttile(4),{'PDO','Sm','Md','Lg'});
lg.Location = 'southoutside';
lg.Orientation = 'horizontal';
print('-dpng',[ppath 'ts_climate_seasonal_FOSI_NPac_PDO_',mod,'size.png'])

%% ENSO
%WHY ARE FISH AMOUNTS DIFF THAN ANN MEAN FIG?!
figure(2)
%CHK
%subplot(4,1,1)
tiledlayout(4,1, 'TileSpacing', 'compact')
nexttile
%climate = ENSO
yyaxis left
plot(CKanom(:,1),manom(4,:),'color',0.5*[1 1 1]); hold on;
ylim([-1.5 1.5])
ylabel('ENSO')
yyaxis right
%Sm
plot(CKanom(:,1),100*CKanom(:,7),'-k'); hold on;
%Md
plot(CKanom(:,1),CKanom(:,8),'-b'); hold on;
%Lg
plot(CKanom(:,1),CKanom(:,9),'-r'); hold on;
ylim([-1.5 1.5])
xlim([CKanom(1,1) CKanom(end,1)])
title('CHK')
ylabel('fish')

nexttile %subplot(4,1,2)
%climate = ENSO
yyaxis left
plot(BSanom(:,1),manom(4,:),'color',0.5*[1 1 1]); hold on;
ylim([-1.5 1.5])
ylabel('ENSO')
yyaxis right
%S
plot(BSanom(:,1),50*BSanom(:,7),'-k'); hold on;
%M
plot(BSanom(:,1),BSanom(:,8),'-b'); hold on;
%L
plot(BSanom(:,1),BSanom(:,9),'-r'); hold on;
ylim([-0.8 0.8])
xlim([BSanom(1,1) BSanom(end,1)])
title('EBS')
ylabel('fish')

% AK
nexttile %subplot(4,1,3)
%climate = ENSO
yyaxis left
plot(AKanom(:,1),manom(4,:),'color',0.5*[1 1 1]); hold on;
ylim([-1.5 1.5])
ylabel('ENSO')
yyaxis right
%S
plot(AKanom(:,1),10*AKanom(:,7),'-k'); hold on;
%M
plot(AKanom(:,1),AKanom(:,8),'-b'); hold on;
%L
plot(AKanom(:,1),AKanom(:,9),'-r'); hold on;
ylim([-0.25 0.25])
xlim([AKanom(1,1) AKanom(end,1)])
title('GAK')
ylabel('fish')

%CCE
nexttile %subplot(4,1,4)
%climate = ENSO
yyaxis left
plot(CCanom(:,1),manom(4,:),'color',0.5*[1 1 1]); hold on;
ylim([-1.5 1.5])
ylabel('ENSO')
yyaxis right
%S
plot(CCanom(:,1),10*CCanom(:,7),'-k'); hold on;
%M
plot(CCanom(:,1),CCanom(:,8),'-b'); hold on;
%L
plot(CCanom(:,1),CCanom(:,9),'-r'); hold on;
ylim([-0.05 0.05])
xlim([CCanom(1,1) CCanom(end,1)])
title('CCE')
ylabel('fish')
%legend({'ENSO','Sm','Md','Lg'},'location','southoutside','orientation','horizontal')
lg = legend(nexttile(4),{'ENSO','Sm','Md','Lg'});
lg.Location = 'southoutside';
lg.Orientation = 'horizontal';
print('-dpng',[ppath 'ts_climate_seasonal_FOSI_NPac_ENSO_',mod,'size.png'])  

%% AMO
%WHY ARE FISH AMOUNTS DIFF THAN ANN MEAN FIG?!
figure(3)
%NE
%subplot(3,1,1)
tiledlayout(3,1, 'TileSpacing', 'compact')
nexttile
yyaxis left
plot(NEanom(:,1),manom(1,:),'color',0.5*[1 1 1]); hold on;
ylim([-0.5 0.5])
ylabel('AMO')
yyaxis right
%Sm
plot(NEanom(:,1),50*NEanom(:,7),'-k'); hold on;
%Md
plot(NEanom(:,1),NEanom(:,8),'-b'); hold on;
%Lg
plot(NEanom(:,1),NEanom(:,9),'-r'); hold on;
ylim([-0.2 0.2])
xlim([NEanom(1,1) NEanom(end,1)])
title('NE')
ylabel('fish')

nexttile %subplot(2,1,2)
% SE
yyaxis left
plot(SEanom(:,1),manom(1,:),'color',0.5*[1 1 1]); hold on;
ylim([-0.5 0.5])
ylabel('AMO')
yyaxis right
%S
plot(SEanom(:,1),25*SEanom(:,7),'-k'); hold on;
%M
plot(SEanom(:,1),SEanom(:,8),'-b'); hold on;
%L
plot(SEanom(:,1),SEanom(:,9),'-r'); hold on;
ylim([-0.08 0.08])
xlim([SEanom(1,1) SEanom(end,1)])
title('SE')
ylabel('fish')

% GMX
nexttile %subplot(4,1,3)
%climate = AMO
yyaxis left
plot(MXanom(:,1),manom(1,:),'color',0.5*[1 1 1]); hold on;
ylim([-0.5 0.5])
ylabel('AMO')
yyaxis right
%S
plot(MXanom(:,1),50*MXanom(:,7),'-k'); hold on;
%M
plot(MXanom(:,1),MXanom(:,8),'-b'); hold on;
%L
plot(MXanom(:,1),MXanom(:,9),'-r'); hold on;
ylim([-0.6 0.6])
xlim([MXanom(1,1) MXanom(end,1)])
title('GMX')
ylabel('fish')

%legend({'AMO','Sm','Md','Lg'},'location','southoutside','orientation','horizontal')
lg = legend(nexttile(3),{'AMO','Sm','Md','Lg'});
lg.Location = 'southoutside';
lg.Orientation = 'horizontal';
print('-dpng',[ppath 'ts_climate_seasonal_FOSI_NAtl_AMO_',mod,'size.png'])

%% NAO
%WHY ARE FISH AMOUNTS DIFF THAN ANN MEAN FIG?!
figure(4)
%NE
%subplot(3,1,1)
tiledlayout(3,1, 'TileSpacing', 'compact')
nexttile
yyaxis left
plot(NEanom(:,1),manom(3,:),'color',0.5*[1 1 1]); hold on;
ylim([-3 3])
ylabel('NAO')
yyaxis right
%Sm
plot(NEanom(:,1),50*NEanom(:,7),'-k'); hold on;
%Md
plot(NEanom(:,1),NEanom(:,8),'-b'); hold on;
%Lg
plot(NEanom(:,1),NEanom(:,9),'-r'); hold on;
ylim([-0.2 0.2])
xlim([NEanom(1,1) NEanom(end,1)])
title('NE')
ylabel('fish')

nexttile %subplot(2,1,2)
% SE
yyaxis left
plot(SEanom(:,1),manom(3,:),'color',0.5*[1 1 1]); hold on;
ylim([-3 3])
ylabel('NAO')
yyaxis right
%S
plot(SEanom(:,1),25*SEanom(:,7),'-k'); hold on;
%M
plot(SEanom(:,1),SEanom(:,8),'-b'); hold on;
%L
plot(SEanom(:,1),SEanom(:,9),'-r'); hold on;
ylim([-0.08 0.08])
xlim([SEanom(1,1) SEanom(end,1)])
title('SE')
ylabel('fish')

% GMX
nexttile %subplot(4,1,3)
%climate = NAO
yyaxis left
plot(MXanom(:,1),manom(3,:),'color',0.5*[1 1 1]); hold on;
ylim([-3 3])
ylabel('NAO')
yyaxis right
%S
plot(MXanom(:,1),50*MXanom(:,7),'-k'); hold on;
%M
plot(MXanom(:,1),MXanom(:,8),'-b'); hold on;
%L
plot(MXanom(:,1),MXanom(:,9),'-r'); hold on;
ylim([-0.6 0.6])
xlim([MXanom(1,1) MXanom(end,1)])
title('GMX')
ylabel('fish')

%legend({'NAO','Sm','Md','Lg'},'location','southoutside','orientation','horizontal')
lg = legend(nexttile(3),{'NAO','Sm','Md','Lg'});
lg.Location = 'southoutside';
lg.Orientation = 'horizontal';
print('-dpng',[ppath 'ts_climate_seasonal_FOSI_NAtl_NAO_',mod,'size.png'])

%% put text of freq of var on fig
vS_EBS = var(BSanom(:,7));
sS_EBS = std(BSanom(:,7));
cvS_EBS = std(BSanom(:,7)) / mean(BSanom(:,7));

vM_EBS = var(BSanom(:,8));
sM_EBS = std(BSanom(:,8));
cvM_EBS = std(BSanom(:,8)) / mean(BSanom(:,8));

vL_EBS = var(BSanom(:,9));
sL_EBS = std(BSanom(:,9));
cvL_EBS = std(BSanom(:,9)) / mean(BSanom(:,9));

vS_CCE = var(CCanom(:,7));
sS_CCE = std(CCanom(:,7));
cvS_CCE = std(CCanom(:,7)) / mean(CCanom(:,7));

vM_CCE = var(CCanom(:,8));
sM_CCE = std(CCanom(:,8));
cvM_CCE = std(CCanom(:,8)) / mean(CCanom(:,8));

vL_CCE = var(CCanom(:,9));
sL_CCE = std(CCanom(:,9));
cvL_CCE = std(CCanom(:,9)) / mean(CCanom(:,9));

xi = BSanom(:,7);
inan = ~isnan(xi);
R = (xi(inan));
T = length(R);
t = (1:T)';
data = [t R];
[m,b] = TheilSen(data); 
tH = m*t + b;
dR = R - tH;
[freq1,p1,b1,int1] = powerspectra_TS(dR,1,0);
Sebs_b = b1;
Sebs_f = freq1(find(p1==max(p1)));
clear b1 freq1 p1

xi = BSanom(:,8);
inan = ~isnan(xi);
R = (xi(inan));
T = length(R);
t = (1:T)';
data = [t R];
[m,b] = TheilSen(data); 
tH = m*t + b;
dR = R - tH;
[freq1,p1,b1,int1] = powerspectra_TS(dR,1,0);
Mebs_b = b1;
Mebs_f = freq1(find(p1==max(p1)));
clear b1 freq1 p1

xi = BSanom(:,9);
inan = ~isnan(xi);
R = (xi(inan));
T = length(R);
t = (1:T)';
data = [t R];
[m,b] = TheilSen(data); 
tH = m*t + b;
dR = R - tH;
[freq1,p1,b1,int1] = powerspectra_TS(dR,1,0);
Lebs_b = b1;
Lebs_f = freq1(find(p1==max(p1)));
clear b1 freq1 p1

xi = CCanom(:,7);
inan = ~isnan(xi);
R = (xi(inan));
T = length(R);
t = (1:T)';
data = [t R];
[m,b] = TheilSen(data); 
tH = m*t + b;
dR = R - tH;
[freq1,p1,b1,int1] = powerspectra_TS(dR,1,0);
Scce_b = b1;
Scce_f = freq1(find(p1==max(p1)));
clear b1 freq1 p1

xi = CCanom(:,8);
inan = ~isnan(xi);
R = (xi(inan));
T = length(R);
t = (1:T)';
data = [t R];
[m,b] = TheilSen(data); 
tH = m*t + b;
dR = R - tH;
[freq1,p1,b1,int1] = powerspectra_TS(dR,1,0);
Mcce_b = b1;
Mcce_f = freq1(find(p1==max(p1)));
clear b1 freq1 p1

xi = CCanom(:,9);
inan = ~isnan(xi);
R = (xi(inan));
T = length(R);
t = (1:T)';
data = [t R];
[m,b] = TheilSen(data); 
tH = m*t + b;
dR = R - tH;
[freq1,p1,b1,int1] = powerspectra_TS(dR,1,0);
Lcce_b = b1;
Lcce_f = freq1(find(p1==max(p1)));
clear b1 freq1 p1

%%
%WHY ARE FISH AMOUNTS DIFF THAN ANN MEAN FIG?!
figure(5)
tiledlayout(3,2, 'TileSpacing', 'compact')
nexttile %S EBS
plot(BSanom(:,1),BSanom(:,7),'-k'); hold on;
xlim([BSanom(1,1) BSanom(end,1)])
ylabel('Sm')
title('EBS')
text(1950, 0.017,['b= ' num2str(Sebs_b)])

nexttile %S CCE
plot(CCanom(:,1),CCanom(:,7),'-k'); hold on;
xlim([CCanom(1,1) CCanom(end,1)])
title('CCE')
text(1950,1.8e-3,['b= ' num2str(Scce_b)])

nexttile %M EBS
plot(BSanom(:,1),BSanom(:,8),'-k'); hold on;
xlim([BSanom(1,1) BSanom(end,1)])
ylabel('Md')
text(1950, 0.4,['b= ' num2str(Mebs_b)])

nexttile %M CCE
plot(CCanom(:,1),CCanom(:,8),'-k'); hold on;
xlim([CCanom(1,1) CCanom(end,1)])
text(1950,0.04,['b= ' num2str(Mcce_b)])

nexttile %L EBS
plot(BSanom(:,1),BSanom(:,9),'-k'); hold on;
xlim([BSanom(1,1) BSanom(end,1)])
ylabel('Lg')
text(1950, 0.8,['b= ' num2str(Lebs_b)])

nexttile %L CCE
plot(CCanom(:,1),CCanom(:,9),'-k'); hold on;
xlim([CCanom(1,1) CCanom(end,1)])
text(1950,0.017,['b= ' num2str(Lcce_b)])

print('-dpng',[ppath 'ts_climate_seasonal_FOSI_EBS_CCE_',mod,'size_only.png'])

%% ---------- PDO & AMO ---------------------------------
%WHY ARE FISH AMOUNTS DIFF THAN ANN MEAN FIG?!
% PDO
f1 = figure('Units','inches','Position',[1 3 8 8]);

%CHK
%subplot(4,1,1)
tiledlayout(4,2, 'TileSpacing', 'compact')
nexttile
%climate = PDO
yyaxis left
plot(CKanom(:,1),manom(5,:),'color',0.5*[1 1 1]); hold on;
ylim([-1.5 1.5])
ylabel('PDO')
yyaxis right
%Sm
plot(CKanom(:,1),100*CKanom(:,7),'-k'); hold on;
%Md
plot(CKanom(:,1),CKanom(:,8),'-b'); hold on;
%Lg
plot(CKanom(:,1),CKanom(:,9),'-r'); hold on;
ylim([-1.5 1.5])
xlim([CKanom(1,1) CKanom(end,1)])
title('CHK')
ylabel('fish')

% AMO
%NE
nexttile
yyaxis left
plot(NEanom(:,1),manom(1,:),'color',0.5*[1 1 1]); hold on;
ylim([-0.5 0.5])
ylabel('AMO')
yyaxis right
%Sm
plot(NEanom(:,1),50*NEanom(:,7),'-k'); hold on;
%Md
plot(NEanom(:,1),NEanom(:,8),'-b'); hold on;
%Lg
plot(NEanom(:,1),NEanom(:,9),'-r'); hold on;
ylim([-0.2 0.2])
xlim([NEanom(1,1) NEanom(end,1)])
title('NE')
ylabel('fish')

nexttile %subplot(4,1,2)
%EBS
yyaxis left
plot(BSanom(:,1),manom(5,:),'color',0.5*[1 1 1]); hold on;
ylim([-1.5 1.5])
ylabel('PDO')
yyaxis right
%S
plot(BSanom(:,1),50*BSanom(:,7),'-k'); hold on;
%M
plot(BSanom(:,1),BSanom(:,8),'-b'); hold on;
%L
plot(BSanom(:,1),BSanom(:,9),'-r'); hold on;
ylim([-0.8 0.8])
xlim([BSanom(1,1) BSanom(end,1)])
title('EBS')
ylabel('fish')

nexttile %subplot(2,1,2)
% SE
yyaxis left
plot(SEanom(:,1),manom(1,:),'color',0.5*[1 1 1]); hold on;
ylim([-0.5 0.5])
ylabel('AMO')
yyaxis right
%S
plot(SEanom(:,1),25*SEanom(:,7),'-k'); hold on;
%M
plot(SEanom(:,1),SEanom(:,8),'-b'); hold on;
%L
plot(SEanom(:,1),SEanom(:,9),'-r'); hold on;
ylim([-0.08 0.08])
xlim([SEanom(1,1) SEanom(end,1)])
title('SE')
ylabel('fish')


% AK
nexttile %subplot(4,1,3)
%climate = PDO
yyaxis left
plot(AKanom(:,1),manom(5,:),'color',0.5*[1 1 1]); hold on;
ylim([-1.5 1.5])
ylabel('PDO')
yyaxis right
%S
plot(AKanom(:,1),10*AKanom(:,7),'-k'); hold on;
%M
plot(AKanom(:,1),AKanom(:,8),'-b'); hold on;
%L
plot(AKanom(:,1),AKanom(:,9),'-r'); hold on;
ylim([-0.25 0.25])
xlim([AKanom(1,1) AKanom(end,1)])
title('GAK')
ylabel('fish')

% GMX
nexttile %subplot(4,1,3)
%climate = AMO
yyaxis left
plot(MXanom(:,1),manom(1,:),'color',0.5*[1 1 1]); hold on;
ylim([-0.5 0.5])
ylabel('AMO')
yyaxis right
%S
plot(MXanom(:,1),50*MXanom(:,7),'-k'); hold on;
%M
plot(MXanom(:,1),MXanom(:,8),'-b'); hold on;
%L
plot(MXanom(:,1),MXanom(:,9),'-r'); hold on;
ylim([-0.6 0.6])
xlim([MXanom(1,1) MXanom(end,1)])
title('GMX')
ylabel('fish')

%CCE
nexttile %subplot(4,1,4)
yyaxis left
plot(CCanom(:,1),manom(5,:),'color',0.5*[1 1 1]); hold on;
ylim([-1.5 1.5])
ylabel('PDO')
yyaxis right
%S
plot(CCanom(:,1),10*CCanom(:,7),'-k'); hold on;
%M
plot(CCanom(:,1),CCanom(:,8),'-b'); hold on;
%L
plot(CCanom(:,1),CCanom(:,9),'-r'); hold on;
ylim([-0.05 0.05])
xlim([CCanom(1,1) CCanom(end,1)])
title('CCE')
ylabel('fish')


lg = legend(nexttile(6),{'climate','Sm','Md','Lg'});
lg.Location = 'southoutside';
lg.Orientation = 'horizontal';
print('-dpng',[ppath 'ts_climate_seasonal_FOSI_PDO_AMO_',mod,'size.png'])

