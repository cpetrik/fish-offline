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
ppath = [pp cfile '/FOSI/'];
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

pid = [54,1:2,3];
pname = {'CHK','EBS','GAK','CCE'};

aid = [7,6,5];
aname = {'NE','SE','GMX'};

% FEISTY outputs grouped
gname = {'S','M','L'};

%%
yr = 1948:2015; %=yanom;

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
tiledlayout(4,1, 'TileSpacing', 'compact')

for i=1:length(pid)
    
    %scale fish -1 to 1 ? - want to keep 0 at 0
    ms = max(abs(as(pid(i),:)));
    mm = max(abs(am(pid(i),:)));
    ml = max(abs(al(pid(i),:)));
    
    sm = as(pid(i),:)./ms;
    md = am(pid(i),:)./mm;
    lg = al(pid(i),:)./ml;
    
    nexttile
    %climate = PDO
    yyaxis left
    plot(yr,manom(5,:),'color',0.5*[1 1 1]); hold on;
    ylim([-2.1 2.1])
    ylabel('PDO')
    yyaxis right
    %Sm
    plot(yr,sm,'-k'); hold on;
    %Md
    plot(yr,md,'-b'); hold on;
    %Lg
    plot(yr,lg,'-r'); hold on;
    ylim([-1 1])
    xlim([yr(1) yr(end)])
    title(pname{i})
    ylabel('fish')
    
end
lg = legend(nexttile(4),{'PDO','Sm','Md','Lg'});
lg.Location = 'southoutside';
lg.Orientation = 'horizontal';
print('-dpng',[ppath 'ts_climate_seasonal_FOSI_NPac_PDO_',mod,'size_v2.png'])

%% ENSO
figure(2)
tiledlayout(4,1, 'TileSpacing', 'compact')
for i=1:length(pid)
    
    %scale fish -1 to 1 ? - want to keep 0 at 0
    ms = max(abs(as(pid(i),:)));
    mm = max(abs(am(pid(i),:)));
    ml = max(abs(al(pid(i),:)));
    
    sm = as(pid(i),:)./ms;
    md = am(pid(i),:)./mm;
    lg = al(pid(i),:)./ml;
    
    nexttile
    %climate = ENSO
    yyaxis left
    plot(yr,manom(4,:),'color',0.5*[1 1 1]); hold on;
    ylim([-1.5 1.5])
    ylabel('ENSO')
    yyaxis right
    %Sm
    plot(yr,sm,'-k'); hold on;
    %Md
    plot(yr,md,'-b'); hold on;
    %Lg
    plot(yr,lg,'-r'); hold on;
    ylim([-1 1])
    xlim([yr(1) yr(end)])
    title(pname{i})
    ylabel('fish')
    
end
lg = legend(nexttile(4),{'ENSO','Sm','Md','Lg'});
lg.Location = 'southoutside';
lg.Orientation = 'horizontal';
print('-dpng',[ppath 'ts_climate_seasonal_FOSI_NPac_ENSO_',mod,'size_v2.png'])

%% AMO
figure(3)
tiledlayout(3,1, 'TileSpacing', 'compact')
for i=1:length(aid)
    
    %scale fish -1 to 1 ? - want to keep 0 at 0
    ms = max(abs(as(aid(i),:)));
    mm = max(abs(am(aid(i),:)));
    ml = max(abs(al(aid(i),:)));
    
    sm = as(aid(i),:)./ms;
    md = am(aid(i),:)./mm;
    lg = al(aid(i),:)./ml;
    
    nexttile
    %climate = AMO
    yyaxis left
    plot(yr,manom(1,:),'color',0.5*[1 1 1]); hold on;
    ylim([-0.45 0.45])
    ylabel('AMO')
    yyaxis right
    %Sm
    plot(yr,sm,'-k'); hold on;
    %Md
    plot(yr,md,'-b'); hold on;
    %Lg
    plot(yr,lg,'-r'); hold on;
    ylim([-1 1])
    xlim([yr(1) yr(end)])
    title(aname{i})
    ylabel('fish')
    
end
lg = legend(nexttile(3),{'AMO','Sm','Md','Lg'});
lg.Location = 'southoutside';
lg.Orientation = 'horizontal';
print('-dpng',[ppath 'ts_climate_seasonal_FOSI_NAtl_AMO_',mod,'size_v2.png'])

%% NAO
figure(4)
tiledlayout(3,1, 'TileSpacing', 'compact')
for i=1:length(aid)
    
    %scale fish -1 to 1 ? - want to keep 0 at 0
    ms = max(abs(as(aid(i),:)));
    mm = max(abs(am(aid(i),:)));
    ml = max(abs(al(aid(i),:)));
    
    sm = as(aid(i),:)./ms;
    md = am(aid(i),:)./mm;
    lg = al(aid(i),:)./ml;
    
    nexttile
    %climate = AMO
    yyaxis left
    plot(yr,manom(3,:),'color',0.5*[1 1 1]); hold on;
    ylim([-3.5 3.5])
    ylabel('NAO')
    yyaxis right
    %Sm
    plot(yr,sm,'-k'); hold on;
    %Md
    plot(yr,md,'-b'); hold on;
    %Lg
    plot(yr,lg,'-r'); hold on;
    ylim([-1 1])
    xlim([yr(1) yr(end)])
    title(aname{i})
    ylabel('fish')
    
end
lg = legend(nexttile(3),{'NAO','Sm','Md','Lg'});
lg.Location = 'southoutside';
lg.Orientation = 'horizontal';
print('-dpng',[ppath 'ts_climate_seasonal_FOSI_NAtl_NAO_',mod,'size_v2.png'])

%% put text of freq of var on fig
% EBS
xi = as(1,:)';
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

xi = am(1,:)';
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

xi = al(1,:)';
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

%CCE
xi = as(3,:)';
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

xi = am(3,:)';
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

xi = al(3,:)';
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

%NE
xi = as(7,:)';
inan = ~isnan(xi);
R = (xi(inan));
T = length(R);
t = (1:T)';
data = [t R];
[m,b] = TheilSen(data);
tH = m*t + b;
dR = R - tH;
[freq1,p1,b1,int1] = powerspectra_TS(dR,1,0);
Sne_b = b1;
Sne_f = freq1(find(p1==max(p1)));
clear b1 freq1 p1

xi = am(7,:)';
inan = ~isnan(xi);
R = (xi(inan));
T = length(R);
t = (1:T)';
data = [t R];
[m,b] = TheilSen(data);
tH = m*t + b;
dR = R - tH;
[freq1,p1,b1,int1] = powerspectra_TS(dR,1,0);
Mne_b = b1;
Mne_f = freq1(find(p1==max(p1)));
clear b1 freq1 p1

xi = al(7,:)';
inan = ~isnan(xi);
R = (xi(inan));
T = length(R);
t = (1:T)';
data = [t R];
[m,b] = TheilSen(data);
tH = m*t + b;
dR = R - tH;
[freq1,p1,b1,int1] = powerspectra_TS(dR,1,0);
Lne_b = b1;
Lne_f = freq1(find(p1==max(p1)));
clear b1 freq1 p1

%%
%WHY ARE FISH AMOUNTS DIFF THAN ANN MEAN FIG?!
figure(5)
tiledlayout(3,3, 'TileSpacing', 'compact')
nexttile %S EBS
plot(yr,as(1,:),'-k'); hold on;
xlim([yr(1) yr(end)])
ylabel('Sm')
title('EBS')
text(1950, 0.017,['b= ' sprintf('%2.2f',Sebs_b)])

nexttile %S CCE
plot(yr,as(3,:),'-k'); hold on;
xlim([yr(1) yr(end)])
title('CCE')
text(1950,0.008,['b= ' sprintf('%2.2f',Scce_b)])

nexttile %S NE
plot(yr,as(7,:),'-k'); hold on;
xlim([yr(1) yr(end)])
title('NE')
text(1950,0.017,['b= ' sprintf('%2.2f',Sne_b)])

nexttile %M EBS
plot(yr,am(1,:),'-k'); hold on;
xlim([yr(1) yr(end)])
ylabel('Md')
text(1950, 0.4,['b= ' sprintf('%2.2f',Mebs_b)])

nexttile %M CCE
plot(yr,am(3,:),'-k'); hold on;
xlim([yr(1) yr(end)])
text(1950,0.4,['b= ' sprintf('%2.2f',Mcce_b)])

nexttile %M NE
plot(yr,am(7,:),'-k'); hold on;
xlim([yr(1) yr(end)])
text(1950,0.4,['b= ' sprintf('%2.2f',Mne_b)])

nexttile %L EBS
plot(yr,al(1,:),'-k'); hold on;
xlim([yr(1) yr(end)])
ylabel('Lg')
text(1950, 0.8,['b= ' sprintf('%2.2f',Lebs_b)])

nexttile %L CCE
plot(yr,al(3,:),'-k'); hold on;
xlim([yr(1) yr(end)])
text(1950,0.16,['b= ' sprintf('%2.2f',Lcce_b)])

nexttile %L NE
plot(yr,al(7,:),'-k'); hold on;
xlim([yr(1) yr(end)])
text(1950,0.5,['b= ' sprintf('%2.2f',Lne_b)])

print('-dpng',[ppath 'ts_climate_seasonal_FOSI_EBS_CCE_',mod,'size_only_v2.png'])

%% ---------- PDO & AMO ---------------------------------
f1 = figure('Units','inches','Position',[1 3 8 8]);
tiledlayout(4,2, 'TileSpacing', 'compact')
% PDO
for i=1:length(pid)
    %scale fish -1 to 1 - want to keep 0 at 0
    ms = max(abs(as(pid(i),:)));
    mm = max(abs(am(pid(i),:)));
    ml = max(abs(al(pid(i),:)));
    
    sm = as(pid(i),:)./ms;
    md = am(pid(i),:)./mm;
    lg = al(pid(i),:)./ml;
    
    nexttile((2*i)-1)
    %climate = PDO
    yyaxis left
    plot(yr,manom(5,:),'color',0.5*[1 1 1]); hold on;
    ylim([-2.1 2.1])
    ylabel('PDO')
    yyaxis right
    %Sm
    plot(yr,sm,'-k'); hold on;
    %Md
    plot(yr,md,'-b'); hold on;
    %Lg
    plot(yr,lg,'-r'); hold on;
    ylim([-1 1])
    xlim([yr(1) yr(end)])
    title(pname{i})
    ylabel('fish')
    
end
% AMO
for i=1:length(aid)
    %scale fish -1 to 1  - want to keep 0 at 0
    ms = max(abs(as(aid(i),:)));
    mm = max(abs(am(aid(i),:)));
    ml = max(abs(al(aid(i),:)));
    
    sm = as(aid(i),:)./ms;
    md = am(aid(i),:)./mm;
    lg = al(aid(i),:)./ml;
    
    nexttile((2*i))
    %climate = AMO
    yyaxis left
    plot(yr,manom(1,:),'color',0.5*[1 1 1]); hold on;
    ylim([-0.45 0.45])
    ylabel('AMO')
    yyaxis right
    %Sm
    plot(yr,sm,'-k'); hold on;
    %Md
    plot(yr,md,'-b'); hold on;
    %Lg
    plot(yr,lg,'-r'); hold on;
    ylim([-1 1])
    xlim([yr(1) yr(end)])
    title(aname{i})
    ylabel('fish')
    
end
lg = legend(nexttile(6),{'climate','Sm','Md','Lg'});
lg.Location = 'southoutside';
lg.Orientation = 'horizontal';
print('-dpng',[ppath 'ts_climate_seasonal_FOSI_PDO_AMO_',mod,'size_v2.png'])

