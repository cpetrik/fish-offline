% Satellite SST data
% 1982-2015

clear 
close all

%%
tpath = '/Volumes/petrik-lab/Feisty/Obs_data/OISST/';

rpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/ImptDocs/Observations/StockAssessments/RAM/RAMLDB v4.66/Excel/';

fpath='/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';

%%
load([tpath 'lme_satellite_sst_interann_means_1982_2020.mat']);

load([fpath 'CESM_FOSI_v15_lme_interann_mean_drivers_1948_2015.mat']);

load([rpath 'pollock_ssb_tb_rec_RAM.mat'])

%% same time period SSB
[Ysst,sid] = intersect(Year,tyr);
[Ytb,fid] = intersect(Year,fyr);
[Ytb,bid] = intersect(fyr,Year);

sst = lme_sst_ts(1,:);
ssb_sst = SSB(sid);

tb = lme_tb_mean(1,bid);
ssb_tb = SSB(fid);

tb_sst = TB(sid);
tb_tb = TB(fid);

rec_sst = Recruit(sid);
rec_tb = Recruit(fid);

%% plots
%colors
% blue 15, 111, 198
%pink 255, 47, 146
blue = [15/255, 111/255, 198/255];
pink = [255/255, 47/255, 146/255];

figure(1)
colororder([blue;pink])
yyaxis left
plot(Ysst,ssb_sst,'color',blue,'LineWidth',2)

figure(2)
colororder([blue;pink])
yyaxis left
plot(Ysst,ssb_sst,'color',blue,'LineWidth',2); hold on
yyaxis right
plot(Ysst,sst,'color',pink,'LineWidth',2);

figure(3)
colororder([blue;pink])
yyaxis left
plot(Ytb,ssb_tb,'color',blue,'LineWidth',2)

figure(4)
colororder([blue;pink])
yyaxis left
plot(Ytb,ssb_tb,'color',blue,'LineWidth',2); hold on
yyaxis right
plot(Ytb,tb,'color',pink,'LineWidth',2);

%%
figure(5)
plot(sst,ssb_sst,'o','MarkerFaceColor',blue,'MarkerSize',10)
axis square

figure(6)
plot(tb,ssb_tb,'o','MarkerFaceColor',blue,'MarkerSize',10)
axis square


%% plots TB
figure(11)
colororder([blue;pink])
yyaxis left
plot(Ysst,tb_sst,'color',blue,'LineWidth',2)

figure(12)
colororder([blue;pink])
yyaxis left
plot(Ysst,tb_sst,'color',blue,'LineWidth',2); hold on
yyaxis right
plot(Ysst,sst,'color',pink,'LineWidth',2);

figure(13)
colororder([blue;pink])
yyaxis left
plot(Ytb,tb_tb,'color',blue,'LineWidth',2)

figure(14)
colororder([blue;pink])
yyaxis left
plot(Ytb,tb_tb,'color',blue,'LineWidth',2); hold on
yyaxis right
plot(Ytb,tb,'color',pink,'LineWidth',2);

figure(15)
plot(sst,tb_sst,'o','MarkerFaceColor',blue,'MarkerSize',10)
axis square

figure(16)
plot(tb,tb_tb,'o','MarkerFaceColor',blue,'MarkerSize',10)
axis square

%% plots Recruit
figure(21)
colororder([blue;pink])
yyaxis left
plot(Ysst,rec_sst,'color',blue,'LineWidth',2)

figure(22)
colororder([blue;pink])
yyaxis left
plot(Ysst,rec_sst,'color',blue,'LineWidth',2); hold on
yyaxis right
plot(Ysst,sst,'color',pink,'LineWidth',2);

figure(23)
colororder([blue;pink])
yyaxis left
plot(Ytb,rec_tb,'color',blue,'LineWidth',2)

figure(24)
colororder([blue;pink])
yyaxis left
plot(Ytb,rec_tb,'color',blue,'LineWidth',2); hold on
yyaxis right
plot(Ytb,tb,'color',pink,'LineWidth',2);

figure(25)
plot(sst,rec_sst,'o','MarkerFaceColor',blue,'MarkerSize',10)
axis square

figure(26)
plot(tb,rec_tb,'o','MarkerFaceColor',blue,'MarkerSize',10)
axis square






