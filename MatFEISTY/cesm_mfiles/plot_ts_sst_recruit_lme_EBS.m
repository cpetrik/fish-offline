% Satellite SST data
% 1982-2015

clear 
close all

%%
tpath = '/Volumes/petrik-lab/Feisty/Obs_data/OISST/';

rpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/ImptDocs/Observations/StockAssessments/RAM/RAMLDB v4.66/Excel/';

fpath='/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';

ppath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/SIO/Research_summary/Birch_talk/';

%%
load([tpath 'lme_satellite_sst_interann_means_1982_2020.mat']);

load([rpath 'pollock_ssb_tb_rec_RAM.mat'])

%% same time period SSB
[Ysst,sid] = intersect(Year,tyr);

sst = lme_sst_ts(1,:);
rec_sst = Recruit(sid);

%% plots
%colors
% blue 15, 111, 198
%pink 255, 47, 146
blue = [15/255, 111/255, 198/255];
pink = [255/255, 47/255, 146/255];


%% plots Recruit
figure(1)
colororder([blue;pink])
yyaxis left
plot(Ysst,rec_sst,'color',blue,'LineWidth',2)
xlim([1980 2021])
xlabel('Year')
ylabel('Recruits')

figure(2)
colororder([blue;pink])
yyaxis left
plot(Ysst,rec_sst,'color',blue,'LineWidth',2); hold on
xlim([1980 2021])
xlabel('Year')
ylabel('Recruits')
yyaxis right
plot(Ysst,sst,'color',pink,'LineWidth',2);
xlim([1980 2021])
ylabel('SST')

%%
fig3 = figure(3);
plot(sst,rec_sst,'o','MarkerFaceColor',blue,'MarkerSize',10)
axis square
xlabel('SST')
ylabel('Recruits')
fontsize(fig3, 16, 'points')
print('-dpng',[ppath 'EBS_pollock_recruit_vs_sst.png'])

%%
mdl = fitlm(sst,rec_sst)

%% Plots as anomalies
msst = mean(sst,'omitnan');
mrec = mean(rec_sst,'omitnan');

arec_sst = rec_sst - mrec;
asst = sst - msst;

fig4=figure(4);
colororder([blue;pink])
yyaxis left
plot(Ysst,arec_sst,'color',blue,'LineWidth',2)
xlim([1980 2021])
ylim([-3e10 3e10])
xlabel('Year')
ylabel('Recruits')
yyaxis right
ylim([-1.5 1.5])
fontsize(fig4, 16, 'points')
print('-dpng',[ppath 'EBS_pollock_recruit_ts.png'])

fig5=figure(5);
colororder([blue;pink])
yyaxis left
plot(Ysst,arec_sst,'color',blue,'LineWidth',2); hold on
xlim([1980 2021])
ylim([-3e10 3e10])
xlabel('Year')
ylabel('Recruits')
yyaxis right
plot(Ysst,asst,'color',pink,'LineWidth',2);
xlim([1980 2021])
ylim([-1.5 1.5])
ylabel('SST')
fontsize(fig5, 16, 'points')
print('-dpng',[ppath 'EBS_pollock_recruit_sst_ts.png'])


% figure(6)
% plot(asst,arec_sst,'o','MarkerFaceColor',blue,'MarkerSize',10)
% axis square
% xlabel('SST')
% ylabel('Recruits')
