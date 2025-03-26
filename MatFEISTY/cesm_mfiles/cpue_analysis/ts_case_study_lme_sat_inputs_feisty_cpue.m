% Plot ts of LME CPUE against most significant driver
% that is not satellite SST or chl
% CESM FOSI
% 1997-2015

clear
close all

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',cfile,'/corrs_cpue/'];
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

mod = 'v15_All_fish03';
mod2 = 'v15_obsfish';

%% Sat
load([fpath 'lme_satellite_sst_chl_ann_mean_anoms_1997_2010_2015.mat'],...
    'achl15','asst15','yyr');

%% FOSI input forcing
% lme means, trend removed, anomaly calc
load([cpath 'CESM_FOSI_v15_lme_interann_mean_forcings_anom_1997_2010_2015.mat'],...
    'adet15','atb15','atp15','azlos15');

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% FEISTY outputs
% Constant effort

% Anoms with linear trend removed
%Biomass
load([fpath 'FEISTY_FOSI_',mod,'_lme_biom_ann_mean_anoms_1997_2010_2015.mat'],...
    'aba15','abd15','abf15','abp15');

% Nu
load([fpath 'FEISTY_FOSI_',mod,'_lme_nu_ann_mean_anoms_1997_2010_2015.mat'],...
    'ana15','and15','anf15','anp15');

%% Obs effort
% 
% % Anoms with linear trend removed
% % Biomass
% load([fpath 'FEISTY_FOSI_',mod2,'_lme_biom_ann_mean_anoms_1997_2015.mat'],...
%     'aba15','abd15','abf15','abp15','eyr');
% 
% % Nu
% load([fpath 'FEISTY_FOSI_',mod2,'_lme_nu_ann_mean_anoms_1997_2015.mat'],...
%     'ana15','and15','anf15','anp15');
% 
% % Yield/Catch
% load([fpath 'FEISTY_FOSI_',mod2,'_lme_yield_ann_mean_anoms_1997_2015.mat'],...
%     'aya15','ayf15','ayp15','ayd15');
% 
% Oaba = aa;
% Oabd = ad;
% Oabf = af;
% Oabp = ap;
% 
% Oana = aa;
% Oand = ad;
% Oanf = af;
% Oanp = ap;

%% Fishing data
% Anoms with linear trend removed
load([ypath 'FishMIP_Phase3a_LME_CPUE_1997-2015_ann_mean_anoms.mat'])

%% colors
load('paul_tol_cmaps.mat')

%colorblind friendly - subselect & re-order drainbow
ctex = {'TP','TB','Det','ZmLoss','SST','Chl','Biom','Prod'};
% orange, dk blue, grey, lt blue, dk purp, lt purp, red, green
mcol(1,:) = drainbow(15,:)/255; %grey
mcol(2,:) = [0 0 0]; %black
mcol(3,:) = drainbow(12,:)/255; % orange
mcol(4,:) = drainbow(4,:)/255; %dk blue
mcol(5,:) = drainbow(6,:)/255; %lt blue
mcol(6,:) = drainbow(14,:)/255; %red
mcol(7,:) = drainbow(7,:)/255; %green
mcol(8,:) = drainbow(3,:)/255; %dk purp
mcol(9,:) = drainbow(1,:)/255; %lt purp

ctp = drainbow(12,:)/255; % orange
ctb = drainbow(4,:)/255; %dk blue
cdet = drainbow(15,:)/255; %grey
czm = drainbow(6,:)/255; %lt blue
csst = drainbow(14,:)/255; %red
cchl = drainbow(7,:)/255; %green
cbio = drainbow(3,:)/255; %dk purp
cnu = drainbow(1,:)/255; %lt purp

colororder(mcol)
close all

%% figure info
axesPosition = [130 40 400 200];  %# Axes position, in pixels
yWidth = 30;                      %# y axes spacing, in pixels
xLimit = [yyr(1) yyr(end)];       %# Range of x values
xOffset = -yWidth*diff(xLimit)/axesPosition(3);

% PvD with ZprodDet, with temp
% figure('Units','pixels','Position',[200 200 330 260]);
% h1 = axes('Units','pixels','Position',axesPosition,...
%           'Color','w','XColor','k','YColor',[0.5 0.5 0.5],...
%           'XLim',xLimit,'YLim',[-0.06 0.02],'NextPlot','add');
% h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
%           'Color','none','XColor','k','YColor',[0.5 0 1.0],...
%           'XLim',xLimit+[xOffset 0],'YLim',[-0.1 0.4],...
%           'XTick',[],'XTickLabel',[],'NextPlot','add');
% h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
%           'Color','none','XColor','k','YColor','k',...
%           'XLim',xLimit+[2*xOffset 0],'YLim',[14.5 17.5],...
%           'XTick',[],'XTickLabel',[],'NextPlot','add');
% xlabel(h1,'Year');
% 
% line(y(19:47),tPD(19:47),'color',cmap_ppt(1,:),'Linewidth',2,'Parent',h1); hold on;
% line(y(19:47),tZD(19:47),'color',cm{6},'Linewidth',2,'Parent',h2); hold on;
% line(y(19:47),mtemp(19:47),'color','k','Linewidth',2,'Parent',h3); hold on;


%% LME 2
% LME slag	sat	dlag	driver	flag	fdriver
% 2	    1	SST	1	    TP	    1	    TP

figure(1)
%plot without moving lags
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',csst,...
          'XLim',xLimit,'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',ctp,...
          'XLim',xLimit+[xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');

line(yyr,asst15(2,:),'color',csst,'Parent',h1); hold on;
line(yyr,atp15(2,:),'color',ctp,'Parent',h2); 
%legend('SST','TP')
%legend('location','northwest')
xlabel(h1,'Year');

line(yyr,aa_cpue97(2,:),'color','k','Parent',h3); 
ylabel(h3,'Anomalies')
%legend('CPUE')
%legend('location','northwest')
title('LME 2, No lags')
print('-dpng',[ppath 'ts_LME2_CPUE_sig_drivers_nolag.png'])

%% plot with lags
dyr1 = yyr(1:(end-1)) + 1;
cyr1 = yyr(2:end);

figure(2)
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',csst,...
          'XLim',xLimit,'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',ctp,...
          'XLim',xLimit+[xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');

line(dyr1,asst15(2,(1:(end-1))),'color',csst,'Parent',h1); hold on;
line(dyr1,atp15(2,(1:(end-1))),'color',ctp,'Parent',h2); 
xlabel(h1,'Year');

line(cyr1,aa_cpue97(2,(2:end)),'color','k','Parent',h3); 
ylabel(h3,'Anomalies')

title('LME 2, Drivers -1 y lag')
print('-dpng',[ppath 'ts_LME2_CPUE_sig_drivers_lags.png'])

%% LME 5
% LME slag	sat	dlag	driver	flag	fdriver
% 5	    0	SST	1	    ZmL	    1	    ZmL

figure(3)
%plot without moving lags
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',csst,...
          'XLim',xLimit,'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',czm,...
          'XLim',xLimit+[xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');

line(yyr,asst15(5,:),'color',csst,'Parent',h1); hold on;
line(yyr,azlos15(5,:),'color',czm,'Parent',h2); 
xlabel(h1,'Year');

line(yyr,aa_cpue97(5,:),'color','k','Parent',h3); 
ylabel(h3,'Anomalies')
title('LME 5, No lags')
print('-dpng',[ppath 'ts_LME5_CPUE_sig_drivers_nolag.png'])

%%
figure(4)
%plot with lags
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',csst,...
          'XLim',xLimit,'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',czm,...
          'XLim',xLimit+[xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');

line(yyr,asst15(5,:),'color',csst,'Parent',h1); hold on;
line(dyr1,azlos15(5,(1:(end-1))),'color',czm,'Parent',h2); 
xlabel(h1,'Year');

line(yyr,aa_cpue97(5,:),'color','k','Parent',h3); 
ylabel(h3,'Anomalies')
title('LME 5, SST0, ZmLoss-1')
print('-dpng',[ppath 'ts_LME5_CPUE_sig_drivers_lags.png'])


%% LME 6
% LME slag	sat	dlag	driver	flag	fdriver
% 6	    3   Chl	2	    TB	    2	    TB
dyr2 = yyr(1:(end-2)) + 2;
cyr2 = yyr(3:end);
dyr3 = yyr(1:(end-3)) + 3;
cyr3 = yyr(4:end);

figure(5)
%plot without moving lags
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',cchl,...
          'XLim',xLimit,'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',ctb,...
          'XLim',xLimit+[xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');

line(yyr,achl15(6,:),'color',cchl,'Parent',h1); hold on;
line(yyr,atb15(6,:),'color',ctb,'Parent',h2); 
xlabel(h1,'Year');

line(yyr,aa_cpue97(6,:),'color','k','Parent',h3); 
ylabel(h3,'Anomalies')
title('LME 6, No lags')
print('-dpng',[ppath 'ts_LME6_CPUE_sig_drivers_nolag.png'])

%% GOOD ONE!!
figure(6)
%plot with lags
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',cchl,...
          'XLim',xLimit,'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',ctb,...
          'XLim',xLimit+[xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');

line(dyr3,achl15(6,(1:(end-3))),'color',cchl,'Parent',h1); hold on;
line(dyr2,atb15(6,(1:(end-2))),'color',ctb,'Parent',h2); 
xlabel(h1,'Year');

line(cyr2,aa_cpue97(6,(3:end)),'color','k','Parent',h3); 
ylabel(h3,'Anomalies')
title('LME 5, Chl-3, TB-2')
print('-dpng',[ppath 'ts_LME6_CPUE_sig_drivers_lags.png'])


%% LME 12
% LME slag	sat	dlag	driver	flag	fdriver
% 12	2	SST	1	    Det	    0	    Biom

figure(7)
%plot without moving lags
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',csst,...
          'XLim',xLimit,'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',cdet,...
          'XLim',xLimit+[xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor',cbio,...
          'XLim',xLimit+[2*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h4 = axes('Units','pixels','Position',axesPosition+yWidth.*[-3 0 3 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[3*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');

line(yyr,asst15(12,:),'color',csst,'Parent',h1); hold on;
line(yyr,adet15(12,:),'color',cdet,'Parent',h2); hold on;
line(yyr,aba15(12,:),'color',cbio,'Parent',h3);
xlabel(h1,'Year');

line(yyr,aa_cpue97(12,:),'color','k','Parent',h4); 
ylabel(h4,'Anomalies')
title('LME 12, No lags')
print('-dpng',[ppath 'ts_LME12_CPUE_sig_drivers_nolag.png'])

%% MAYBE
% LME slag	sat	dlag	driver	flag	fdriver
% 12	2	SST	1	    Det	    0	    Biom

figure(8)
%plot with lags
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',csst,...
          'XLim',xLimit,'YTick',-0.4:0.1:0.3,'YTickLabel',[0.3,0.2,0.1,0,-0.1,-0.2,-0.3,-0.4],...
          'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',cdet,...
          'XLim',xLimit+[xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor',cbio,...
          'XLim',xLimit+[2*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h4 = axes('Units','pixels','Position',axesPosition+yWidth.*[-3 0 3 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[3*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');

line(dyr2,-1*asst15(12,(1:(end-2))),'color',csst,'Parent',h1); hold on;
line(dyr1,adet15(12,(1:(end-1))),'color',cdet,'Parent',h2); hold on;
line(yyr,aba15(12,:),'color',cbio,'Parent',h3);
xlabel(h1,'Year');

line(yyr,aa_cpue97(12,:),'color','k','Parent',h4); 
ylabel(h4,'Anomalies')
title('LME 5, SST-2, Det-1, Biom0')
print('-dpng',[ppath 'ts_LME12_CPUE_sig_drivers_lags.png'])


%% LME 15
% LME slag	sat	dlag	driver	flag	fdriver
% 15	0	SST	1	    Det	    0	    Biom

figure(9)
%plot without moving lags
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',csst,...
          'XLim',xLimit,'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',cdet,...
          'XLim',xLimit+[xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor',cbio,...
          'XLim',xLimit+[2*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h4 = axes('Units','pixels','Position',axesPosition+yWidth.*[-3 0 3 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[3*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');

line(yyr,asst15(15,:),'color',csst,'Parent',h1); hold on;
line(yyr,adet15(15,:),'color',cdet,'Parent',h2); hold on;
line(yyr,aba15(15,:),'color',cbio,'Parent',h3);
xlabel(h1,'Year');

line(yyr,aa_cpue97(12,:),'color','k','Parent',h4); 
ylabel(h4,'Anomalies')
title('LME 15, No lags')
print('-dpng',[ppath 'ts_LME15_CPUE_sig_drivers_nolag.png'])

%% 
% LME slag	sat	dlag	driver	flag	fdriver
% 15	0	SST	1	    Det	    0	    Biom

figure(10)
%plot with lags
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',csst,...
          'XLim',xLimit,'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',cdet,...
          'XLim',xLimit+[xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor',cbio,...
          'XLim',xLimit+[2*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h4 = axes('Units','pixels','Position',axesPosition+yWidth.*[-3 0 3 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[3*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');

line(yyr,asst15(15,:),'color',csst,'Parent',h1); hold on;
line(dyr1,adet15(15,(1:(end-1))),'color',cdet,'Parent',h2); hold on;
line(yyr,aba15(15,:),'color',cbio,'Parent',h3);
xlabel(h1,'Year');

line(yyr,aa_cpue97(15,:),'color','k','Parent',h4); 
ylabel(h4,'Anomalies')
title('LME 15, SST0, Det-1, Biom0')
print('-dpng',[ppath 'ts_LME15_CPUE_sig_drivers_lags.png'])


%% LME 18
% LME slag	sat	dlag	driver	flag	fdriver
% 18	0	SST	0	    ZmL	    0	    ZmL

figure(11)
%plot without moving lags
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',csst,...
          'XLim',xLimit,'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',czm,...
          'XLim',xLimit+[xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');

line(yyr,asst15(18,:),'color',csst,'Parent',h1); hold on;
line(yyr,azlos15(18,:),'color',czm,'Parent',h2); 
xlabel(h1,'Year');

line(yyr,aa_cpue97(18,:),'color','k','Parent',h3); 
ylabel(h3,'Anomalies')
title('LME 18, No lags')
print('-dpng',[ppath 'ts_LME18_CPUE_sig_drivers_nolag.png'])





