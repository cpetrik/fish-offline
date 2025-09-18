% Plot ts of LME CPUE against most significant driver
% that is not satellite SST or chl
% CESM FOSI
% 1997-2015
% See if significance related to t.s. length, and SST does better with
% longer record

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

%% time ranges for lags
% move drivers forward in time
dyr1 = yyr(1:(end-1)) + 1;
dyr2 = yyr(1:(end-2)) + 2;
dyr3 = yyr(1:(end-3)) + 3;

% shorten CPUE to match
cyr1 = yyr(2:end);
cyr2 = yyr(3:end);
cyr3 = yyr(4:end);


%% LME 6
% GOOD ONE!!
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

line(dyr3,achl15(6,(1:(end-3))),'color',cchl,'LineWidth',2,'Parent',h1); hold on;
line(dyr2,atb15(6,(1:(end-2))),'color',ctb,'LineWidth',2,'Parent',h2); 
xlabel(h1,'Year');

line(cyr2,aa_cpue97(6,(3:end)),'color','k','LineWidth',2,'Parent',h3); 
ylabel(h3,'Anomalies')
title('LME 6, Chl-3, TB-2')
print('-dpng',[ppath 'ts_LME6_CPUE_sig_drivers_chl15_bestlags.png'])

%% LME 12
% MAYBE
% LME slag	sat	dlag	driver	flag	fdriver
% 12	2	SST	1	    Det	    0	    Biom

figure(12)
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

line(dyr2,-1*asst15(12,(1:(end-2))),'color',csst,'LineWidth',2,'Parent',h1); hold on;
line(dyr1,adet15(12,(1:(end-1))),'color',cdet,'LineWidth',2,'Parent',h2); hold on;
line(yyr,aba15(12,:),'color',cbio,'LineWidth',2,'Parent',h3);
xlabel(h1,'Year');

line(yyr,aa_cpue97(12,:),'color','k','LineWidth',2,'Parent',h4); 
ylabel(h4,'Anomalies')
title('LME 12, SST-2, Det-1, Biom0')
print('-dpng',[ppath 'ts_LME12_CPUE_sig_drivers_chl15_bestlags.png'])

%% LME 26
% MAYBE
% LME slag	sat	dlag	driver	flag	fdriver
% 26	3	Chl	1	    ZmL	    1	    ZmL

figure(26)
%plot with lags
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',cchl,...
          'XLim',xLimit,'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',czm,...
          'XLim',xLimit+[xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');

line(dyr3,achl15(26,(1:(end-3))),'color',cchl,'LineWidth',2,'Parent',h1); hold on;
line(dyr1,azlos15(26,(1:(end-1))),'color',czm,'LineWidth',2,'Parent',h2); hold on;
xlabel(h1,'Year');

line(cyr1,aa_cpue97(26,(2:end)),'color','k','LineWidth',2,'Parent',h3); 
ylabel(h3,'Anomalies')
title('LME 26, Chl-3, ZmL-1')
print('-dpng',[ppath 'ts_LME26_CPUE_sig_drivers_chl15_bestlags.png'])

%% LME 31
% Pretty good
% LME slag	sat	dlag	driver	flag	fdriver
% 31	3	SST	1	    TP	    1	    TP

figure(31)
%plot with lags
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

line(dyr3,asst15(31,(1:(end-3))),'color',csst,'LineWidth',2,'Parent',h1); hold on;
line(dyr1,atp15(31,(1:(end-1))),'color',ctp,'LineWidth',2,'Parent',h2); hold on;
xlabel(h1,'Year');

line(cyr1,aa_cpue97(31,(2:end)),'color','k','LineWidth',2,'Parent',h3); 
ylabel(h3,'Anomalies')
title('LME 31, SST-3, TP-1')
print('-dpng',[ppath 'ts_LME31_CPUE_sig_drivers_chl15_bestlags.png'])

%% LME 48 
% Maybe
% LME slag	sat	dlag	driver	flag	fdriver
% 48	1	SST	3	    Det	    2	    Biom

figure(48)
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

line(dyr1,asst15(48,(1:(end-1))),'color',csst,'LineWidth',2,'Parent',h1); hold on;
line(dyr3,adet15(48,(1:(end-3))),'color',cdet,'LineWidth',2,'Parent',h2); hold on;
line(dyr2,aba15(48,(1:(end-2))),'color',cbio,'LineWidth',2,'Parent',h3); hold on;
xlabel(h1,'Year');

line(cyr1,aa_cpue97(48,(2:end)),'color','k','LineWidth',2,'Parent',h4); 
ylabel(h3,'Anomalies')
title('LME 48, SST-1, Det-3, Bio-2')
print('-dpng',[ppath 'ts_LME48_CPUE_sig_drivers_chl15_bestlags.png'])

%% LME 54
% GOOD ONE
% LME slag	sat	dlag	driver	flag	fdriver
% 54	3	Chl	1	    Det	    1	    Det

figure(54)
%plot with lags
h1 = axes('Units','pixels','Position',axesPosition,...
          'Color','w','XColor','k','YColor',cchl,...
          'XLim',xLimit,'NextPlot','add');
h2 = axes('Units','pixels','Position',axesPosition+yWidth.*[-1 0 1 0],...
          'Color','none','XColor','k','YColor',cdet,...
          'XLim',xLimit+[xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');
h3 = axes('Units','pixels','Position',axesPosition+yWidth.*[-2 0 2 0],...
          'Color','none','XColor','k','YColor','k',...
          'XLim',xLimit+[2*xOffset 0],'XTick',[],'XTickLabel',[],...
          'NextPlot','add');

line(dyr3,achl15(54,(1:(end-3))),'color',cchl,'LineWidth',2,'Parent',h1); hold on;
line(dyr1,adet15(54,(1:(end-1))),'color',cdet,'LineWidth',2,'Parent',h2); hold on;
xlabel(h1,'Year');

line(cyr1,aa_cpue97(54,(2:end)),'color','k','LineWidth',2,'Parent',h3); 
ylabel(h3,'Anomalies')
title('LME 54, Chl-3, Det-1')
print('-dpng',[ppath 'ts_LME54_CPUE_sig_drivers_chl15_bestlags.png'])


%% LME 58
% Pretty good
% LME slag	sat	dlag	driver	flag	fdriver
% 58	2	Chl	    0	-TB	    0	    -TB

figure(58)
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

line(dyr2,achl15(58,(1:(end-2))),'color',cchl,'LineWidth',2,'Parent',h1); hold on;
line(yyr,-1*atb15(58,:),'color',ctb,'LineWidth',2,'Parent',h2); hold on;
xlabel(h1,'Year');

line(yyr,aa_cpue97(58,:),'color','k','LineWidth',2,'Parent',h3); 
ylabel(h3,'Anomalies')
title('LME 58, Chl-2, TB0')
print('-dpng',[ppath 'ts_LME58_CPUE_sig_drivers_chl15_bestlags.png'])

%% LME 61
% MAYBE
% LME slag	sat	dlag	driver	flag	fdriver
% 61	3	SST	0	    ZmL	    0	    ZmL

figure(61)
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

line(dyr3,asst15(61,(1:(end-3))),'color',csst,'LineWidth',2,'Parent',h1); hold on;
line(yyr,azlos15(61,:),'color',czm,'LineWidth',2,'Parent',h2); 
xlabel(h1,'Year');

line(yyr,aa_cpue97(61,:),'color','k','LineWidth',2,'Parent',h3); 
ylabel(h3,'Anomalies')
title('LME 61, SST-3, ZmL0')
print('-dpng',[ppath 'ts_LME61_CPUE_sig_drivers_chl15_bestlags.png'])

