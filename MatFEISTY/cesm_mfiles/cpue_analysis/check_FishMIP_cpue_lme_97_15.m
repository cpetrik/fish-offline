% Fish-MIP Phase 3a CPUE reconstruction 1861-2017
% functional type CPUE by lme
% Julia asked to check if <=2010 is the same as old file to 2010

clear
close all

%% Fishing data
%fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/Fish-MIP/Phase3/fishing/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

load([fpath 'FishMIP_Phase3a_LME_Catch_annual_1948-2015.mat']);
load([fpath 'FishMIP_Phase3a_LME_Effort_annual_1961-2010.mat']);

%% Cols: Year, LME, F, P, D
yr = unique(LMECatch(:,1));
lid = unique(LMECatch(:,2));

ni = length(yr);
nj = length(lid);

% Not even, can't reshape
LMECatch(:,6) = sum(LMECatch(:,3:5),2,'omitnan');

%% Brute force
lme_f_ts  = zeros(66,ni);
lme_p_ts  = zeros(66,ni);
lme_d_ts  = zeros(66,ni);
lme_a_ts  = zeros(66,ni);

for L=1:66
    for i = 1:ni
        Y = yr(i);
        yid = find(LMECatch(:,1)==Y);
        lid = find(LMECatch(:,2)==L);
        id = intersect(yid,lid);
        if(~isempty(id))
            lme_f_ts(L,i) = LMECatch(id,3);
            lme_p_ts(L,i) = LMECatch(id,4);
            lme_d_ts(L,i) = LMECatch(id,5);
            lme_a_ts(L,i) = LMECatch(id,6);
        end
    end
end

%% Effort
eyr = unique(Effort(:,1));
elid = unique(Effort(:,2));

ni = length(eyr);

Effort(:,6) = sum(Effort(:,3:5),2,'omitnan');

% Brute force
eff_f_ts  = zeros(66,ni);
eff_p_ts  = zeros(66,ni);
eff_d_ts  = zeros(66,ni);
eff_a_ts  = zeros(66,ni);

for L=1:66
    for i = 1:ni
        Y = eyr(i);
        yid = find(Effort(:,1)==Y);
        lid = find(Effort(:,2)==L);
        id = intersect(yid,lid);
        if(~isempty(id))
            eff_f_ts(L,i) = Effort(id,3);
            eff_p_ts(L,i) = Effort(id,4);
            eff_d_ts(L,i) = Effort(id,5);
            eff_a_ts(L,i) = Effort(id,6);
        end
    end
end

%% Calc Catch/Effort
% FOSI is 1948 to 2015
% catches are 1948 to 2015
% effort is 1961 to 2010
[yrs,yid] = intersect(yr,eyr);

cpue_f_ts  = lme_f_ts(:,yid) ./ eff_f_ts;
cpue_p_ts  = lme_p_ts(:,yid) ./ eff_p_ts;
cpue_d_ts  = lme_d_ts(:,yid) ./ eff_d_ts;
cpue_a_ts  = lme_a_ts(:,yid) ./ eff_a_ts;

cpue_f_ts(isnan(cpue_f_ts)) = 0;
cpue_p_ts(isnan(cpue_p_ts)) = 0;
cpue_d_ts(isnan(cpue_d_ts)) = 0;
cpue_a_ts(isnan(cpue_a_ts)) = 0;

cpue_f_ts(isinf(cpue_f_ts)) = 0;
cpue_p_ts(isinf(cpue_p_ts)) = 0;
cpue_d_ts(isinf(cpue_d_ts)) = 0;
cpue_a_ts(isinf(cpue_a_ts)) = 0;

%%
cpue_f_ts10 = cpue_f_ts;
cpue_p_ts10 = cpue_p_ts;
cpue_d_ts10 = cpue_d_ts;
cpue_a_ts10 = cpue_a_ts;

clear cpue_p_ts cpue_f_ts cpue_d_ts cpue_a_ts

%% 2017
load([fpath 'FishMIP_Phase3a_LME_CPUE_1948-2015_annual.mat']);

%% sum all LMEs
F10 = sum(cpue_f_ts10);
P10 = sum(cpue_p_ts10);
D10 = sum(cpue_d_ts10);
A10 = sum(cpue_a_ts10);

F15 = sum(cpue_f_ts);
P15 = sum(cpue_p_ts);
D15 = sum(cpue_d_ts);
A15 = sum(cpue_a_ts);

%% Plot time series

figure(1)
plot(eyr,F10,'k'); hold on
plot(yrs,F15,'--b');
title('Forage')
xlim([1960 2015])

figure(2)
plot(eyr,P10,'k'); hold on
plot(yrs,P15,'--b');
title('Lg Pel')
xlim([1960 2015])

figure(3)
plot(eyr,D10,'k'); hold on
plot(yrs,D15,'--b');
title('Demersal')
xlim([1960 2015])

figure(4)
plot(eyr,A10,'k'); hold on
plot(yrs,A15,'--b');
title('All')
xlim([1960 2015])

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',...
    cfile,'/corrs_cpue/'];

figure(5)
subplot(2,2,1)
plot(eyr,F10,'k','LineWidth',1.5); hold on
plot(yrs,F15,'--b','LineWidth',2);
title('Forage')
xlim([1948 2017])

subplot(2,2,2)
plot(eyr,P10,'k','LineWidth',1.5); hold on
plot(yrs,P15,'--b','LineWidth',2);
title('Lg Pel')
xlim([1948 2017])

subplot(2,2,3)
plot(eyr,D10,'k','LineWidth',1.5); hold on
plot(yrs,D15,'--b','LineWidth',2);
title('Demersal')
xlim([1948 2017])

subplot(2,2,4)
plot(eyr,A10,'k','LineWidth',1.5); hold on
plot(yrs,A15,'--b','LineWidth',2);
title('All')
xlim([1948 2017])
legend('2010','2017')
print('-dpng',[ppath 'Timeseries_effort_LME_1948_2010vs2017.png'])

%%
figure(6)
subplot(2,2,1)
plot(eyr,F10,'k','LineWidth',1.5); hold on
plot(yrs,F15,'--b','LineWidth',2);
title('Forage')
xlim([1961 2017])

subplot(2,2,2)
plot(eyr,P10,'k','LineWidth',1.5); hold on
plot(yrs,P15,'--b','LineWidth',2);
title('Lg Pel')
xlim([1961 2017])

subplot(2,2,3)
plot(eyr,D10,'k','LineWidth',1.5); hold on
plot(yrs,D15,'--b','LineWidth',2);
title('Demersal')
xlim([1961 2017])

subplot(2,2,4)
plot(eyr,A10,'k','LineWidth',1.5); hold on
plot(yrs,A15,'--b','LineWidth',2);
title('All')
xlim([1961 2017])
legend('2010','2017')
print('-dpng',[ppath 'Timeseries_effort_LME_1961_2010vs2017.png'])

%%
figure(7)
subplot(2,1,1)
plot(eyr,F10,'k','LineWidth',1.5); hold on
plot(yrs,F15,'--b','LineWidth',2);
title('Forage')
xlim([1961 2017])
ylim([0 15])

subplot(2,1,2)
plot(eyr,P10,'k','LineWidth',1.5); hold on
plot(yrs,P15,'--b','LineWidth',2);
title('Lg Pel')
xlim([1961 2017])
ylim([0 15])
legend('2010','2017')
print('-dpng',[ppath 'Timeseries_effort_LME_FPzoom_1961_2010vs2017.png'])
