% Fish-MIP Phase 3a CPUE reconstruction
% functional type CPUE by lme
% calc anom ts
% Calculate anomaly time series for different ranges
% 1982-2015 

clear
close all

%% Fishing data
%fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/Fish-MIP/Phase3/fishing/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

% load([fpath 'FishMIP_Phase3a_LME_Catch_annual_1948-2015.mat']);
% load([fpath 'FishMIP_Phase3a_LME_Effort_annual_1961-2010.mat']);
load([fpath 'FishMIP_Phase3a_LME_CPUE_1948-2015_annual.mat']);

%% Next select specific yrs
cyr = 1982:2015; % subset effort years
[~,fidc] = intersect(yrs,cyr);

%% mean & std by lme

% Mean
cpue_f82_mean  = mean(cpue_f_ts(:,fidc),2,'omitnan');
cpue_p82_mean  = mean(cpue_p_ts(:,fidc),2,'omitnan');
cpue_d82_mean  = mean(cpue_d_ts(:,fidc),2,'omitnan');
cpue_a82_mean  = mean(cpue_a_ts(:,fidc),2,'omitnan');

% Std dev
cpue_f82_std  = std(cpue_f_ts(:,fidc),0,2,'omitnan');
cpue_p82_std  = std(cpue_p_ts(:,fidc),0,2,'omitnan');
cpue_d82_std  = std(cpue_d_ts(:,fidc),0,2,'omitnan');
cpue_a82_std  = std(cpue_a_ts(:,fidc),0,2,'omitnan');

% Coefficient of variance
cpue_f82_cv = cpue_f82_std ./ cpue_f82_mean;
cpue_p82_cv = cpue_p82_std ./ cpue_p82_mean;
cpue_d82_cv = cpue_d82_std ./ cpue_d82_mean;
cpue_a82_cv = cpue_a82_std ./ cpue_a82_mean;

%% ANOMALIES -------------------------------------------------

% remove linear trend
nyr = length(fidc);
F = zeros(66,nyr);
P = zeros(66,nyr);
D = zeros(66,nyr);
A = zeros(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62

    %TYPES
    xi = cpue_f_ts(i,fidc);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        F(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = cpue_p_ts(i,fidc);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        P(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = cpue_d_ts(i,fidc);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        D(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = cpue_a_ts(i,fidc);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        A(i,:) = dR;
    end
    clear R T t b m tH dR data

end

%% anomalies
af_cpue82 = F - mean(F,2,'omitnan');
ap_cpue82 = P - mean(P,2,'omitnan');
ad_cpue82 = D - mean(D,2,'omitnan');
aa_cpue82 = A - mean(A,2,'omitnan');

%% var of anomalies by lme
vf_cpue82 = var(af_cpue82,0,2,'omitnan');
vp_cpue82 = var(ap_cpue82,0,2,'omitnan');
vd_cpue82 = var(ad_cpue82,0,2,'omitnan');
va_cpue82 = var(aa_cpue82,0,2,'omitnan');

%% save
%units = 'MT?';
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
save([fpath 'FishMIP_Phase3a_LME_CPUE_1982-2015_interann_var.mat'],...
    'cpue_f82_std','cpue_p82_std','cpue_d82_std','cpue_a82_std',...
    'cpue_f82_mean','cpue_p82_mean','cpue_d82_mean','cpue_a82_mean',...
    'cpue_f_ts','cpue_p_ts','cpue_d_ts','cpue_a_ts',...
    'cpue_f82_cv','cpue_p82_cv','cpue_d82_cv','cpue_a82_cv');

%%
save([fpath 'FishMIP_Phase3a_LME_CPUE_1982-2015_ann_mean_anoms.mat'],...
    'af_cpue82','ap_cpue82','ad_cpue82','vf_cpue82','vp_cpue82','vd_cpue82','aa_cpue82','va_cpue82');
