% Fish-MIP Phase 3a catch reconstruction
% functional type catch by lme
% calc anom ts
% Calculate anomaly time series for different ranges
% 1997-2015 

clear
close all

%% Fishing data
%fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/Fish-MIP/Phase3/fishing/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

load([fpath 'FishMIP_Phase3a_LME_Catch_annual_1948-2015.mat']);
load([fpath 'FishMIP_Phase3a_LME_Effort_annual_1948-2017.mat']);

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

cpue_f_ts  = lme_f_ts(:,yid) ./ eff_f_ts(:,yid);
cpue_p_ts  = lme_p_ts(:,yid) ./ eff_p_ts(:,yid);
cpue_d_ts  = lme_d_ts(:,yid) ./ eff_d_ts(:,yid);
cpue_a_ts  = lme_a_ts(:,yid) ./ eff_a_ts(:,yid);

cpue_f_ts(isnan(cpue_f_ts)) = 0;
cpue_p_ts(isnan(cpue_p_ts)) = 0;
cpue_d_ts(isnan(cpue_d_ts)) = 0;
cpue_a_ts(isnan(cpue_a_ts)) = 0;

cpue_f_ts(isinf(cpue_f_ts)) = 0;
cpue_p_ts(isinf(cpue_p_ts)) = 0;
cpue_d_ts(isinf(cpue_d_ts)) = 0;
cpue_a_ts(isinf(cpue_a_ts)) = 0;

%% save CPUE ts
save([fpath 'FishMIP_Phase3a_LME_CPUE_1948-2015_annual.mat'],...
    'cpue_f_ts','cpue_p_ts','cpue_d_ts','cpue_a_ts','yrs');

%% Next select specific yrs
cyr = 1997:2015; % subset effort years
[~,fidc] = intersect(eyr,cyr);

%% mean & std by lme

% Mean
cpue_f97_mean  = mean(cpue_f_ts(:,fidc),2,'omitnan');
cpue_p97_mean  = mean(cpue_p_ts(:,fidc),2,'omitnan');
cpue_d97_mean  = mean(cpue_d_ts(:,fidc),2,'omitnan');
cpue_a97_mean  = mean(cpue_a_ts(:,fidc),2,'omitnan');

% Std dev
cpue_f97_std  = std(cpue_f_ts(:,fidc),0,2,'omitnan');
cpue_p97_std  = std(cpue_p_ts(:,fidc),0,2,'omitnan');
cpue_d97_std  = std(cpue_d_ts(:,fidc),0,2,'omitnan');
cpue_a97_std  = std(cpue_a_ts(:,fidc),0,2,'omitnan');

% Coefficient of variance
cpue_f97_cv = cpue_f97_std ./ cpue_f97_mean;
cpue_p97_cv = cpue_p97_std ./ cpue_p97_mean;
cpue_d97_cv = cpue_d97_std ./ cpue_d97_mean;
cpue_a97_cv = cpue_a97_std ./ cpue_a97_mean;

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
af_cpue97 = F - mean(F,2,'omitnan');
ap_cpue97 = P - mean(P,2,'omitnan');
ad_cpue97 = D - mean(D,2,'omitnan');
aa_cpue97 = A - mean(A,2,'omitnan');

%% var of anomalies by lme
vf_cpue97 = var(af_cpue97,0,2,'omitnan');
vp_cpue97 = var(ap_cpue97,0,2,'omitnan');
vd_cpue97 = var(ad_cpue97,0,2,'omitnan');
va_cpue97 = var(aa_cpue97,0,2,'omitnan');

%% save
%units = 'MT?';
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
save([fpath 'FishMIP_Phase3a_LME_CPUE_1997-2015_interann_var.mat'],...
    'cpue_f97_std','cpue_p97_std','cpue_d97_std','cpue_a97_std',...
    'cpue_f97_mean','cpue_p97_mean','cpue_d97_mean','cpue_a97_mean',...
    'cpue_f_ts','cpue_p_ts','cpue_d_ts','cpue_a_ts',...
    'cpue_f97_cv','cpue_p97_cv','cpue_d97_cv','cpue_a97_cv');

%%
save([fpath 'FishMIP_Phase3a_LME_CPUE_1997-2015_ann_mean_anoms.mat'],...
    'af_cpue97','ap_cpue97','ad_cpue97','vf_cpue97','vp_cpue97','vd_cpue97','aa_cpue97','va_cpue97');
