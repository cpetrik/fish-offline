% Fish-MIP Phase 3a catch reconstruction
% functional type catch by lme
% calc anom ts

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

%% 
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

%% mean & std by lme

% Mean
cpue_f_mean  = mean(cpue_f_ts,2,'omitnan');
cpue_p_mean  = mean(cpue_p_ts,2,'omitnan');
cpue_d_mean  = mean(cpue_d_ts,2,'omitnan');
cpue_a_mean  = mean(cpue_a_ts,2,'omitnan');

% Std dev
cpue_f_std  = std(cpue_f_ts,0,2,'omitnan');
cpue_p_std  = std(cpue_p_ts,0,2,'omitnan');
cpue_d_std  = std(cpue_d_ts,0,2,'omitnan');
cpue_a_std  = std(cpue_a_ts,0,2,'omitnan');

% Coefficient of variance
cpue_f_cv = cpue_f_std ./ cpue_f_mean;
cpue_p_cv = cpue_p_std ./ cpue_p_mean;
cpue_d_cv = cpue_d_std ./ cpue_d_mean;
cpue_a_cv = cpue_a_std ./ cpue_a_mean;

%% ANOMALIES -------------------------------------------------

% remove linear trend
nyr = ni;
F = zeros(66,nyr);
P = zeros(66,nyr);
D = zeros(66,nyr);
A = zeros(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62

    %TYPES
    xi = cpue_f_ts(i,:);
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

    xi = cpue_p_ts(i,:);
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

    xi = cpue_d_ts(i,:);
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

    xi = cpue_a_ts(i,:);
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
af = F - mean(F,2,'omitnan');
ap = P - mean(P,2,'omitnan');
ad = D - mean(D,2,'omitnan');
aa = A - mean(A,2,'omitnan');

%% var of anomalies by lme
vf = var(af,0,2,'omitnan');
vp = var(ap,0,2,'omitnan');
vd = var(ad,0,2,'omitnan');
va = var(aa,0,2,'omitnan');

%% save
%units = 'MT?';
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
save([fpath 'FishMIP_Phase3a_LME_CPUE_1961-2010_interann_var.mat'],...
    'cpue_f_std','cpue_p_std','cpue_d_std','cpue_a_std',...
    'cpue_f_mean','cpue_p_mean','cpue_d_mean','cpue_a_mean',...
    'cpue_f_ts','cpue_p_ts','cpue_d_ts','cpue_a_ts',...
    'cpue_f_cv','cpue_p_cv','cpue_d_cv','cpue_a_cv');

%%
aall = aa;
vall = va;
save([fpath 'FishMIP_Phase3a_LME_CPUE_1961-2010_ann_mean_anoms.mat'],...
    'af','ap','ad','vf','vp','vd','aall','vall');
