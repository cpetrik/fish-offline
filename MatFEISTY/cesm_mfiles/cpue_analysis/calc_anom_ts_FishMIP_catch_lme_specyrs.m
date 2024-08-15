% Fish-MIP Phase 3a catch reconstruction
% functional type catch by lme
% calc anom ts
% Calculate anomaly time series for different ranges
% 1997-2015

clear
close all

%% Fishing data
%fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

load([fpath 'FishMIP_Phase3a_LME_Catch_annual_1948-2015.mat']);

%% Cols: Year, LME, F, P, D
yr = unique(LMECatch(:,1));
lid = unique(LMECatch(:,2));

ni = length(yr);
nj = length(lid);

% mf = LMECatch(:,[1:3]);
% lp = LMECatch(:,[1:2,4]);
% ld = LMECatch(:,[1:2,5]);
%
% % Not even, can't reshape
% mf1 = reshape(mf(:,1),ni,nj);
% mf2 = reshape(mf(:,2),ni,nj);
% mf3 = reshape(mf(:,3),ni,nj);

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

%% Next select specific yrs
cyr = 1997:2015; % subset effort years
[~,fidc] = intersect(yr,cyr);

%% mean & std by lme

% Mean
lme_f97_mean  = mean(lme_f_ts(:,fidc),2,'omitnan');
lme_p97_mean  = mean(lme_p_ts(:,fidc),2,'omitnan');
lme_d97_mean  = mean(lme_d_ts(:,fidc),2,'omitnan');
lme_a97_mean  = mean(lme_a_ts(:,fidc),2,'omitnan');

% Std dev
lme_f97_std  = std(lme_f_ts(:,fidc),0,2,'omitnan');
lme_p97_std  = std(lme_p_ts(:,fidc),0,2,'omitnan');
lme_d97_std  = std(lme_d_ts(:,fidc),0,2,'omitnan');
lme_a97_std  = std(lme_a_ts(:,fidc),0,2,'omitnan');

% Coefficient of variance
lme_f97_cv = lme_f97_std ./ lme_f97_mean;
lme_p97_cv = lme_p97_std ./ lme_p97_mean;
lme_d97_cv = lme_d97_std ./ lme_d97_mean;
lme_a97_cv = lme_a97_std ./ lme_a97_mean;

%% ANOMALIES -------------------------------------------------

% FOSI is 1948 to 2015
% catches are 1948 to 2015

%% remove linear trend
nyr = length(cyr);
F = zeros(66,nyr);
P = zeros(66,nyr);
D = zeros(66,nyr);
A = zeros(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62

    %TYPES
    xi = lme_f_ts(i,fidc);
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

    xi = lme_p_ts(i,fidc);
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

    xi = lme_d_ts(i,fidc);
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

    xi = lme_a_ts(i,fidc);
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
af_catch97 = F - mean(F,2,'omitnan');
ap_catch97 = P - mean(P,2,'omitnan');
ad_catch97 = D - mean(D,2,'omitnan');
aa_catch97 = A - mean(A,2,'omitnan');

%% var of anomalies by lme
vf_catch97 = var(af_catch97,0,2,'omitnan');
vp_catch97 = var(ap_catch97,0,2,'omitnan');
vd_catch97 = var(ad_catch97,0,2,'omitnan');
va_catch97 = var(aa_catch97,0,2,'omitnan');

%% save
units = 'MT?';
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
save([fpath 'FishMIP_Phase3a_LME_Catch_1997-2015_interann_var.mat'],...
    'lme_f97_std','lme_p97_std','lme_d97_std','lme_a97_std',...
    'lme_f97_mean','lme_p97_mean','lme_d97_mean','lme_a97_mean',...
    'lme_f_ts','lme_p_ts','lme_d_ts','lme_a_ts',...
    'lme_f97_cv','lme_p97_cv','lme_d97_cv','lme_a97_cv','units');

%%
save([fpath 'FishMIP_Phase3a_LME_Catch_1997-2015_ann_mean_anoms.mat'],...
    'af_catch97','ap_catch97','ad_catch97','vf_catch97','vp_catch97','vd_catch97','aa_catch97','va_catch97','units');
