% Fish-MIP Phase 3a catch reconstruction
% functional type catch by lme
% calc anom ts

clear
close all

%% Fishing data
%fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

load([fpath 'FishMIP_Phase3a_LME_Catch_annual.mat']);


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

%% Brute force
lme_f_ts  = nan*ones(66,ni);
lme_p_ts  = nan*ones(66,ni);
lme_d_ts  = nan*ones(66,ni);

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
        end
    end
end

%% mean & std by lme

% Mean
lme_f_mean  = mean(lme_f_ts,2,'omitnan');
lme_p_mean  = mean(lme_p_ts,2,'omitnan');
lme_d_mean  = mean(lme_d_ts,2,'omitnan');

% Std dev
lme_f_std  = std(lme_f_ts,0,2,'omitnan');
lme_p_std  = std(lme_p_ts,0,2,'omitnan');
lme_d_std  = std(lme_d_ts,0,2,'omitnan');

% Coefficient of variance
lme_f_cv = lme_f_std ./ lme_f_mean;
lme_p_cv = lme_p_std ./ lme_p_mean;
lme_d_cv = lme_d_std ./ lme_d_mean;

%% ANOMALIES -------------------------------------------------

% FOSI is 1948 to 2015
% catches are 1961 to 2005

%% remove linear trend
F = NaN*ones(66,nyr);
P = NaN*ones(66,nyr);
D = NaN*ones(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62

    %TYPES
    xi = lme_f_ts(i,:);
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

    xi = lme_p_ts(i,:);
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

    xi = lme_d_ts(i,:);
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

    xi = lme_a_mean_ts(i,:);
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
units = 'per year';
%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
save([fpath 'FEISTY_FOSI_',mod,'lme_nu_interann_var.mat'],...
    'lme_f_std','lme_p_std','lme_d_std','lme_a_std',...
    'lme_f_mean','lme_p_mean','lme_d_mean','lme_a_mean',...
    'lme_f_ts','lme_p_ts','lme_d_ts','lme_a_mean_ts',...
    'lme_f_cv','lme_p_cv','lme_d_cv','lme_a_cv','units');

%%
save([fpath 'FEISTY_FOSI_',mod,'lme_nu_ann_mean_anoms.mat'],...
    'af','ap','ad','aa',...
    'vf','vp','vd','va','units');
