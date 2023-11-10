% Fish-MIP Phase 3a catch & effort reconstruction
% functional type by lme
% regress catch on effort and remove the effect
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
eyr = unique(Effort(:,1));
elid = unique(Effort(:,2));

ni = length(eyr);

% Not even, can't reshape
LMECatch(:,6) = sum(LMECatch(:,3:5),2,'omitnan');

%% Brute force
lme_f_ts  = zeros(66,ni);
lme_p_ts  = zeros(66,ni);
lme_d_ts  = zeros(66,ni);
lme_a_ts  = zeros(66,ni);

for L=1:66
    for i = 1:length(yr)
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

%%
[yrs,yid] = intersect(yr,eyr);
lme_f_ts  = lme_f_ts(:,yid);
lme_p_ts  = lme_p_ts(:,yid);
lme_d_ts  = lme_d_ts(:,yid);
lme_a_ts  = lme_a_ts(:,yid);

%% Effort
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

%% Remove effort trend
nyr = ni;
cme_F = zeros(66,nyr);
cme_P = zeros(66,nyr);
cme_D = zeros(66,nyr);
cme_A = zeros(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62

    %TYPES
    xi = lme_f_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        t = eff_f_ts(i,:)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        cme_F(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = lme_p_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        t = eff_p_ts(i,:)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        cme_P(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = lme_d_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        t = eff_d_ts(i,:)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        cme_D(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = lme_a_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        t = eff_a_ts(i,:)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        cme_A(i,:) = dR;
    end
    clear R T t b m tH dR data

end


%% mean & std by lme

% Mean
cme_f_mean  = mean(cme_F,2,'omitnan');
cme_p_mean  = mean(cme_P,2,'omitnan');
cme_d_mean  = mean(cme_D,2,'omitnan');
cme_a_mean  = mean(cme_A,2,'omitnan');

% Std dev
cme_f_std  = std(cme_F,0,2,'omitnan');
cme_p_std  = std(cme_P,0,2,'omitnan');
cme_d_std  = std(cme_D,0,2,'omitnan');
cme_a_std  = std(cme_A,0,2,'omitnan');

% Coefficient of variance
cme_f_cv = cme_f_std ./ cme_f_mean;
cme_p_cv = cme_p_std ./ cme_p_mean;
cme_d_cv = cme_d_std ./ cme_d_mean;
cme_a_cv = cme_a_std ./ cme_a_mean;

%% ANOMALIES -------------------------------------------------

% remove linear trend
F = zeros(66,nyr);
P = zeros(66,nyr);
D = zeros(66,nyr);
A = zeros(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62

    %TYPES
    xi = cme_F(i,:);
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

    xi = cme_P(i,:);
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

    xi = cme_D(i,:);
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

    xi = cme_A(i,:);
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
save([fpath 'FishMIP_Phase3a_LME_catch_minus_effort_1961-2010_interann_var.mat'],...
    'cme_f_std','cme_p_std','cme_d_std','cme_a_std',...
    'cme_f_mean','cme_p_mean','cme_d_mean','cme_a_mean',...
    'cme_F','cme_P','cme_D','cme_A',...
    'cme_f_cv','cme_p_cv','cme_d_cv','cme_a_cv');

%%
aall = aa;
vall = va;
save([fpath 'FishMIP_Phase3a_LME_catch_minus_effort_1961-2010_ann_mean_anoms.mat'],...
    'af','ap','ad','vf','vp','vd','aall','vall');
