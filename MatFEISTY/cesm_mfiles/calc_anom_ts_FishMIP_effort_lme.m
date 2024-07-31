% Fish-MIP Phase 3a effort reconstruction
% functional type by lme
% calc anom ts

clear
close all

%% Fishing data
%fpath='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/Fish-MIP/Phase3/fishing/';
fpath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

load([fpath 'FishMIP_Phase3a_LME_Effort_annual_1948-2010.mat']);

%% Cols: Year, LME, F, P, D
eyr = unique(Effort(:,1));
elid = unique(Effort(:,2));

ni = length(eyr);

%% Brute force reshape
% Effort
Effort(:,6) = sum(Effort(:,3:5),2,'omitnan');

% Brute force
efrt_f_ts  = zeros(66,ni);
efrt_p_ts  = zeros(66,ni);
efrt_d_ts  = zeros(66,ni);
efrt_a_ts  = zeros(66,ni);

for L=1:66
    for i = 1:ni
        Y = eyr(i);
        yid = find(Effort(:,1)==Y);
        lid = find(Effort(:,2)==L);
        id = intersect(yid,lid);
        if(~isempty(id))
            efrt_f_ts(L,i) = Effort(id,3);
            efrt_p_ts(L,i) = Effort(id,4);
            efrt_d_ts(L,i) = Effort(id,5);
            efrt_a_ts(L,i) = Effort(id,6);
        end
    end
end

%% Remove effort trend
nyr = ni;
eff_F = zeros(66,nyr);
eff_P = zeros(66,nyr);
eff_D = zeros(66,nyr);
eff_A = zeros(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62

    %TYPES
    xi = efrt_f_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        eff_F(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = efrt_p_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        eff_P(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = efrt_d_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        eff_D(i,:) = dR;
    end
    clear R T t b m tH dR data

    xi = efrt_a_ts(i,:);
    inan = ~isnan(xi);
    if(sum(inan)~=0)
        R = (xi(inan))';
        T = length(R);
        t = (1:T)';
        data = [t R];
        [m,b] = TheilSen(data);
        tH = m*t + b;
        dR = R - tH;
        eff_A(i,:) = dR;
    end
    clear R T t b m tH dR data

end

%% mean & std by lme

% Mean
eff_f_mean  = mean(eff_F,2,'omitnan');
eff_p_mean  = mean(eff_P,2,'omitnan');
eff_d_mean  = mean(eff_D,2,'omitnan');
eff_a_mean  = mean(eff_A,2,'omitnan');

% Std dev
eff_f_std  = std(eff_F,0,2,'omitnan');
eff_p_std  = std(eff_P,0,2,'omitnan');
eff_d_std  = std(eff_D,0,2,'omitnan');
eff_a_std  = std(eff_A,0,2,'omitnan');

% Coefficient of variance
eff_f_cv = eff_f_std ./ eff_f_mean;
eff_p_cv = eff_p_std ./ eff_p_mean;
eff_d_cv = eff_d_std ./ eff_d_mean;
eff_a_cv = eff_a_std ./ eff_a_mean;

%% ANOMALIES -------------------------------------------------

% remove linear trend
F = zeros(66,nyr);
P = zeros(66,nyr);
D = zeros(66,nyr);
A = zeros(66,nyr);

for i = 1:66 %in future should remove interior seas 23, 33, 62

    %TYPES
    xi = eff_F(i,:);
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

    xi = eff_P(i,:);
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

    xi = eff_D(i,:);
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

    xi = eff_A(i,:);
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
save([fpath 'FishMIP_Phase3a_LME_effort_1948-2010_interann_var.mat'],...
    'eff_f_std','eff_p_std','eff_d_std','eff_a_std',...
    'eff_f_mean','eff_p_mean','eff_d_mean','eff_a_mean',...
    'eff_F','eff_P','eff_D','eff_A',...
    'eff_f_cv','eff_p_cv','eff_d_cv','eff_a_cv');

%%
aall = aa;
vall = va;
save([fpath 'FishMIP_Phase3a_LME_effort_1948-2010_ann_mean_anoms.mat'],...
    'af','ap','ad','vf','vp','vd','aall','vall');
