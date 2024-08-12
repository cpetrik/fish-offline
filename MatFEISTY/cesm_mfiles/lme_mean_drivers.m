% Make table of mean drivers by LME
% Satellite, FOSI inputs, effort
% Without linear trend removed
% Area-weighted

clear
close all

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
cpath='/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
epath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');
load([cpath 'LME-mask-POP_gx1v6.mat'],'lme_mask','TAREA');

AREA_OCN = TAREA * 1e-4; %grid cell area in m2
tlme = double(lme_mask);
tlme(tlme<0) = nan;
olme = tlme(GRD.ID);

%% Sat
load([fpath 'lme_satellite_sst_chl_interann_var.mat'],...
    'lme_sst_mean','lme_chl_mean',...
    'lme_sst_ts','lme_chl_ts','cyr','tyr');

%Means over whole t.s.

%% FOSI inputs
% load([cpath 'CESM_FOSI_v15_interann_mean_forcings_anom.mat'],...
%     'tp','tb','dety','zlosy') 
% 
% [ni,nj,nyr] = size(dety);
fyr = 1948:2015;

%% lme area-weighted means
% lme_tp_mean = NaN*ones(66,nyr);
% lme_tb_mean = NaN*ones(66,nyr);
% lme_dety_mean = NaN*ones(66,nyr);
% lme_mzly_mean = NaN*ones(66,nyr);
% 
% %vectorize
% tp2 = reshape(tp,ni*nj,nyr);
% tb2 = reshape(tb,ni*nj,nyr);
% dety = reshape(dety,ni*nj,nyr);
% mzly = reshape(zlosy,ni*nj,nyr);
% area_vec = reshape(AREA_OCN,ni*nj,1);
% area_vec = repmat(area_vec,1,nyr);
% 
% %% Area-weighted LME t.s.
% for L=1:66
%     lid = find(tlme==L);
%     if ~isempty(lid)
% 
%         lme_tp_mean(L,:)  = (sum(tp2(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
%         lme_tb_mean(L,:)  = (sum(tb2(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
%         lme_dety_mean(L,:) = (sum(dety(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
%         lme_mzly_mean(L,:) = (sum(mzly(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
% 
%     end
% end
% 
% %% Save
% save([cpath 'CESM_FOSI_v15_interann_mean_forcings_anom.mat'],...
%     'lme_tp_mean','lme_tb_mean','lme_dety_mean','lme_mzly_mean','-append') 

load([cpath 'CESM_FOSI_v15_interann_mean_forcings_anom.mat'],...
    'lme_tp_mean','lme_tb_mean','lme_dety_mean','lme_mzly_mean') 

%% Effort
load([epath 'FishMIP_Phase3a_LME_Effort_annual_1948-2010.mat']);

% Cols: Year, LME, F, P, D
eyr = unique(Effort(:,1));
elid = unique(Effort(:,2));

ni = length(eyr);

%% Brute force reshape
% Effort
% Effort(:,6) = sum(Effort(:,3:5),2,'omitnan');
% 
% % Brute force
% efrt_f_ts  = zeros(66,ni);
% efrt_p_ts  = zeros(66,ni);
% efrt_d_ts  = zeros(66,ni);
% efrt_a_ts  = zeros(66,ni);
% 
% for L=1:66
%     for i = 1:ni
%         Y = eyr(i);
%         yid = find(Effort(:,1)==Y);
%         lid = find(Effort(:,2)==L);
%         id = intersect(yid,lid);
%         if(~isempty(id))
%             efrt_f_ts(L,i) = Effort(id,3);
%             efrt_p_ts(L,i) = Effort(id,4);
%             efrt_d_ts(L,i) = Effort(id,5);
%             efrt_a_ts(L,i) = Effort(id,6);
%         end
%     end
% end
% 
% %%
% save([epath 'FishMIP_Phase3a_LME_Effort_annual_1948-2010.mat'],...
%     'efrt_f_ts','efrt_p_ts','efrt_d_ts','efrt_a_ts',...
%     'eyr','elid','-append') 

%% Means over whole time series avail
% add ZmLoss and Det to make total secondary prod
sprod = lme_mzly_mean + lme_dety_mean;
zdrat = lme_mzly_mean ./ lme_dety_mean;

ctex = {'LME','SST','TP','TB','Chl','2ndProd','ZDratio',...
    'Aeffort','Feffort','Peffort','Deffort'};

mall(:,1) = [1:66]';
mall(:,2) = lme_sst_mean;
mall(:,3) = mean(lme_tp_mean,2);
mall(:,4) = mean(lme_tb_mean,2);
mall(:,5) = lme_chl_mean;
mall(:,6) = mean(sprod,2);
mall(:,7) = mean(zdrat,2);
mall(:,8) = mean(efrt_a_ts,2);
mall(:,9) = mean(efrt_f_ts,2);
mall(:,10) = mean(efrt_p_ts,2);
mall(:,11) = mean(efrt_d_ts,2);

%% Means over 1997-2010
yr = 1997:2010;
[~,cid] = intersect(cyr,yr);
[~,sid] = intersect(tyr,yr);
[~,eid] = intersect(eyr,yr);
[~,fid] = intersect(fyr,yr);

mlim(:,1) = [1:66]';
mlim(:,2) = mean(lme_sst_ts(:,sid),2);
mlim(:,3) = mean(lme_tp_mean(:,fid),2);
mlim(:,4) = mean(lme_tb_mean(:,fid),2);
mlim(:,5) = mean(lme_chl_ts(:,cid),2);
mlim(:,6) = mean(sprod(:,fid),2);
mlim(:,7) = mean(zdrat(:,fid),2);
mlim(:,8) = mean(efrt_a_ts(:,eid),2);
mlim(:,9) = mean(efrt_f_ts(:,eid),2);
mlim(:,10) = mean(efrt_p_ts(:,eid),2);
mlim(:,11) = mean(efrt_d_ts(:,eid),2);

%%
TabAll = array2table(mall,"VariableNames",ctex);

TabLim = array2table(mlim,"VariableNames",ctex);

%%
writetable(TabAll,[fpath,'LME_means_sat_driver_effort_fullts.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(TabLim,[fpath,'LME_means_sat_driver_effort_satys.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([fpath,'LME_means_sat_driver_effort.mat'],...
    'mall','mlim','TabAll','TabLim');





