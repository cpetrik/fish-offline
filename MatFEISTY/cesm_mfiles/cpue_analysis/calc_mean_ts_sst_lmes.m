% Satellite SST data
% LME means


clear 
close all

tpath = '/Volumes/petrik-lab/Feisty/Obs_data/OISST/';

%%
load([tpath 'oisst.annual.mean.1deg_1982_2020.mat']);
tyr = yrs;
clear yrs

%% LME data for regular 1 deg grid
opath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/OneDeg/';

load([opath 'lme_gfdl-mom6-cobalt2_onedeg_temporary.mat'])
load([opath 'cellarea_onedeg.mat']) %area in m2
tlme = fliplr(tlme);
cell_area = fliplr(cell_area);

%% LME only cells
[ni,nj,nts] = size(tos);
sst2 = reshape(tos,ni*nj,nts);

oid = find(~isnan(tlme(:)));
olme = tlme(oid);
sst2 = sst2(oid,:);
area_vec = cell_area(oid);

%% First create area-weighted mean time series for each LME
lme_sst_ts = NaN*ones(66,length(tyr));

for L=1:66
    lid = find(olme==L);
    
    lme_sst_ts(L,:)  = (sum(sst2(lid,:).*area_vec(lid,:),1,'omitnan')) ./ sum(area_vec(lid,:),'omitnan');
    
end


%% save
save([tpath 'lme_satellite_sst_interann_means_1982_2020.mat'],...
    'lme_sst_ts','tyr');

