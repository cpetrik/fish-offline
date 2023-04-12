% SAUP new F (squids) top 10 calc

clear
close all

spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/SAUP/';
cpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/grid_cobalt/';
dp = '/Volumes/petrik-lab/feisty/NC/Matlab_new_size/';

%% ESM2.4 grid info
Pdir = '/Volumes/petrik-lab/Feisty/GCM_Data/ESM26_hist/';
load([Pdir 'ESM26_1deg_5yr_clim_191_195_gridspec.mat']);
load([cpath 'esm26_lme_mask_onedeg_SAU_66.mat']);
load([cpath 'esm26_area_1deg.mat']);
load([cpath 'LME_clim_temp.mat']);

AREA_OCN = max(area,1);

%%
lpath = '/Volumes/petrik-lab/Feisty/NC/Matlab_new_size/Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100/Climatology/';
load([lpath 'LME_clim_fished_All_fish03_Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100.mat'],'lme_area');
lme_area_km2 = lme_area * 1e-6;

%% NEW forage data

load([spath 'SAUP_LME_Catch_annual_newF.mat'],'Year','LMEID','FwCatch',...
    'PwCatch','DwCatch');

%% arrange in matrix LME x Year
yr = min(Year):max(Year); 

FCmat = nan*ones(length(yr),66);
PCmat = nan*ones(length(yr),66);
DCmat = nan*ones(length(yr),66);

for t=1:length(yr)
    for l=1:66
        tid = find(Year==yr(t));
        lid = find(LMEID==l);
        id = intersect(tid,lid);

        if(~isempty(id))
            FCmat(t,l) = FwCatch(id);
            PCmat(t,l) = PwCatch(id);
            DCmat(t,l) = DwCatch(id);
        end
    end
end

TCmat = FCmat + PCmat + DCmat;

%% 1956-2005 SAUP average
id = find(yr>1955 & yr<=2005);

mTlme = nanmean(TCmat(id,:))';
mFlme = nanmean(FCmat(id,:))';
mPlme = nanmean(PCmat(id,:))';
mDlme = nanmean(DCmat(id,:))';

%% Find 10 years of max catch
Tlme_m10 = NaN*ones(size(mTlme));
Flme_m10 = NaN*ones(size(mTlme));
Plme_m10 = NaN*ones(size(mTlme));
Dlme_m10 = NaN*ones(size(mTlme));
%Top 10 yrs SAUP
for i=1:66
    [sort_Tcatch,ix] = sort(TCmat(:,i),'descend');
    sort_Fcatch = FCmat(ix,i);
    sort_Pcatch = PCmat(ix,i);
    sort_Dcatch = DCmat(ix,i);
    Tlme_m10(i) = nanmean(sort_Tcatch(1:10));
    Flme_m10(i) = nanmean(sort_Fcatch(1:10));
    Plme_m10(i) = nanmean(sort_Pcatch(1:10));
    Dlme_m10(i) = nanmean(sort_Dcatch(1:10));
end

%% MT/km2
slme_mcatch10 = Tlme_m10 ./ lme_area_km2;
Flme_mcatch10 = Flme_m10 ./ lme_area_km2;
Plme_mcatch10 = Plme_m10 ./ lme_area_km2;
Dlme_mcatch10 = Dlme_m10 ./ lme_area_km2;


%% save new dataset
units = 'tonnes/km2';
save([spath 'SAUP_LME_Catch_top10_Stock_newF.mat'],'slme_mcatch10','Flme_mcatch10',...
    'Plme_mcatch10','Dlme_mcatch10','units');




