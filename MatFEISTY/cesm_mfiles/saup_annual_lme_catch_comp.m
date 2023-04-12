% Comp new SAUP LME cathces to old

clear
close all

spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/SAUP/';
cpath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY_other/grid_cobalt/';
dp = '/Volumes/petrik-lab/feisty/NC/Matlab_new_size/';
pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';

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

% plot info
[ni,nj]=size(lon);
geolon_t = double(lon);
geolat_t = double(lat);
plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon]; %[-255 -60] = Pac

land=-999*ones(ni,nj);
land(ID)=NaN*ones(size(ID));

%% SAUP data
% use weighted catches
load([spath 'SAUP_LME_Catch_annual.mat'],'yr','totcatch','lme_catch',...
    'Flme_wcatch','Dlme_wcatch','Plme_wcatch');

Flme_catch_all = nansum(Flme_wcatch,3);
Plme_catch_all = nansum(Plme_wcatch,3);
Dlme_catch_all = nansum(Dlme_wcatch,3);

%1956-2005 SAUP average
id = find(yr>1955 & yr<=2005);

slme_mcatch = nanmean(lme_catch(id,:));
slme_mcatch = slme_mcatch';
Fslme_mcatch = nanmean(Flme_catch_all(id,:));
Fslme_mcatch = Fslme_mcatch';
Pslme_mcatch = nanmean(Plme_catch_all(id,:));
Pslme_mcatch = Pslme_mcatch';
Dslme_mcatch = nanmean(Dlme_catch_all(id,:));
Dslme_mcatch = Dslme_mcatch';

%% Find 10 years of max catch
slme_mcatch10 = NaN*ones(size(slme_mcatch));
Flme_mcatch10 = NaN*ones(size(slme_mcatch));
Plme_mcatch10 = NaN*ones(size(slme_mcatch));
Dlme_mcatch10 = NaN*ones(size(slme_mcatch));
%Top 10 yrs SAUP
for i=1:66
    [sort_lme_catch,ix] = sort(lme_catch(:,i),'descend');
    sort_Flme_catch = Flme_catch_all(ix,i);
    sort_Plme_catch = Plme_catch_all(ix,i);
    sort_Dlme_catch = Dlme_catch_all(ix,i);
    slme_mcatch10(i) = nanmean(sort_lme_catch(1:10));
    Flme_mcatch10(i) = nanmean(sort_Flme_catch(1:10));
    Plme_mcatch10(i) = nanmean(sort_Plme_catch(1:10));
    Dlme_mcatch10(i) = nanmean(sort_Dlme_catch(1:10));
end

% MT/km2
slme_mcatch10 = slme_mcatch10 ./ lme_area_km2;
Flme_mcatch10 = Flme_mcatch10 ./ lme_area_km2;
Plme_mcatch10 = Plme_mcatch10 ./ lme_area_km2;
Dlme_mcatch10 = Dlme_mcatch10 ./ lme_area_km2;

%% NEW forage data
% save([spath 'SAUP_LME_Catch_annual_newF.mat'],'Year','LMEID','FwCatch',...
%     'PwCatch','DwCatch');

load([spath 'SAUP_LME_Catch_annual_newF.mat'],'Year','LMEID','FwCatch',...
    'PwCatch','DwCatch');

%% arrange in matrix LME x Year
yrs_new = min(Year):max(Year); %same as yr

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

%% P & D are the same except for one LME
figure(1)
subplot(2,2,1)
plot(log10(slme_mcatch+1),log10(mTlme+1),'.')
title('All')

subplot(2,2,2)
plot(log10(Fslme_mcatch+1),log10(mFlme+1),'.')
title('F')

subplot(2,2,3)
plot(log10(Pslme_mcatch+1),log10(mPlme+1),'.')
title('P')

subplot(2,2,4)
plot(log10(Dslme_mcatch+1),log10(mDlme+1),'.')
title('D')

%% Find 10 years of max catch
Tlme_m10 = NaN*ones(size(slme_mcatch));
Flme_m10 = NaN*ones(size(slme_mcatch));
Plme_m10 = NaN*ones(size(slme_mcatch));
Dlme_m10 = NaN*ones(size(slme_mcatch));
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
Tlme_mc10 = Tlme_m10 ./ lme_area_km2;
Flme_mc10 = Flme_m10 ./ lme_area_km2;
Plme_mc10 = Plme_m10 ./ lme_area_km2;
Dlme_mc10 = Dlme_m10 ./ lme_area_km2;

%% b/c total is diff, P & D a little diff (only 2-4 noticable)
figure(2)
subplot(2,2,1)
plot(log10(slme_mcatch10+1),log10(Tlme_mc10+1),'.')
title('All')

subplot(2,2,2)
plot(log10(Flme_mcatch10+1),log10(Flme_mc10+1),'.')
title('F')

subplot(2,2,3)
plot(log10(Plme_mcatch10+1),log10(Plme_mc10+1),'.')
title('P')

subplot(2,2,4)
plot(log10(Dlme_mcatch10+1),log10(Dlme_mc10+1),'.')
title('D')

%% save new dataset
save([spath 'SAUP_LME_Catch_top10_Stock_newF.mat'],'Tlme_mc10','Flme_mc10',...
    'Plme_mc10','Dlme_mc10');




