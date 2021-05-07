% Calc seasonal climatology - Step 1 of powerspec slope
% ln-transform abund
% Spatially average across biomes
% Then create anomaly ts - Step 2 of powerspec slope

% Steps
% 1. Calc seasonal climatology
% 2. Create anomaly time series by removing climatol
% 3. Calc Theil-Sen linear trend
% 4. Remove trend if significant
% 5. Calc spectral slopes

clear all
close all

cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/grid_cobalt/';
load([cpath 'hindcast_gridspec.mat'],'geolon_t','geolat_t');
grid = csvread([cpath 'grid_csv.csv']);

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/Matlab_New_sizes/';
ppath = [pp cfile '/'];

%% Biomes 
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';
load([cpath 'COBALT_hist_biomes_1950_2005.mat']);

biome = {'LC','ECCS','ECSS','Coast'};

lid = find(biome4_hist==1);
uid = find(biome4_hist==2);
sid = find(biome4_hist==3);
cid = find(biome4_hist==4);


%% Preindust 400 yr output
load([fpath 'Monthly_Pre400_' cfile '.mat']);

%%
All = sp_bio+sf_bio+sd_bio+mp_bio+mf_bio+md_bio+lp_bio+ld_bio;
AllF = sf_bio+mf_bio;
AllP = sp_bio+mp_bio+lp_bio;
AllD = sd_bio+md_bio+ld_bio;
AllS = sp_bio+sf_bio+sd_bio;
AllM = mp_bio+mf_bio+md_bio;
AllL = lp_bio+ld_bio;

%% Biome means
% [ni,nj,nt] = size(All);
% sf_bio = reshape(sf_bio,ni*nj,nt);
% sp_bio = reshape(sp_bio,ni*nj,nt);
% sd_bio = reshape(sd_bio,ni*nj,nt);
% mf_bio = reshape(mf_bio,ni*nj,nt);
% mp_bio = reshape(mp_bio,ni*nj,nt);
% md_bio = reshape(md_bio,ni*nj,nt);
% lp_bio = reshape(lp_bio,ni*nj,nt);
% ld_bio = reshape(ld_bio,ni*nj,nt);
% b_bio  = reshape(b_bio,ni*nj,nt);
% 
% All  = reshape(All,ni*nj,nt);
% AllF = reshape(AllF,ni*nj,nt);
% AllP = reshape(AllP,ni*nj,nt);
% AllD = reshape(AllD,ni*nj,nt);
% AllS = reshape(AllS,ni*nj,nt);
% AllM = reshape(AllM,ni*nj,nt);
% AllL = reshape(AllL,ni*nj,nt);

[lg,id1] = intersect(grid(:,1),lid);
[ug,id2] = intersect(grid(:,1),uid);
[sg,id3] = intersect(grid(:,1),sid);
[cg,id4] = intersect(grid(:,1),cid);

sf_biome(1,:) = log(nanmean(sf_bio(id1,:),1));
sf_biome(2,:) = log(nanmean(sf_bio(id2,:),1));
sf_biome(3,:) = log(nanmean(sf_bio(id3,:),1));
sf_biome(4,:) = log(nanmean(sf_bio(id4,:),1));

sp_biome(1,:) = log(nanmean(sp_bio(id1,:),1));
sp_biome(2,:) = log(nanmean(sp_bio(id2,:),1));
sp_biome(3,:) = log(nanmean(sp_bio(id3,:),1));
sp_biome(4,:) = log(nanmean(sp_bio(id4,:),1));

sd_biome(1,:) = log(nanmean(sd_bio(id1,:),1));
sd_biome(2,:) = log(nanmean(sd_bio(id2,:),1));
sd_biome(3,:) = log(nanmean(sd_bio(id3,:),1));
sd_biome(4,:) = log(nanmean(sd_bio(id4,:),1));

mf_biome(1,:) = log(nanmean(mf_bio(id1,:),1));
mf_biome(2,:) = log(nanmean(mf_bio(id2,:),1));
mf_biome(3,:) = log(nanmean(mf_bio(id3,:),1));
mf_biome(4,:) = log(nanmean(mf_bio(id4,:),1));

mp_biome(1,:) = log(nanmean(mp_bio(id1,:),1));
mp_biome(2,:) = log(nanmean(mp_bio(id2,:),1));
mp_biome(3,:) = log(nanmean(mp_bio(id3,:),1));
mp_biome(4,:) = log(nanmean(mp_bio(id4,:),1));

md_biome(1,:) = log(nanmean(md_bio(id1,:),1));
md_biome(2,:) = log(nanmean(md_bio(id2,:),1));
md_biome(3,:) = log(nanmean(md_bio(id3,:),1));
md_biome(4,:) = log(nanmean(md_bio(id4,:),1));

lp_biome(1,:) = log(nanmean(lp_bio(id1,:),1));
lp_biome(2,:) = log(nanmean(lp_bio(id2,:),1));
lp_biome(3,:) = log(nanmean(lp_bio(id3,:),1));
lp_biome(4,:) = log(nanmean(lp_bio(id4,:),1));

ld_biome(1,:) = log(nanmean(ld_bio(id1,:),1));
ld_biome(2,:) = log(nanmean(ld_bio(id2,:),1));
ld_biome(3,:) = log(nanmean(ld_bio(id3,:),1));
ld_biome(4,:) = log(nanmean(ld_bio(id4,:),1));

b_biome(1,:) = log(nanmean(b_bio(id1,:),1));
b_biome(2,:) = log(nanmean(b_bio(id2,:),1));
b_biome(3,:) = log(nanmean(b_bio(id3,:),1));
b_biome(4,:) = log(nanmean(b_bio(id4,:),1));

All_biome(1,:) = log(nanmean(All(id1,:),1));
All_biome(2,:) = log(nanmean(All(id2,:),1));
All_biome(3,:) = log(nanmean(All(id3,:),1));
All_biome(4,:) = log(nanmean(All(id4,:),1));

AllF_biome(1,:) = log(nanmean(AllF(id1,:),1));
AllF_biome(2,:) = log(nanmean(AllF(id2,:),1));
AllF_biome(3,:) = log(nanmean(AllF(id3,:),1));
AllF_biome(4,:) = log(nanmean(AllF(id4,:),1));

AllP_biome(1,:) = log(nanmean(AllP(id1,:),1));
AllP_biome(2,:) = log(nanmean(AllP(id2,:),1));
AllP_biome(3,:) = log(nanmean(AllP(id3,:),1));
AllP_biome(4,:) = log(nanmean(AllP(id4,:),1));

AllD_biome(1,:) = log(nanmean(AllD(id1,:),1));
AllD_biome(2,:) = log(nanmean(AllD(id2,:),1));
AllD_biome(3,:) = log(nanmean(AllD(id3,:),1));
AllD_biome(4,:) = log(nanmean(AllD(id4,:),1));

AllS_biome(1,:) = log(nanmean(AllS(id1,:),1));
AllS_biome(2,:) = log(nanmean(AllS(id2,:),1));
AllS_biome(3,:) = log(nanmean(AllS(id3,:),1));
AllS_biome(4,:) = log(nanmean(AllS(id4,:),1));

AllM_biome(1,:) = log(nanmean(AllM(id1,:),1));
AllM_biome(2,:) = log(nanmean(AllM(id2,:),1));
AllM_biome(3,:) = log(nanmean(AllM(id3,:),1));
AllM_biome(4,:) = log(nanmean(AllM(id4,:),1));

AllL_biome(1,:) = log(nanmean(AllL(id1,:),1));
AllL_biome(2,:) = log(nanmean(AllL(id2,:),1));
AllL_biome(3,:) = log(nanmean(AllL(id3,:),1));
AllL_biome(4,:) = log(nanmean(AllL(id4,:),1));

%% quick plot
yid = (time/12)+1000;
nmo = length(yid);
nyr = length(yid)/12;

figure
subplot(2,2,1)
plot(AllF_biome(1,:))
subplot(2,2,2)
plot(AllF_biome(2,:))
subplot(2,2,3)
plot(AllF_biome(3,:))
subplot(2,2,4)
plot(AllF_biome(4,:))

figure
subplot(2,2,1)
plot(AllP_biome(1,:))
subplot(2,2,2)
plot(AllP_biome(2,:))
subplot(2,2,3)
plot(AllP_biome(3,:))
subplot(2,2,4)
plot(AllP_biome(4,:))

figure
subplot(2,2,1)
plot(yid,AllF_biome)
subplot(2,2,2)
plot(yid,AllP_biome)
subplot(2,2,3)
plot(yid,AllD_biome)
subplot(2,2,4)
plot(yid,b_biome)

%% By hand
all_clim = nan*ones(4,12);
sf_clim = all_clim;
sp_clim = all_clim;
sd_clim = all_clim;
mf_clim = all_clim;
mp_clim = all_clim;
md_clim = all_clim;
lp_clim = all_clim;
ld_clim = all_clim;
B_clim = all_clim;
F_clim = all_clim;
P_clim = all_clim;
D_clim = all_clim;
S_clim = all_clim;
M_clim = all_clim;
L_clim = all_clim;
for m = 1:12
    %mo = (m+600):12:nmo; %exclude 1st 50 yrs spinup 
    mo = (m+1200):12:nmo; %exclude 1st 100 yrs
    sf_clim(:,m) = nanmean(sf_biome(:,mo),2);
    sp_clim(:,m) = nanmean(sp_biome(:,mo),2);
    sd_clim(:,m) = nanmean(sd_biome(:,mo),2);
    mf_clim(:,m) = nanmean(mf_biome(:,mo),2);
    mp_clim(:,m) = nanmean(mp_biome(:,mo),2);
    md_clim(:,m) = nanmean(md_biome(:,mo),2);
    lp_clim(:,m) = nanmean(lp_biome(:,mo),2);
    ld_clim(:,m) = nanmean(ld_biome(:,mo),2);
    B_clim(:,m)  = nanmean(b_biome(:,mo),2);
    F_clim(:,m)  = nanmean(AllF_biome(:,mo),2);
    P_clim(:,m)  = nanmean(AllP_biome(:,mo),2);
    D_clim(:,m)  = nanmean(AllD_biome(:,mo),2);
    S_clim(:,m)  = nanmean(AllS_biome(:,mo),2);
    M_clim(:,m)  = nanmean(AllM_biome(:,mo),2);
    L_clim(:,m)  = nanmean(AllL_biome(:,mo),2);
    all_clim(:,m) = nanmean(All_biome(:,mo),2);
end

%% quick plot
figure
subplot(2,2,1)
plot(F_clim(1,:))
subplot(2,2,2)
plot(F_clim(2,:))
subplot(2,2,3)
plot(F_clim(3,:))
subplot(2,2,4)
plot(F_clim(4,:))

figure
subplot(2,2,1)
plot(1:12,F_clim)
subplot(2,2,2)
plot(1:12,P_clim)
subplot(2,2,3)
plot(1:12,D_clim)
subplot(2,2,4)
plot(1:12,B_clim)

%% Calc anomaly ts
all_anom = nan*ones(4,nmo);
sf_anom = all_anom;
sp_anom = all_anom;
sd_anom = all_anom;
mf_anom = all_anom;
mp_anom = all_anom;
md_anom = all_anom;
lp_anom = all_anom;
ld_anom = all_anom;
B_anom = all_anom;
F_anom = all_anom;
P_anom = all_anom;
D_anom = all_anom;
S_anom = all_anom;
M_anom = all_anom;
L_anom = all_anom;

for m = 1:12
    mo = m:12:nmo;
    sf_anom(:,mo) = sf_biome(:,mo) - sf_clim(:,m);
    sp_anom(:,mo) = sp_biome(:,mo) - sp_clim(:,m);
    sd_anom(:,mo) = sd_biome(:,mo) - sd_clim(:,m);
    mf_anom(:,mo) = mf_biome(:,mo) - mf_clim(:,m);
    mp_anom(:,mo) = mp_biome(:,mo) - mp_clim(:,m);
    md_anom(:,mo) = md_biome(:,mo) - md_clim(:,m);
    lp_anom(:,mo) = lp_biome(:,mo) - lp_clim(:,m);
    ld_anom(:,mo) = ld_biome(:,mo) - ld_clim(:,m);
    B_anom(:,mo)  = b_biome(:,mo)  - B_clim(:,m);
    F_anom(:,mo)  = AllF_biome(:,mo) - F_clim(:,m);
    P_anom(:,mo)  = AllP_biome(:,mo) - P_clim(:,m);
    D_anom(:,mo)  = AllD_biome(:,mo) - D_clim(:,m);
    S_anom(:,mo)  = AllS_biome(:,mo) - S_clim(:,m);
    M_anom(:,mo)  = AllM_biome(:,mo) - M_clim(:,m);
    L_anom(:,mo)  = AllL_biome(:,mo) - L_clim(:,m);
    all_anom(:,mo) = All_biome(:,mo) - all_clim(:,m);
end

%% quick plot
% figure
% subplot(2,2,1)
% plot(yid(600:end),F_anom(1,600:end))
% subplot(2,2,2)
% plot(yid(600:end),F_anom(2,600:end))
% subplot(2,2,3)
% plot(yid(600:end),F_anom(3,600:end))
% subplot(2,2,4)
% plot(yid(600:end),F_anom(4,600:end))
% 
% figure
% subplot(2,2,1)
% plot(yid(600:end),F_anom(:,600:end))
% subplot(2,2,2)
% plot(yid(600:end),P_anom(:,600:end))
% subplot(2,2,3)
% plot(yid(600:end),D_anom(:,600:end))
% subplot(2,2,4)
% plot(yid(600:end),B_anom(:,600:end))

figure
subplot(2,2,1)
plot(yid(1200:end),F_anom(1,1200:end))
subplot(2,2,2)
plot(yid(1200:end),F_anom(2,1200:end))
subplot(2,2,3)
plot(yid(1200:end),F_anom(3,1200:end))
subplot(2,2,4)
plot(yid(1200:end),F_anom(4,1200:end))

figure
subplot(2,2,1)
plot(yid(1200:end),F_anom(:,1200:end))
subplot(2,2,2)
plot(yid(1200:end),P_anom(:,1200:end))
subplot(2,2,3)
plot(yid(1200:end),D_anom(:,1200:end))
subplot(2,2,4)
plot(yid(1200:end),B_anom(:,1200:end))

%% save climatologies
save([fpath 'fesity_pi400_biomes_climatol300_ln.mat'],...
    'geolat_t','geolon_t','yid','grid',...
    'sf_clim','sp_clim','sd_clim','mf_clim','mp_clim','md_clim',...
    'lp_clim','ld_clim','B_clim','F_clim','P_clim','D_clim',...
    'S_clim','M_clim','L_clim','all_clim');

%% Save anom
save([fpath 'feisty_pi400_biomes_anom300_ln.mat'],...
    'geolat_t','geolon_t','yid','grid',...
    'sf_anom','sp_anom','sd_anom','mf_anom','mp_anom','md_anom',...
    'lp_anom','ld_anom','B_anom','F_anom','P_anom','D_anom',...
    'S_anom','M_anom','L_anom','all_anom');

