% Calc seasonal climatology - Step 1 of powerspec slope
% ln-transform abund
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

[ni,nj]=size(geolon_t);

%% Preindust 400 yr output
load([fpath 'Monthly_Pre400_' cfile '.mat']);

% % SP
% ncid = netcdf.open([fpath 'Core_All_fish03_sml_p.nc'],'NC_NOWRITE');
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% for i = 1:nvars
%     varname = netcdf.inqVar(ncid, i-1);
%     eval([ varname ' = netcdf.getVar(ncid,i-1);']);
%     eval([ varname '(' varname ' == 99999) = NaN;']);
% end
% netcdf.close(ncid);
% 
% [nid,nt] = size(biomass);
% SP.bio = biomass;
% clear biomass 
% 
% % SF
% ncid = netcdf.open([fpath 'Core_All_fish03_sml_f.nc'],'NC_NOWRITE');
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% for i = 1:nvars
%     varname = netcdf.inqVar(ncid, i-1);
%     eval([ varname ' = netcdf.getVar(ncid,i-1);']);
%     eval([ varname '(' varname ' == 99999) = NaN;']);
% end
% netcdf.close(ncid);
% 
% SF.bio = biomass(:,1:nt);
% clear biomass
% 
% % SD
% ncid = netcdf.open([fpath 'Core_All_fish03_sml_d.nc'],'NC_NOWRITE');
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% for i = 1:nvars
%     varname = netcdf.inqVar(ncid, i-1);
%     eval([ varname ' = netcdf.getVar(ncid,i-1);']);
%     eval([ varname '(' varname ' == 99999) = NaN;']);
% end
% netcdf.close(ncid);
% 
% SD.bio = biomass;
% clear biomass 
% 
% % MP
% ncid = netcdf.open([fpath 'Core_All_fish03_med_p.nc'],'NC_NOWRITE');
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% for i = 1:nvars
%     varname = netcdf.inqVar(ncid, i-1);
%     eval([ varname ' = netcdf.getVar(ncid,i-1);']);
%     eval([ varname '(' varname ' == 99999) = NaN;']);
% end
% netcdf.close(ncid);
% 
% MP.bio = biomass;
% clear biomass 
% 
% % MF
% ncid = netcdf.open([fpath 'Core_All_fish03_med_f.nc'],'NC_NOWRITE');
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% for i = 1:nvars
%     varname = netcdf.inqVar(ncid, i-1);
%     eval([ varname ' = netcdf.getVar(ncid,i-1);']);
%     eval([ varname '(' varname ' == 99999) = NaN;']);
% end
% netcdf.close(ncid);
% 
% MF.bio = biomass;
% clear biomass
% 
% % MD
% ncid = netcdf.open([fpath 'Core_All_fish03_med_d.nc'],'NC_NOWRITE');
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% for i = 1:nvars
%     varname = netcdf.inqVar(ncid, i-1);
%     eval([ varname ' = netcdf.getVar(ncid,i-1);']);
%     eval([ varname '(' varname ' == 99999) = NaN;']);
% end
% netcdf.close(ncid);
% 
% MD.bio = biomass;
% clear biomass 
% 
% % LP
% ncid = netcdf.open([fpath 'Core_All_fish03_lrg_p.nc'],'NC_NOWRITE');
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% for i = 1:nvars
%     varname = netcdf.inqVar(ncid, i-1);
%     eval([ varname ' = netcdf.getVar(ncid,i-1);']);
%     eval([ varname '(' varname ' == 99999) = NaN;']);
% end
% netcdf.close(ncid);
% 
% LP.bio = biomass;
% clear biomass 
% 
% % LD
% ncid = netcdf.open([fpath 'Core_All_fish03_lrg_d.nc'],'NC_NOWRITE');
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% for i = 1:nvars
%     varname = netcdf.inqVar(ncid, i-1);
%     eval([ varname ' = netcdf.getVar(ncid,i-1);']);
%     eval([ varname '(' varname ' == 99999) = NaN;']);
% end
% netcdf.close(ncid);
% 
% LD.bio = biomass;
% clear biomass 
% 
% % Benthic material
% ncid = netcdf.open([fpath 'Core_All_fish03_bent.nc'],'NC_NOWRITE');
% [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
% for i = 1:nvars
%     varname = netcdf.inqVar(ncid, i-1);
%     eval([ varname ' = netcdf.getVar(ncid,i-1);']);
%     eval([ varname '(' varname ' == 99999) = NaN;']);
% end
% netcdf.close(ncid);
% 
% Bent.bio = biomass;
% clear biomass 


%%
% All = SP.bio+SF.bio+SD.bio+MP.bio+MF.bio+MD.bio+LP.bio+LD.bio;
% AllF = SF.bio+MF.bio;
% AllP = SP.bio+MP.bio+LP.bio;
% AllD = SD.bio+MD.bio+LD.bio;
% AllS = SP.bio+SF.bio+SD.bio;
% AllM = MP.bio+MF.bio+MD.bio;
% AllL = LP.bio+LD.bio;

All = sp_bio+sf_bio+sd_bio+mp_bio+mf_bio+md_bio+lp_bio+ld_bio;
AllF = sf_bio+mf_bio;
AllP = sp_bio+mp_bio+lp_bio;
AllD = sd_bio+md_bio+ld_bio;
AllS = sp_bio+sf_bio+sd_bio;
AllM = mp_bio+mf_bio+md_bio;
AllL = lp_bio+ld_bio;

%% By hand
yid = (time/12)+1000;
nmo = length(yid);
nyr = length(yid)/12;

%%
[nid,nt] = size(sf_bio);
all_clim = nan*ones(nid,12);
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
    mo = (m+600):12:nmo; %exclude 1st 50 yrs spinup
%     sf_clim(:,m) = nanmean(log(SF.bio(:,mo)),2);
%     sp_clim(:,m) = nanmean(log(SP.bio(:,mo)),2);
%     sd_clim(:,m) = nanmean(log(SD.bio(:,mo)),2);
%     mf_clim(:,m) = nanmean(log(MF.bio(:,mo)),2);
%     mp_clim(:,m) = nanmean(log(MP.bio(:,mo)),2);
%     md_clim(:,m) = nanmean(log(MD.bio(:,mo)),2);
%     lp_clim(:,m) = nanmean(log(LP.bio(:,mo)),2);
%     ld_clim(:,m) = nanmean(log(LD.bio(:,mo)),2);
%     B_clim(:,m)  = nanmean(log(Bent.bio(:,mo)),2);
    sf_clim(:,m) = nanmean(log(sf_bio(:,mo)),2);
    sp_clim(:,m) = nanmean(log(sp_bio(:,mo)),2);
    sd_clim(:,m) = nanmean(log(sd_bio(:,mo)),2);
    mf_clim(:,m) = nanmean(log(mf_bio(:,mo)),2);
    mp_clim(:,m) = nanmean(log(mp_bio(:,mo)),2);
    md_clim(:,m) = nanmean(log(md_bio(:,mo)),2);
    lp_clim(:,m) = nanmean(log(lp_bio(:,mo)),2);
    ld_clim(:,m) = nanmean(log(ld_bio(:,mo)),2);
    B_clim(:,m)  = nanmean(log(b_bio(:,mo)),2);
    F_clim(:,m)  = nanmean(log(AllF(:,mo)),2);
    P_clim(:,m)  = nanmean(log(AllP(:,mo)),2);
    D_clim(:,m)  = nanmean(log(AllD(:,mo)),2);
    S_clim(:,m)  = nanmean(log(AllS(:,mo)),2);
    M_clim(:,m)  = nanmean(log(AllM(:,mo)),2);
    L_clim(:,m)  = nanmean(log(AllL(:,mo)),2);
    all_clim(:,m) = nanmean(log(All(:,mo)),2);
end

%% quick plot
figure
subplot(2,2,1)
plot(F_clim(1,:))
subplot(2,2,2)
plot(F_clim(1000,:))
subplot(2,2,3)
plot(F_clim(10000,:))
subplot(2,2,4)
plot(F_clim(48000,:))

figure
subplot(2,2,1)
plot(F_clim(10,:))
subplot(2,2,2)
plot(P_clim(10,:))
subplot(2,2,3)
plot(D_clim(10,:))
subplot(2,2,4)
plot(B_clim(10,:))

%% Calc anomaly ts
all_anom = nan*ones(nid,nmo);
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
    sf_anom(:,mo) = log(sf_bio(:,mo)) - sf_clim(:,m);
    sp_anom(:,mo) = log(sp_bio(:,mo)) - sp_clim(:,m);
    sd_anom(:,mo) = log(sd_bio(:,mo)) - sd_clim(:,m);
    mf_anom(:,mo) = log(mf_bio(:,mo)) - mf_clim(:,m);
    mp_anom(:,mo) = log(mp_bio(:,mo)) - mp_clim(:,m);
    md_anom(:,mo) = log(md_bio(:,mo)) - md_clim(:,m);
    lp_anom(:,mo) = log(lp_bio(:,mo)) - lp_clim(:,m);
    ld_anom(:,mo) = log(ld_bio(:,mo)) - ld_clim(:,m);
    B_anom(:,mo)  = log(b_bio(:,mo))  - B_clim(:,m);
    F_anom(:,mo)  = log(AllF(:,mo))   - F_clim(:,m);
    P_anom(:,mo)  = log(AllP(:,mo))   - P_clim(:,m);
    D_anom(:,mo)  = log(AllD(:,mo))   - D_clim(:,m);
    S_anom(:,mo)  = log(AllS(:,mo))   - S_clim(:,m);
    M_anom(:,mo)  = log(AllM(:,mo))   - M_clim(:,m);
    L_anom(:,mo)  = log(AllL(:,mo))   - L_clim(:,m);
    all_anom(:,mo) = log(All(:,mo))   - all_clim(:,m);
end

%% quick plot
figure
subplot(2,2,1)
plot(yid(600:end),F_anom(1,600:end))
subplot(2,2,2)
plot(yid(600:end),F_anom(1000,600:end))
subplot(2,2,3)
plot(yid(600:end),F_anom(10000,600:end))
subplot(2,2,4)
plot(yid(600:end),F_anom(48000,600:end))

figure
subplot(2,2,1)
plot(yid(600:end),F_anom(10,600:end))
subplot(2,2,2)
plot(yid(600:end),P_anom(10,600:end))
subplot(2,2,3)
plot(yid(600:end),D_anom(10,600:end))
subplot(2,2,4)
plot(yid(600:end),B_anom(10,600:end))

%% save climatologies
save([fpath 'fesity_pi400_climatol_ln.mat'],...
    'geolat_t','geolon_t','yid','grid',...
    'sf_clim','sp_clim','sd_clim','mf_clim','mp_clim','md_clim',...
    'lp_clim','ld_clim','B_clim','F_clim','P_clim','D_clim',...
    'S_clim','M_clim','L_clim','all_clim');

%% Save anom
save([fpath 'feisty_pi400_anom_ln.mat'],...
    'geolat_t','geolon_t','yid','grid',...
    'sf_anom','sp_anom','sd_anom','mf_anom','mp_anom','md_anom',...
    'lp_anom','ld_anom','B_anom','F_anom','P_anom','D_anom',...
    'S_anom','M_anom','L_anom','all_anom');

