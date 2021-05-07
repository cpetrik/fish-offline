% Calc seasonal climatology - Step 1 of powerspec slope
% Log-transform abund
% Then create anomaly ts - Step 2 of powerspec slope

% Steps
% 1. Calc seasonal climatology
% 2. Create anomaly time series by removing climatol
% 3. Calc Theil-Sen linear trend
% 4. Remove trend if significant
% 5. Calc spectral slopes

clear all
close all

spath='/Volumes/MIP/GCM_DATA/CORE-forced/';
cpath = '/Users/cpetrik/Dropbox/Princeton/POEM_other/cobalt_data/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/FEISTY/NC/Matlab_new_size/' cfile '/'];

%% CORE-forced output

% SP
ncid = netcdf.open([fpath 'Biome_exper_TP_All_fish03_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[nid,nt] = size(biomass);
SPb = biomass;
clear biomass 

% SF
ncid = netcdf.open([fpath 'Biome_exper_TP_All_fish03_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SFb = biomass(:,1:nt);
clear biomass

% SD
ncid = netcdf.open([fpath 'Biome_exper_TP_All_fish03_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SDb = biomass;
clear biomass 

% MP
ncid = netcdf.open([fpath 'Biome_exper_TP_All_fish03_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MPb = biomass;
clear biomass 

% MF
ncid = netcdf.open([fpath 'Biome_exper_TP_All_fish03_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MFb = biomass;
clear biomass

% MD
ncid = netcdf.open([fpath 'Biome_exper_TP_All_fish03_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MDb = biomass;
clear biomass 

% LP
ncid = netcdf.open([fpath 'Biome_exper_TP_All_fish03_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LPb = biomass;
clear biomass 

% LD
ncid = netcdf.open([fpath 'Biome_exper_TP_All_fish03_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LDb = biomass;
clear biomass 

% Benthic material
ncid = netcdf.open([fpath 'Biome_exper_TP_All_fish03_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bb = biomass;
clear biomass 

%%
All = SPb+SFb+SDb+MPb+MFb+MDb+LPb+LDb;
AllF = SFb+MFb;
AllP = SPb+MPb+LPb;
AllD = SDb+MDb+LDb;
AllS = SPb+SFb+SDb;
AllM = MPb+MFb+MDb;
AllL = LPb+LDb;

%% By hand
yid = (time/12);
nmo = length(yid);
nyr = length(yid)/12;

%%
[nid,nt] = size(All);
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
    mo = m:12:nyr;
    sf_clim(:,m) = nanmean(log10(SFb(:,mo)),2);
    sp_clim(:,m) = nanmean(log10(SPb(:,mo)),2);
    sd_clim(:,m) = nanmean(log10(SDb(:,mo)),2);
    mf_clim(:,m) = nanmean(log10(MFb(:,mo)),2);
    mp_clim(:,m) = nanmean(log10(MPb(:,mo)),2);
    md_clim(:,m) = nanmean(log10(MDb(:,mo)),2);
    lp_clim(:,m) = nanmean(log10(LPb(:,mo)),2);
    ld_clim(:,m) = nanmean(log10(LDb(:,mo)),2);
    B_clim(:,m)  = nanmean(log10(Bb(:,mo)),2);
    F_clim(:,m)  = nanmean(log10(AllF(:,mo)),2);
    P_clim(:,m)  = nanmean(log10(AllP(:,mo)),2);
    D_clim(:,m)  = nanmean(log10(AllD(:,mo)),2);
    S_clim(:,m)  = nanmean(log10(AllS(:,mo)),2);
    M_clim(:,m)  = nanmean(log10(AllM(:,mo)),2);
    L_clim(:,m)  = nanmean(log10(AllL(:,mo)),2);
    all_clim(:,m) = nanmean(log10(All(:,mo)),2);
end

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
    sf_anom(:,mo) = log10(SFb(:,mo)) - sf_clim(:,m);
    sp_anom(:,mo) = log10(SPb(:,mo)) - sp_clim(:,m);
    sd_anom(:,mo) = log10(SDb(:,mo)) - sd_clim(:,m);
    mf_anom(:,mo) = log10(MFb(:,mo)) - mf_clim(:,m);
    mp_anom(:,mo) = log10(MPb(:,mo)) - mp_clim(:,m);
    md_anom(:,mo) = log10(MDb(:,mo)) - md_clim(:,m);
    lp_anom(:,mo) = log10(LPb(:,mo)) - lp_clim(:,m);
    ld_anom(:,mo) = log10(LDb(:,mo)) - ld_clim(:,m);
    B_anom(:,mo)  = log10(Bb(:,mo)) - B_clim(:,m);
    F_anom(:,mo)  = log10(AllF(:,mo))  - F_clim(:,m);
    P_anom(:,mo)  = log10(AllP(:,mo))  - P_clim(:,m);
    D_anom(:,mo)  = log10(AllD(:,mo))  - D_clim(:,m);
    S_anom(:,mo)  = log10(AllS(:,mo))  - S_clim(:,m);
    M_anom(:,mo)  = log10(AllM(:,mo))  - M_clim(:,m);
    L_anom(:,mo)  = log10(AllL(:,mo))  - L_clim(:,m);
    all_anom(:,mo) = log10(All(:,mo)) - all_clim(:,m);
 end

%% Save anom and climatologies
save([fpath 'fesity_biome_tp_exper_climatol_anom.mat'],...
    'sf_clim','sp_clim','sd_clim','mf_clim','mp_clim','md_clim',...
    'lp_clim','ld_clim','B_clim','F_clim','P_clim','D_clim',...
    'S_clim','M_clim','L_clim','all_clim',...
    'sf_anom','sp_anom','sd_anom','mf_anom','mp_anom','md_anom',...
    'lp_anom','ld_anom','B_anom','F_anom','P_anom','D_anom',...
    'S_anom','M_anom','L_anom','all_anom');

