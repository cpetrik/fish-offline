% Calc temporal variance of FEISTY with
% CORE-forced Cobalt inputs

clear all
close all

spath='/Volumes/MIP/GCM_DATA/CORE-forced/';

cfile = 'Dc_enc70-b200_m4-b175-k086_c20-b250_D075_J100_A050_Sm025_nmort1_BE08_noCC_RE00100';
fpath=['/Volumes/MIP/NC/Matlab_new_size/' cfile '/'];

load('/Volumes/MIP/GCM_DATA/CORE-forced/ocean_cobalt_grid.mat',...
    'geolon_t','geolat_t');
load('/Volumes/MIP/GCM_DATA/CORE-forced/Data_grid_ocean_cobalt_ESM2Mcore.mat',...
    'GRD');

[ni,nj]=size(geolon_t);

%% CORE-forced output

% SP
ncid = netcdf.open([fpath 'Core_All_fish03_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[nid,nt] = size(biomass);
SP.bio = biomass;
clear biomass 

% SF
ncid = netcdf.open([fpath 'Core_All_fish03_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.bio = biomass(:,1:nt);
clear biomass

% SD
ncid = netcdf.open([fpath 'Core_All_fish03_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.bio = biomass;
clear biomass 

% MP
ncid = netcdf.open([fpath 'Core_All_fish03_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.bio = biomass;
clear biomass 

% MF
ncid = netcdf.open([fpath 'Core_All_fish03_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.bio = biomass;
clear biomass

% MD
ncid = netcdf.open([fpath 'Core_All_fish03_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.bio = biomass;
clear biomass 

% LP
ncid = netcdf.open([fpath 'Core_All_fish03_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.bio = biomass;
clear biomass 

% LD
ncid = netcdf.open([fpath 'Core_All_fish03_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.bio = biomass;
clear biomass 

% Benthic material
ncid = netcdf.open([fpath 'Core_All_fish03_bent.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

Bent.bio = biomass;
clear biomass 


%%
All = SP.bio+SF.bio+SD.bio+MP.bio+MF.bio+MD.bio+LP.bio+LD.bio;
AllF = SF.bio+MF.bio;
AllP = SP.bio+MP.bio+LP.bio;
AllD = SD.bio+MD.bio+LD.bio;
AllS = SP.bio+SF.bio+SD.bio;
AllM = MP.bio+MF.bio+MD.bio;
AllL = LP.bio+LD.bio;

%% %% variance
sf_var = var(SF.bio,0,2);
sp_var = var(SP.bio,0,2);
sd_var = var(SD.bio,0,2);
mf_var = var(MF.bio,0,2);
mp_var = var(MP.bio,0,2);
md_var = var(MD.bio,0,2);
lp_var = var(LP.bio,0,2);
ld_var = var(LD.bio,0,2);
F_var = var(AllF,0,2);
P_var = var(AllP,0,2);
D_var = var(AllD,0,2);
S_var = var(AllS,0,2);
M_var = var(AllM,0,2);
L_var = var(AllL,0,2);
All_var = var(All,0,2);
B_var = var(Bent.bio,0,2);

var_vec(:,1) = sf_var;
var_vec(:,2) = sp_var;
var_vec(:,3) = sd_var;
var_vec(:,4) = mf_var;
var_vec(:,5) = mp_var;
var_vec(:,6) = md_var;
var_vec(:,7) = lp_var;
var_vec(:,8) = ld_var;
var_vec(:,9) = F_var;
var_vec(:,10) = P_var;
var_vec(:,11) = D_var;
var_vec(:,12) = S_var;
var_vec(:,13) = M_var;
var_vec(:,14) = L_var;
var_vec(:,15) = All_var;
var_vec(:,16) = B_var;

var_vec_tex = {'SF','SP','SD','MF','MP','MD','LP','LD',...
    'F','P','D','S','M','L','All','B'};

%%
save([fpath 'fesity_core_variance_1950_2007.mat'],...
    'geolat_t','geolon_t','GRD',...
    'sf_var','sp_var','sd_var','mf_var','mp_var','md_var',...
    'lp_var','ld_var','B_var','F_var','P_var','D_var',...
    'S_var','M_var','L_var','All_var',...
    'var_vec','var_vec_tex')

