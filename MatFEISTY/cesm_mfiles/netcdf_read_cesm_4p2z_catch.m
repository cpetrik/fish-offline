% FEISTY output at all locations

clear 
close all

%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b200-k086_c20-b250_D075_A065_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100/4P2Z

gpath='/Volumes/petrik-lab/Feisty/GCM_Data/CESM/4P2Z/';
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/4P2Z/'];
harv = 'All_fish03';

%% MP
ncid = netcdf.open([fpath '4P2Z_' harv '_catch_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);
%%
[ni,nt] = size(yield);

MP.yield = yield;
clear yield

% MF
ncid = netcdf.open([fpath '4P2Z_' harv '_catch_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MF.yield = yield;
clear yield

% MD
ncid = netcdf.open([fpath '4P2Z_' harv '_catch_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.yield = yield;
clear yield

% LP
ncid = netcdf.open([fpath '4P2Z_' harv '_catch_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.yield = yield;
clear yield

% LD
ncid = netcdf.open([fpath '4P2Z_' harv '_catch_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.yield = yield;
clear yield

%% Take means and totals
% Totals only in lmes
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);
ID = GRD.ID;

tlme = double(lme_mask);
tlme(tlme<0) = nan;
tlme(~isnan(tlme)) = 1;
lme_grid = tlme(ID);

%TAREA units 'cm^2'
AREA_OCN = TAREA * 1e-4;
area = AREA_OCN(ID);
area_km2 = area * 1e-6;

MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
nyr = nt/12;
mos = repmat(MNTH,ni,nyr);
mns = repmat(MNTH,ni,1);
area_mat = repmat(area_km2,1,nt);
lme_mat = repmat(lme_grid,1,nt);

%% Units
units_yield = 'g_m2_day';
units_catch = 'g_km2_mo';

MF.catch = MF.yield .*mos .*area_mat .*lme_mat;
MP.catch = MP.yield .*mos .*area_mat .*lme_mat;
MD.catch = MD.yield .*mos .*area_mat .*lme_mat;
LP.catch = LP.yield .*mos .*area_mat .*lme_mat;
LD.catch = LD.yield .*mos .*area_mat .*lme_mat;

%% Time
%mean yield per mo
mf_tmy=nanmean(MF.yield,1);
mp_tmy=nanmean(MP.yield,1);
md_tmy=nanmean(MD.yield,1);
lp_tmy=nanmean(LP.yield,1);
ld_tmy=nanmean(LD.yield,1);
%mean catch per mo
mf_tmc=nanmean(MF.catch,1);
mp_tmc=nanmean(MP.catch,1);
md_tmc=nanmean(MD.catch,1);
lp_tmc=nanmean(LP.catch,1);
ld_tmc=nanmean(LD.catch,1);

%total yield per mo
mf_tty=nansum(MF.yield,1);
mp_tty=nansum(MP.yield,1);
md_tty=nansum(MD.yield,1);
lp_tty=nansum(LP.yield,1);
ld_tty=nansum(LD.yield,1);
%total catch per mo
mf_ttc=nansum(MF.catch,1);
mp_ttc=nansum(MP.catch,1);
md_ttc=nansum(MD.catch,1);
lp_ttc=nansum(LP.catch,1);
ld_ttc=nansum(LD.catch,1);

%% Spatially
%mean yield per mo
mf_smy=nanmean(MF.yield,2);
mp_smy=nanmean(MP.yield,2);
md_smy=nanmean(MD.yield,2);
lp_smy=nanmean(LP.yield,2);
ld_smy=nanmean(LD.yield,2);
%mean catch per mo
mf_smc=nanmean(MF.catch,2);
mp_smc=nanmean(MP.catch,2);
md_smc=nanmean(MD.catch,2);
lp_smc=nanmean(LP.catch,2);
ld_smc=nanmean(LD.catch,2);

%total yield per mo
mf_sty=nansum(MF.yield,2);
mp_sty=nansum(MP.yield,2);
md_sty=nansum(MD.yield,2);
lp_sty=nansum(LP.yield,2);
ld_sty=nansum(LD.yield,2);
%total catch per mo
mf_stc=nansum(MF.catch,2);
mp_stc=nansum(MP.catch,2);
md_stc=nansum(MD.catch,2);
lp_stc=nansum(LP.catch,2);
ld_stc=nansum(LD.catch,2);

%% Every year
st=1:12:nt;
en=12:12:nt;

mf_tac = nan*ones(ni,nyr);
mp_tac = nan*ones(ni,nyr);
md_tac = nan*ones(ni,nyr);
lp_tac = nan*ones(ni,nyr);
ld_tac = nan*ones(ni,nyr);
for n=1:length(st)

    mp_tac(:,n)=nansum(MP.catch(:,st(n):en(n)),2);
    mf_tac(:,n)=nansum(MF.catch(:,st(n):en(n)),2);
    md_tac(:,n)=nansum(MD.catch(:,st(n):en(n)),2);
    lp_tac(:,n)=nansum(LP.catch(:,st(n):en(n)),2);
    ld_tac(:,n)=nansum(LD.catch(:,st(n):en(n)),2);
end

tmn = mf_tac + mp_tac + md_tac + lp_tac + ld_tac;
stmn = sum(tmn);

mp_tsac = nansum(mp_tac);
mf_tsac = nansum(mf_tac);
md_tsac = nansum(md_tac);
lp_tsac = nansum(lp_tac);
ld_tsac = nansum(ld_tac);

%%
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/4P2Z/'];
save([fpath 'Time_Means_4P2Z_1deg_' cfile '.mat'],'units_yield','units_catch',...
    'mf_tmy','mp_tmy','md_tmy','lp_tmy','ld_tmy',...
    'mf_tmc','mp_tmc','md_tmc','lp_tmc','ld_tmc',...
    'mf_tty','mp_tty','md_tty','lp_tty','ld_tty',...
    'mf_ttc','mp_ttc','md_ttc','lp_ttc','ld_ttc','-append');

save([fpath 'Space_Means_4P2Z_1deg_' cfile '.mat'],'units_yield','units_catch',...
    'mf_smy','mp_smy','md_smy','lp_smy','ld_smy',...
    'mf_smc','mp_smc','md_smc','lp_smc','ld_smc',...
    'mf_sty','mp_sty','md_sty','lp_sty','ld_sty',...
    'mf_stc','mp_stc','md_stc','lp_stc','ld_stc','-append');

save([fpath 'Annual_Means_4P2Z_1deg_' cfile '.mat'],...
    'mf_tac','mp_tac','md_tac','lp_tac','ld_tac',...
    'mf_tsac','mp_tsac','md_tsac','lp_tsac','ld_tsac',...
    'units_yield','units_catch','-append');

%%
