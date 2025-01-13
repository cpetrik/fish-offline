% FEISTY output at all locations

clear 
close all

%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
harv = 'v15_All_fish03';

%% MP
ncid = netcdf.open([spath 'FOSI_' harv '_catch_med_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([spath 'FOSI_' harv '_catch_med_f.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([spath 'FOSI_' harv '_catch_med_d.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([spath 'FOSI_' harv '_catch_lrg_p.nc'],'NC_NOWRITE');
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
ncid = netcdf.open([spath 'FOSI_' harv '_catch_lrg_d.nc'],'NC_NOWRITE');
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
% Catch totals only in lmes
fpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([fpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([fpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([fpath 'LME-mask-POP_gx1v6.mat']);
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
mf_tmy=mean(MF.yield,1,'omitnan');
mp_tmy=mean(MP.yield,1,'omitnan');
md_tmy=mean(MD.yield,1,'omitnan');
lp_tmy=mean(LP.yield,1,'omitnan');
ld_tmy=mean(LD.yield,1,'omitnan');
%mean catch per mo
mf_tmc=mean(MF.catch,1,'omitnan');
mp_tmc=mean(MP.catch,1,'omitnan');
md_tmc=mean(MD.catch,1,'omitnan');
lp_tmc=mean(LP.catch,1,'omitnan');
ld_tmc=mean(LD.catch,1,'omitnan');

%total yield per mo
mf_tty=sum(MF.yield,1,'omitnan');
mp_tty=sum(MP.yield,1,'omitnan');
md_tty=sum(MD.yield,1,'omitnan');
lp_tty=sum(LP.yield,1,'omitnan');
ld_tty=sum(LD.yield,1,'omitnan');
%total catch per mo
mf_ttc=sum(MF.catch,1,'omitnan');
mp_ttc=sum(MP.catch,1,'omitnan');
md_ttc=sum(MD.catch,1,'omitnan');
lp_ttc=sum(LP.catch,1,'omitnan');
ld_ttc=sum(LD.catch,1,'omitnan');


%% Spatially
%mean yield per mo
mf_smy=mean(MF.yield,2,'omitnan');
mp_smy=mean(MP.yield,2,'omitnan');
md_smy=mean(MD.yield,2,'omitnan');
lp_smy=mean(LP.yield,2,'omitnan');
ld_smy=mean(LD.yield,2,'omitnan');
%mean catch per mo
mf_smc=mean(MF.catch,2,'omitnan');
mp_smc=mean(MP.catch,2,'omitnan');
md_smc=mean(MD.catch,2,'omitnan');
lp_smc=mean(LP.catch,2,'omitnan');
ld_smc=mean(LD.catch,2,'omitnan');

%total yield per mo
mf_sty=sum(MF.yield,2,'omitnan');
mp_sty=sum(MP.yield,2,'omitnan');
md_sty=sum(MD.yield,2,'omitnan');
lp_sty=sum(LP.yield,2,'omitnan');
ld_sty=sum(LD.yield,2,'omitnan');
%total catch per mo
mf_stc=sum(MF.catch,2,'omitnan');
mp_stc=sum(MP.catch,2,'omitnan');
md_stc=sum(MD.catch,2,'omitnan');
lp_stc=sum(LP.catch,2,'omitnan');
ld_stc=sum(LD.catch,2,'omitnan');


%% Every year
st=1:12:nt;
en=12:12:nt;

mf_tac = nan*ones(ni,nyr);
mp_tac = nan*ones(ni,nyr);
md_tac = nan*ones(ni,nyr);
lp_tac = nan*ones(ni,nyr);
ld_tac = nan*ones(ni,nyr);
mf_may = nan*ones(ni,nyr);
mp_may = nan*ones(ni,nyr);
md_may = nan*ones(ni,nyr);
lp_may = nan*ones(ni,nyr);
ld_may = nan*ones(ni,nyr);

for n=1:length(st)

    mp_tac(:,n)=sum(MP.catch(:,st(n):en(n)),2,'omitnan');
    mf_tac(:,n)=sum(MF.catch(:,st(n):en(n)),2,'omitnan');
    md_tac(:,n)=sum(MD.catch(:,st(n):en(n)),2,'omitnan');
    lp_tac(:,n)=sum(LP.catch(:,st(n):en(n)),2,'omitnan');
    ld_tac(:,n)=sum(LD.catch(:,st(n):en(n)),2,'omitnan');

    mp_may(:,n)=mean(MP.yield(:,st(n):en(n)),2,'omitnan');
    mf_may(:,n)=mean(MF.yield(:,st(n):en(n)),2,'omitnan');
    md_may(:,n)=mean(MD.yield(:,st(n):en(n)),2,'omitnan');
    lp_may(:,n)=mean(LP.yield(:,st(n):en(n)),2,'omitnan');
    ld_may(:,n)=mean(LD.yield(:,st(n):en(n)),2,'omitnan');
end

tmn = mf_tac + mp_tac + md_tac + lp_tac + ld_tac;
stmn = sum(tmn,'omitnan');

mp_tsac = sum(mp_tac,'omitnan');
mf_tsac = sum(mf_tac,'omitnan');
md_tsac = sum(md_tac,'omitnan');
lp_tsac = sum(lp_tac,'omitnan');
ld_tsac = sum(ld_tac,'omitnan');

%%
save([spath 'Time_Means_FOSI_' harv '_' cfile '.mat'],'units_yield','units_catch',...
    'mf_tmy','mp_tmy','md_tmy','lp_tmy','ld_tmy',...
    'mf_tmc','mp_tmc','md_tmc','lp_tmc','ld_tmc',...
    'mf_tty','mp_tty','md_tty','lp_tty','ld_tty',...
    'mf_ttc','mp_ttc','md_ttc','lp_ttc','ld_ttc',...
    '-append');
    %'mf_tnu','lp_tnu','ld_tnu','-append');

save([spath 'Space_Means_FOSI_' harv '_' cfile '.mat'],'units_yield','units_catch',...
    'mf_smy','mp_smy','md_smy','lp_smy','ld_smy',...
    'mf_smc','mp_smc','md_smc','lp_smc','ld_smc',...
    'mf_sty','mp_sty','md_sty','lp_sty','ld_sty',...
    'mf_stc','mp_stc','md_stc','lp_stc','ld_stc',...
    '-append');
    %'mf_snu','lp_snu','ld_snu','-append');

save([spath 'Annual_Means_FOSI_' harv '_' cfile '.mat'],...
    'mf_tac','mp_tac','md_tac','lp_tac','ld_tac',...
    'mf_tsac','mp_tsac','md_tsac','lp_tsac','ld_tsac',...
    'units_yield','units_catch',...
    'mf_may','mp_may','md_may','lp_may','ld_may','-append');
    %'mf_anu','lp_anu','ld_anu','-append');

%%
% mo = (1:nt)/12;
% figure
% plot(mo,(lp_tnu),'b'); hold on;
% plot(mo,(mf_tnu),'r'); hold on;
% plot(mo,(ld_tnu),'k'); hold on;
% 
% figure
% plot(mo,5.5902e3*(lp_tnu),'b'); hold on;
% plot(mo,11.1803*(mf_tnu),'r'); hold on;
% plot(mo,5.5902e3*(ld_tnu),'k'); hold on;
% 

