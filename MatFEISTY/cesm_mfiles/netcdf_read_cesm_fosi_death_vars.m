% FEISTY output at all locations

clear 
close all

%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
harv = 'v15_All_fish03_';

%% SP
ncid = netcdf.open([fpath 'FOSI_' harv 'die_nmort_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
[ni,nt] = size(die);

SP.die = die;
SP.mort = mort;
clear die mort

%% SF
ncid = netcdf.open([fpath 'FOSI_' harv 'die_nmort_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.die = die;
SF.mort = mort;
clear die mort

% SD
ncid = netcdf.open([fpath 'FOSI_' harv 'die_nmort_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.die = die;
SD.mort = mort;
clear die mort

% MP
ncid = netcdf.open([fpath 'FOSI_' harv 'die_nmort_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.die = die;
MP.mort = mort;
clear die mort

% MD
ncid = netcdf.open([fpath 'FOSI_' harv 'die_nmort_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.die = die;
MD.mort = mort;
clear die mort

%% MF
ncid = netcdf.open([fpath 'FOSI_' harv 'die_nmort_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(die);

MF.die = die;
MF.mort = mort;
clear die mort

% LP
ncid = netcdf.open([fpath 'FOSI_' harv 'die_nmort_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.die = die;
LP.mort = mort;
clear die mort

% LD
ncid = netcdf.open([fpath 'FOSI_' harv 'die_nmort_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.die = die;
LD.mort = mort;
clear die mort

%% Take means for visualization

%Time
sf_tdie = mean(SF.die,1,"omitnan");
sp_tdie = mean(SP.die,1,"omitnan");
sd_tdie = mean(SD.die,1,"omitnan");
mf_tdie = mean(MF.die,1,"omitnan");
mp_tdie = mean(MP.die,1,"omitnan");
md_tdie = mean(MD.die,1,"omitnan");
lp_tdie = mean(LP.die,1,"omitnan");
ld_tdie = mean(LD.die,1,"omitnan");

sf_tmort = mean(SF.mort,1,"omitnan");
sp_tmort = mean(SP.mort,1,"omitnan");
sd_tmort = mean(SD.mort,1,"omitnan");
mf_tmort = mean(MF.mort,1,"omitnan");
mp_tmort = mean(MP.mort,1,"omitnan");
md_tmort = mean(MD.mort,1,"omitnan");
lp_tmort = mean(LP.mort,1,"omitnan");
ld_tmort = mean(LD.mort,1,"omitnan");

%% Space
sf_sdie = mean(SF.die,2,"omitnan");
sp_sdie = mean(SP.die,2,"omitnan");
sd_sdie = mean(SD.die,2,"omitnan");
mf_sdie = mean(MF.die,2,"omitnan");
mp_sdie = mean(MP.die,2,"omitnan");
md_sdie = mean(MD.die,2,"omitnan");
lp_sdie = mean(LP.die,2,"omitnan");
ld_sdie = mean(LD.die,2,"omitnan");

sf_smort = mean(SF.mort,2,"omitnan");
sp_smort = mean(SP.mort,2,"omitnan");
sd_smort = mean(SD.mort,2,"omitnan");
mf_smort = mean(MF.mort,2,"omitnan");
mp_smort = mean(MP.mort,2,"omitnan");
md_smort = mean(MD.mort,2,"omitnan");
lp_smort = mean(LP.mort,2,"omitnan");
ld_smort = mean(LD.mort,2,"omitnan");

%% Annual means
nyr = nt/12;
st=1:12:nt;
en=12:12:nt;

for n=1:length(st)
    % mean die
    sf_adie(:,n)=mean(SF.die(:,st(n):en(n)),2,"omitnan");
    sp_adie(:,n)=mean(SP.die(:,st(n):en(n)),2,"omitnan");
    sd_adie(:,n)=mean(SD.die(:,st(n):en(n)),2,"omitnan");
    mf_adie(:,n)=mean(MF.die(:,st(n):en(n)),2,"omitnan");
    mp_adie(:,n)=mean(MP.die(:,st(n):en(n)),2,"omitnan");
    md_adie(:,n)=mean(MD.die(:,st(n):en(n)),2,"omitnan");
    lp_adie(:,n)=mean(LP.die(:,st(n):en(n)),2,"omitnan");
    ld_adie(:,n)=mean(LD.die(:,st(n):en(n)),2,"omitnan");
    
    % mean mort
    sf_amort(:,n)=mean(SF.mort(:,st(n):en(n)),2,"omitnan");
    sp_amort(:,n)=mean(SP.mort(:,st(n):en(n)),2,"omitnan");
    sd_amort(:,n)=mean(SD.mort(:,st(n):en(n)),2,"omitnan");
    mf_amort(:,n)=mean(MF.mort(:,st(n):en(n)),2,"omitnan");
    mp_amort(:,n)=mean(MP.mort(:,st(n):en(n)),2,"omitnan");
    md_amort(:,n)=mean(MD.mort(:,st(n):en(n)),2,"omitnan");
    lp_amort(:,n)=mean(LP.mort(:,st(n):en(n)),2,"omitnan");
    ld_amort(:,n)=mean(LD.mort(:,st(n):en(n)),2,"omitnan");

end

%%
save([fpath 'Time_Means_FOSI_' harv cfile '.mat'],...
    'sf_tdie','sp_tdie','sd_tdie',...
    'mf_tdie','mp_tdie','md_tdie',...
    'lp_tdie','ld_tdie',...
    'sf_tmort','sp_tmort','sd_tmort',...
    'mf_tmort','mp_tmort','md_tmort',...
    'lp_tmort','ld_tmort','-append')

save([fpath 'Space_Means_FOSI_' harv cfile '.mat'],...
    'sf_sdie','sp_sdie','sd_sdie',...
    'mf_sdie','mp_sdie','md_sdie',...
    'lp_sdie','ld_sdie',...
    'sf_smort','sp_smort','sd_smort',...
    'mf_smort','mp_smort','md_smort',...
    'lp_smort','ld_smort','-append')

save([fpath 'Annual_Means_FOSI_' harv cfile '.mat'],...
    'sf_adie','sp_adie','sd_adie',...
    'mf_adie','mp_adie','md_adie',...
    'lp_adie','ld_adie',...
    'sf_amort','sp_amort','sd_amort',...
    'mf_amort','mp_amort','md_amort',...
    'lp_amort','ld_amort','-append')

%%
mo = (1:nt)/12;
figure
plot(mo,(mp_tdie),'b'); hold on;
plot(mo,(mf_tdie),'r'); hold on;
plot(mo,(md_tdie),'k'); hold on;
title('die')

figure
plot(mo,(mp_tmort),'b'); hold on;
plot(mo,(mf_tmort),'r'); hold on;
plot(mo,(md_tmort),'k'); hold on;
title('mort')


