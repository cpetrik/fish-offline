% FEISTY output at all locations

clear 
close all

%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
harv = 'v15_All_fish03_';

%% SF
ncid = netcdf.open([fpath 'FOSI_' harv 'nu_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(nu);

SF.nu = nu;
clear nu gamma

% SP
ncid = netcdf.open([fpath 'FOSI_' harv 'nu_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.nu = nu;
clear nu gamma

% SD
ncid = netcdf.open([fpath 'FOSI_' harv 'nu_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.nu = nu;
clear nu gamma

%% MF
ncid = netcdf.open([fpath 'FOSI_' harv 'nu_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(nu);

MF.nu = nu;
clear nu gamma

% MP
ncid = netcdf.open([fpath 'FOSI_' harv 'nu_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.nu = nu;
clear nu gamma

% MD
ncid = netcdf.open([fpath 'FOSI_' harv 'nu_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.nu = nu;
clear nu gamma


%% LP
ncid = netcdf.open([fpath 'FOSI_' harv 'nu_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.nu = nu;
clear nu gamma

% LD
ncid = netcdf.open([fpath 'FOSI_' harv 'nu_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.nu = nu;
clear nu gamma

%% Take means for visualization

%Time
sf_tnu = mean(SF.nu,1,"omitnan");
sp_tnu = mean(SP.nu,1,"omitnan");
sd_tnu = mean(SD.nu,1,"omitnan");
mf_tnu = mean(MF.nu,1,"omitnan");
mp_tnu = mean(MP.nu,1,"omitnan");
md_tnu = mean(MD.nu,1,"omitnan");
lp_tnu = mean(LP.nu,1,"omitnan");
ld_tnu = mean(LD.nu,1,"omitnan");

%% Space
sf_snu = mean(SF.nu,2,"omitnan");
sp_snu = mean(SP.nu,2,"omitnan");
sd_snu = mean(SD.nu,2,"omitnan");
mf_snu = mean(MF.nu,2,"omitnan");
mp_snu = mean(MP.nu,2,"omitnan");
md_snu = mean(MD.nu,2,"omitnan");
lp_snu = mean(LP.nu,2,"omitnan");
ld_snu = mean(LD.nu,2,"omitnan");

%% Annual means
nyr = nt/12;
st=1:12:nt;
en=12:12:nt;

for n=1:length(st)
    % mean nu
    sf_anu(:,n)=mean(SF.nu(:,st(n):en(n)),2,"omitnan");
    sp_anu(:,n)=mean(SP.nu(:,st(n):en(n)),2,"omitnan");
    sd_anu(:,n)=mean(SD.nu(:,st(n):en(n)),2,"omitnan");
    mf_anu(:,n)=mean(MF.nu(:,st(n):en(n)),2,"omitnan");
    mp_anu(:,n)=mean(MP.nu(:,st(n):en(n)),2,"omitnan");
    md_anu(:,n)=mean(MD.nu(:,st(n):en(n)),2,"omitnan");
    lp_anu(:,n)=mean(LP.nu(:,st(n):en(n)),2,"omitnan");
    ld_anu(:,n)=mean(LD.nu(:,st(n):en(n)),2,"omitnan");

end

%%
save([fpath 'Time_Means_FOSI_' harv cfile '.mat'],...
    'sf_tnu','sp_tnu','sd_tnu',...
    'mf_tnu','mp_tnu','md_tnu',...
    'lp_tnu','ld_tnu','-append')

save([fpath 'Space_Means_FOSI_' harv cfile '.mat'],...
    'sf_snu','sp_snu','sd_snu',...
    'mf_snu','mp_snu','md_snu',...
    'lp_snu','ld_snu','-append')

save([fpath 'Annual_Means_FOSI_' harv cfile '.mat'],...
    'sf_anu','sp_anu','sd_anu',...
    'mf_anu','mp_anu','md_anu',...
    'lp_anu','ld_anu','-append')

%%
mo = (1:nt)/12;
figure
plot(mo,(lp_tnu),'b'); hold on;
plot(mo,(mf_tnu),'r'); hold on;
plot(mo,(ld_tnu),'k'); hold on;
title('nu')




