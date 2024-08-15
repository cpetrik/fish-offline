% FEISTY output at all locations

clear 
close all

%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
harv = 'v15_All_fish03_';


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
clear nu

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
clear nu

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
clear nu

%% Take means for visualization

%Time
mf_tnu = nanmean(MF.nu,1);
lp_tnu = nanmean(LP.nu,1);
ld_tnu = nanmean(LD.nu,1);

%% Space
mf_snu = nanmean(MF.nu,2);
lp_snu = nanmean(LP.nu,2);
ld_snu = nanmean(LD.nu,2);

%% Annual means
nyr = nt/12;
st=1:12:nt;
en=12:12:nt;

for n=1:length(st)
    % mean nu
    mf_anu(:,n)=nanmean(MF.nu(:,st(n):en(n)),2);
    lp_anu(:,n)=nanmean(LP.nu(:,st(n):en(n)),2);
    ld_anu(:,n)=nanmean(LD.nu(:,st(n):en(n)),2);

end

%%
save([fpath 'Time_Means_FOSI_' harv cfile '.mat'],...
    'mf_tnu','lp_tnu','ld_tnu','-append')

save([fpath 'Space_Means_FOSI_' harv cfile '.mat'],...
    'mf_snu','lp_snu','ld_snu','-append')

save([fpath 'Annual_Means_FOSI_' harv cfile '.mat'],...
    'mf_anu','lp_anu','ld_anu','-append')

%%
mo = (1:nt)/12;
figure
plot(mo,(lp_tnu),'b'); hold on;
plot(mo,(mf_tnu),'r'); hold on;
plot(mo,(ld_tnu),'k'); hold on;

figure
plot(mo,5.5902e3*(lp_tnu),'b'); hold on;
plot(mo,11.1803*(mf_tnu),'r'); hold on;
plot(mo,5.5902e3*(ld_tnu),'k'); hold on;


