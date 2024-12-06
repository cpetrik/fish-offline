% FEISTY output at all locations

clear 
close all

%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
harv = 'v15_All_fish03_';

%% SF
ncid = netcdf.open([fpath 'FOSI_' harv 'rec_gamma_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(rec);

SF.rec = rec;
clear rec 

% SP
ncid = netcdf.open([fpath 'FOSI_' harv 'rec_gamma_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SP.rec = rec;
clear rec 

% SD
ncid = netcdf.open([fpath 'FOSI_' harv 'rec_gamma_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.rec = rec;
clear rec 

%% MF
ncid = netcdf.open([fpath 'FOSI_' harv 'rec_gamma_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(rec);

MF.rec = rec;
MF.gamma = gamma;
clear rec gamma

% MP
ncid = netcdf.open([fpath 'FOSI_' harv 'rec_gamma_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.rec = rec;
MP.gamma = gamma;
clear rec gamma

% MD
ncid = netcdf.open([fpath 'FOSI_' harv 'rec_gamma_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.rec = rec;
MD.gamma = gamma;
clear rec gamma


%% LP
ncid = netcdf.open([fpath 'FOSI_' harv 'rec_gamma_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.rec = rec;
LP.gamma = gamma;
clear rec gamma

% LD
ncid = netcdf.open([fpath 'FOSI_' harv 'rec_gamma_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.rec = rec;
LD.gamma = gamma;
clear rec gamma

%% Take means for visualization

%Time
sf_trec = mean(SF.rec,1,"omitnan");
sp_trec = mean(SP.rec,1,"omitnan");
sd_trec = mean(SD.rec,1,"omitnan");
mf_trec = mean(MF.rec,1,"omitnan");
mp_trec = mean(MP.rec,1,"omitnan");
md_trec = mean(MD.rec,1,"omitnan");
lp_trec = mean(LP.rec,1,"omitnan");
ld_trec = mean(LD.rec,1,"omitnan");

mf_tgam = mean(MF.gamma,1,"omitnan");
mp_tgam = mean(MP.gamma,1,"omitnan");
md_tgam = mean(MD.gamma,1,"omitnan");
lp_tgam = mean(LP.gamma,1,"omitnan");
ld_tgam = mean(LD.gamma,1,"omitnan");

%% Space
sf_srec = mean(SF.rec,2,"omitnan");
sp_srec = mean(SP.rec,2,"omitnan");
sd_srec = mean(SD.rec,2,"omitnan");
mf_srec = mean(MF.rec,2,"omitnan");
mp_srec = mean(MP.rec,2,"omitnan");
md_srec = mean(MD.rec,2,"omitnan");
lp_srec = mean(LP.rec,2,"omitnan");
ld_srec = mean(LD.rec,2,"omitnan");

mf_sgam = mean(MF.gamma,2,"omitnan");
mp_sgam = mean(MP.gamma,2,"omitnan");
md_sgam = mean(MD.gamma,2,"omitnan");
lp_sgam = mean(LP.gamma,2,"omitnan");
ld_sgam = mean(LD.gamma,2,"omitnan");

%% Annual means
nyr = nt/12;
st=1:12:nt;
en=12:12:nt;

for n=1:length(st)
    % mean rec
    sf_arec(:,n)=mean(SF.rec(:,st(n):en(n)),2,"omitnan");
    sp_arec(:,n)=mean(SP.rec(:,st(n):en(n)),2,"omitnan");
    sd_arec(:,n)=mean(SD.rec(:,st(n):en(n)),2,"omitnan");
    mf_arec(:,n)=mean(MF.rec(:,st(n):en(n)),2,"omitnan");
    mp_arec(:,n)=mean(MP.rec(:,st(n):en(n)),2,"omitnan");
    md_arec(:,n)=mean(MD.rec(:,st(n):en(n)),2,"omitnan");
    lp_arec(:,n)=mean(LP.rec(:,st(n):en(n)),2,"omitnan");
    ld_arec(:,n)=mean(LD.rec(:,st(n):en(n)),2,"omitnan");
    
    % mean gam
    mf_agam(:,n)=mean(MF.gamma(:,st(n):en(n)),2,"omitnan");
    mp_agam(:,n)=mean(MP.gamma(:,st(n):en(n)),2,"omitnan");
    md_agam(:,n)=mean(MD.gamma(:,st(n):en(n)),2,"omitnan");
    lp_agam(:,n)=mean(LP.gamma(:,st(n):en(n)),2,"omitnan");
    ld_agam(:,n)=mean(LD.gamma(:,st(n):en(n)),2,"omitnan");

end

%%
save([fpath 'Time_Means_FOSI_' harv cfile '.mat'],...
    'sf_trec','sp_trec','sd_trec',...
    'mf_trec','mp_trec','md_trec',...
    'lp_trec','ld_trec',...
    'mf_tgam','mp_tgam','md_tgam',...
    'lp_tgam','ld_tgam','-append')

save([fpath 'Space_Means_FOSI_' harv cfile '.mat'],...
    'sf_srec','sp_srec','sd_srec',...
    'mf_srec','mp_srec','md_srec',...
    'lp_srec','ld_srec',...
    'mf_sgam','mp_sgam','md_sgam',...
    'lp_sgam','ld_sgam','-append')

save([fpath 'Annual_Means_FOSI_' harv cfile '.mat'],...
    'sf_arec','sp_arec','sd_arec',...
    'mf_arec','mp_arec','md_arec',...
    'lp_arec','ld_arec',...
    'mf_agam','mp_agam','md_agam',...
    'lp_agam','ld_agam','-append')

%%
mo = (1:nt)/12;
figure
plot(mo,(lp_trec),'b'); hold on;
plot(mo,(mf_trec),'r'); hold on;
plot(mo,(ld_trec),'k'); hold on;
title('rec')

figure
plot(mo,(lp_tgam),'b'); hold on;
plot(mo,(mf_tgam),'r'); hold on;
plot(mo,(ld_tgam),'k'); hold on;
title('gamma')



