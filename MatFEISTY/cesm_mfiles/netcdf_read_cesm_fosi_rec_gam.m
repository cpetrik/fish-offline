% FEISTY output at all locations

clear 
close all

%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
harv = 'v15_All_fish03_';


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
mf_trec = mean(MF.rec,1,"omitnan");
lp_trec = mean(LP.rec,1,"omitnan");
ld_trec = mean(LD.rec,1,"omitnan");
mf_tgam = mean(MF.gamma,1,"omitnan");
lp_tgam = mean(LP.gamma,1,"omitnan");
ld_tgam = mean(LD.gamma,1,"omitnan");

%% Space
mf_srec = mean(MF.rec,2,"omitnan");
lp_srec = mean(LP.rec,2,"omitnan");
ld_srec = mean(LD.rec,2,"omitnan");
mf_sgam = mean(MF.gamma,2,"omitnan");
lp_sgam = mean(LP.gamma,2,"omitnan");
ld_sgam = mean(LD.gamma,2,"omitnan");

%% Annual means
nyr = nt/12;
st=1:12:nt;
en=12:12:nt;

for n=1:length(st)
    % mean rec
    mf_arec(:,n)=mean(MF.rec(:,st(n):en(n)),2,"omitnan");
    lp_arec(:,n)=mean(LP.rec(:,st(n):en(n)),2,"omitnan");
    ld_arec(:,n)=mean(LD.rec(:,st(n):en(n)),2,"omitnan");
    % mean gamma
    mf_agam(:,n)=mean(MF.gamma(:,st(n):en(n)),2,"omitnan");
    lp_agam(:,n)=mean(LP.gamma(:,st(n):en(n)),2,"omitnan");
    ld_agam(:,n)=mean(LD.gamma(:,st(n):en(n)),2,"omitnan");

end

%%
save([fpath 'Time_Means_FOSI_' harv cfile '.mat'],...
    'mf_trec','lp_trec','ld_trec',...
    'mf_tgam','lp_tgam','ld_tgam','-append')

save([fpath 'Space_Means_FOSI_' harv cfile '.mat'],...
    'mf_srec','lp_srec','ld_srec',...
    'mf_sgam','lp_sgam','ld_sgam','-append')

save([fpath 'Annual_Means_FOSI_' harv cfile '.mat'],...
    'mf_arec','lp_arec','ld_arec',...
    'mf_agam','lp_agam','ld_agam','-append')

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



