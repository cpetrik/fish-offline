% FEISTY output at all locations

clear 
close all

%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
harv = 'v15_All_fish03_';

%% SP
ncid = netcdf.open([fpath 'FOSI_' harv 'con_rep_sml_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
[ni,nt] = size(con);

SP.con = con;
clear con

%% SF
ncid = netcdf.open([fpath 'FOSI_' harv 'con_rep_sml_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SF.con = con;
clear con

% SD
ncid = netcdf.open([fpath 'FOSI_' harv 'con_rep_sml_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

SD.con = con;
clear con

% MP
ncid = netcdf.open([fpath 'FOSI_' harv 'con_rep_med_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MP.con = con;
clear con

% MD
ncid = netcdf.open([fpath 'FOSI_' harv 'con_rep_med_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

MD.con = con;
clear con

%% MF
ncid = netcdf.open([fpath 'FOSI_' harv 'con_rep_med_f.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

[ni,nt] = size(rep);

MF.rep = rep;
MF.con = con;
clear rep con

% LP
ncid = netcdf.open([fpath 'FOSI_' harv 'con_rep_lrg_p.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LP.rep = rep;
LP.con = con;
clear rep con

% LD
ncid = netcdf.open([fpath 'FOSI_' harv 'con_rep_lrg_d.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

LD.rep = rep;
LD.con = con;
clear rep con

%% Take means for visualization

%Time
mf_trep = mean(MF.rep,1,"omitnan");
lp_trep = mean(LP.rep,1,"omitnan");
ld_trep = mean(LD.rep,1,"omitnan");

sf_tcon = mean(SF.con,1,"omitnan");
sp_tcon = mean(SP.con,1,"omitnan");
sd_tcon = mean(SD.con,1,"omitnan");
mf_tcon = mean(MF.con,1,"omitnan");
mp_tcon = mean(MP.con,1,"omitnan");
md_tcon = mean(MD.con,1,"omitnan");
lp_tcon = mean(LP.con,1,"omitnan");
ld_tcon = mean(LD.con,1,"omitnan");

%% Space
mf_srep = mean(MF.rep,2,"omitnan");
lp_srep = mean(LP.rep,2,"omitnan");
ld_srep = mean(LD.rep,2,"omitnan");

sf_scon = mean(SF.con,2,"omitnan");
sp_scon = mean(SP.con,2,"omitnan");
sd_scon = mean(SD.con,2,"omitnan");
mf_scon = mean(MF.con,2,"omitnan");
mp_scon = mean(MP.con,2,"omitnan");
md_scon = mean(MD.con,2,"omitnan");
lp_scon = mean(LP.con,2,"omitnan");
ld_scon = mean(LD.con,2,"omitnan");

%% Annual means
nyr = nt/12;
st=1:12:nt;
en=12:12:nt;

for n=1:length(st)
    % mean rep
    mf_arep(:,n)=mean(MF.rep(:,st(n):en(n)),2,"omitnan");
    lp_arep(:,n)=mean(LP.rep(:,st(n):en(n)),2,"omitnan");
    ld_arep(:,n)=mean(LD.rep(:,st(n):en(n)),2,"omitnan");
    
    % mean con
    sf_acon(:,n)=mean(SF.con(:,st(n):en(n)),2,"omitnan");
    sp_acon(:,n)=mean(SP.con(:,st(n):en(n)),2,"omitnan");
    sd_acon(:,n)=mean(SD.con(:,st(n):en(n)),2,"omitnan");
    mf_acon(:,n)=mean(MF.con(:,st(n):en(n)),2,"omitnan");
    mp_acon(:,n)=mean(MP.con(:,st(n):en(n)),2,"omitnan");
    md_acon(:,n)=mean(MD.con(:,st(n):en(n)),2,"omitnan");
    lp_acon(:,n)=mean(LP.con(:,st(n):en(n)),2,"omitnan");
    ld_acon(:,n)=mean(LD.con(:,st(n):en(n)),2,"omitnan");

end

%%
save([fpath 'Time_Means_FOSI_' harv cfile '.mat'],...
    'mf_trep','lp_trep','ld_trep',...
    'sf_tcon','sp_tcon','sd_tcon',...
    'mf_tcon','mp_tcon','md_tcon',...
    'lp_tcon','ld_tcon','-append')

save([fpath 'Space_Means_FOSI_' harv cfile '.mat'],...
    'mf_srep','lp_srep','ld_srep',...
    'sf_scon','sp_scon','sd_scon',...
    'mf_scon','mp_scon','md_scon',...
    'lp_scon','ld_scon','-append')

save([fpath 'Annual_Means_FOSI_' harv cfile '.mat'],...
    'mf_arep','lp_arep','ld_arep',...
    'sf_acon','sp_acon','sd_acon',...
    'mf_acon','mp_acon','md_acon',...
    'lp_acon','ld_acon','-append')

%%
mo = (1:nt)/12;
figure
plot(mo,(lp_trep),'b'); hold on;
plot(mo,(mf_trep),'r'); hold on;
plot(mo,(ld_trep),'k'); hold on;
title('rep')

figure
plot(mo,(lp_tcon),'b'); hold on;
plot(mo,(mf_tcon),'r'); hold on;
plot(mo,(ld_tcon),'k'); hold on;
title('con')



