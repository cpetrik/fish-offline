% FEISTY output at all locations
% 1st 4 yrs run on Python code

clear all
close all

%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ100_mMZ045_nmort1_BE08_noCC_RE00100';

fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
harv = 'cesm.4_yr';

%%
ncdisp([fpath 'FOSI_' harv '.nc'])

group = {'Sf','Sp','Sd','Mf','Mp','Md','Lp','Ld','B'};

%% All variables together
ncid = netcdf.open([fpath 'FOSI_' harv '.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
% skip variable 2 = group
for i = 3:nvars
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 99999) = NaN;']);
end
netcdf.close(ncid);

%%
[ni,ng,nt] = size(biomass);

SF.bio = squeeze(biomass(:,1,:));
SP.bio = squeeze(biomass(:,2,:));
SD.bio = squeeze(biomass(:,3,:));
MF.bio = squeeze(biomass(:,4,:));
MP.bio = squeeze(biomass(:,5,:));
MD.bio = squeeze(biomass(:,6,:));
LP.bio = squeeze(biomass(:,7,:));
LD.bio = squeeze(biomass(:,8,:));
Bent.bio = squeeze(biomass(:,9,:));

%clear biomass

%% Take means for visualization

%Time
sp_tmean = nanmean(SP.bio,1);
sf_tmean = nanmean(SF.bio,1);
sd_tmean = nanmean(SD.bio,1);
mp_tmean = nanmean(MP.bio,1);
mf_tmean = nanmean(MF.bio,1);
md_tmean = nanmean(MD.bio,1);
lp_tmean = nanmean(LP.bio,1);
ld_tmean = nanmean(LD.bio,1);
b_tmean  = nanmean(Bent.bio,1);

%% Space
sp_sbio = nanmean(SP.bio,2);
sf_sbio = nanmean(SF.bio,2);
sd_sbio = nanmean(SD.bio,2);
mp_sbio = nanmean(MP.bio,2);
mf_sbio = nanmean(MF.bio,2);
md_sbio = nanmean(MD.bio,2);
lp_sbio = nanmean(LP.bio,2);
ld_sbio = nanmean(LD.bio,2);
b_sbio  = nanmean(Bent.bio,2);

%% Annual means
nyr = nt/365;
st=1:365:length(time);
en=365:365:length(time);

for n=1:length(st)
    % mean biomass
    sp_abio(:,n)=nanmean(SP.bio(:,st(n):en(n)),2);
    sf_abio(:,n)=nanmean(SF.bio(:,st(n):en(n)),2);
    sd_abio(:,n)=nanmean(SD.bio(:,st(n):en(n)),2);
    mp_abio(:,n)=nanmean(MP.bio(:,st(n):en(n)),2);
    mf_abio(:,n)=nanmean(MF.bio(:,st(n):en(n)),2);
    md_abio(:,n)=nanmean(MD.bio(:,st(n):en(n)),2);
    lp_abio(:,n)=nanmean(LP.bio(:,st(n):en(n)),2);
    ld_abio(:,n)=nanmean(LD.bio(:,st(n):en(n)),2);
    b_abio(:,n)=nanmean(Bent.bio(:,st(n):en(n)),2);

end

%%
save([fpath 'Time_Means_FOSI_' harv,'_',cfile '.mat'],'time',...
    'sf_tmean','sp_tmean','sd_tmean',...
    'mf_tmean','mp_tmean','md_tmean',...
    'lp_tmean','ld_tmean','b_tmean')

save([fpath 'Space_Means_FOSI_' harv,'_',cfile '.mat'],'time',...
    'sf_sbio','sp_sbio','sd_sbio',...
    'mf_sbio','mp_sbio','md_sbio',...
    'lp_sbio','ld_sbio','b_sbio')

save([fpath 'Annual_Means_FOSI_' harv,'_',cfile '.mat'],'time',...
    'sf_abio','sp_abio','sd_abio',...
    'mf_abio','mp_abio','md_abio',...
    'lp_abio','ld_abio','b_abio')

%%
mo = time/365;
figure
plot(mo,log10(lp_tmean),'b'); hold on;
plot(mo,log10(mf_tmean),'r'); hold on;
plot(mo,log10(ld_tmean),'k'); hold on;


