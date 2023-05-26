% Calc LME production of FEISTY
% CESM FOSI

clear 
close all

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod = 'v15_All_fish03_';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
dpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([dpath 'Space_Means_FOSI_' mod cfile '.mat'],...
    'sf_sprod','sp_sprod','sd_sprod',...
    'mf_sprod','mp_sprod','md_sprod',...
    'lp_sprod','ld_sprod');

%% Map data
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%% Put biomass on grid
Csf=NaN*ones(ni,nj);
Csp=NaN*ones(ni,nj);
Csd=NaN*ones(ni,nj);
Cmf=NaN*ones(ni,nj);
Cmp=NaN*ones(ni,nj);
Cmd=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);

Csf(ID)=sf_sprod;
Csp(ID)=sp_sprod;
Csd(ID)=sd_sprod;
Cmf(ID)=mf_sprod;
Cmp(ID)=mp_sprod;
Cmd(ID)=md_sprod;
Clp(ID)=lp_sprod;
Cld(ID)=ld_sprod;

CF = Csf+Cmf;
CP = Csp+Cmp+Clp;
CD = Csd+Cmd+Cld;
CS = Csp+Csf+Csd;
CM = Cmp+Cmf+Cmd;
CL = Clp+Cld;

% g/m2 --> total g
AREA_OCN = TAREA * 1e-4;
Asf_mean = Csf .* AREA_OCN;
Asp_mean = Csp .* AREA_OCN;
Asd_mean = Csd .* AREA_OCN;
Amf_mean = Cmf .* AREA_OCN;
Amp_mean = Cmp .* AREA_OCN;
Amd_mean = Cmd .* AREA_OCN;
Alp_mean = Clp .* AREA_OCN;
Ald_mean = Cld .* AREA_OCN;

%% Calc LMEs
tlme = double(lme_mask);
tlme(tlme<0) = nan;

lme_sprod = NaN*ones(66,8);
lme_area = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    %total prod for area-weighted mean
    lme_sprod(L,1) = sum(Asf_mean(lid),'omitnan');
    lme_sprod(L,2) = sum(Asp_mean(lid),'omitnan');
    lme_sprod(L,3) = sum(Asd_mean(lid),'omitnan');
    lme_sprod(L,4) = sum(Amf_mean(lid),'omitnan');
    lme_sprod(L,5) = sum(Amp_mean(lid),'omitnan');
    lme_sprod(L,6) = sum(Amd_mean(lid),'omitnan');
    lme_sprod(L,7) = sum(Alp_mean(lid),'omitnan');
    lme_sprod(L,8) = sum(Ald_mean(lid),'omitnan');

    %LME area
    lme_area(L,1) = sum(AREA_OCN(lid),'omitnan');
end

%% Change to area-weighted means
lme_area_mat = repmat(lme_area,1,8);
lme_mprod = lme_sprod ./ lme_area_mat;

%%
save([dpath 'LME_fosi_fished_',mod,cfile '.mat'],...
    'lme_mprod','lme_area','-append');
