% Calc LME gamma & rec of FEISTY
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
load([dpath 'Space_Means_FOSI_' mod, cfile '.mat'],'time',...
    'mf_srec','lp_srec','ld_srec',...
    'mf_sgam','lp_sgam','ld_sgam');

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
Cmf=NaN*ones(ni,nj);
Clp=NaN*ones(ni,nj);
Cld=NaN*ones(ni,nj);

Rmf=NaN*ones(ni,nj);
Rlp=NaN*ones(ni,nj);
Rld=NaN*ones(ni,nj);

Cmf(ID)=mf_sgam;
Clp(ID)=lp_sgam;
Cld(ID)=ld_sgam;

Rmf(ID)=mf_srec;
Rlp(ID)=lp_srec;
Rld(ID)=ld_srec;

% g/g/m2 --> total g/g (for area-weighted mean)
AREA_OCN = TAREA * 1e-4;
Amf_mean = Cmf .* AREA_OCN;
Alp_mean = Clp .* AREA_OCN;
Ald_mean = Cld .* AREA_OCN;
Rmf_mean = Rmf .* AREA_OCN;
Rlp_mean = Rlp .* AREA_OCN;
Rld_mean = Rld .* AREA_OCN;

%% Calc LMEs
tlme = double(lme_mask);
tlme(tlme<0) = nan;

lme_sgam = NaN*ones(66,3);
lme_srec = NaN*ones(66,3);
lme_area = NaN*ones(66,1);

for L=1:66
    lid = find(tlme==L);
    %total gamma for area-weighted mean
    lme_sgam(L,1) = sum(Amf_mean(lid),'omitnan');
    lme_sgam(L,2) = sum(Alp_mean(lid),'omitnan');
    lme_sgam(L,3) = sum(Ald_mean(lid),'omitnan');

    lme_srec(L,1) = sum(Rmf_mean(lid),'omitnan');
    lme_srec(L,2) = sum(Rlp_mean(lid),'omitnan');
    lme_srec(L,3) = sum(Rld_mean(lid),'omitnan');

    %LME area
    lme_area(L,1) = sum(AREA_OCN(lid),'omitnan');
end

%% Change to area-weighted means
lme_area_mat = repmat(lme_area,1,3);
lme_mgam = lme_sgam ./ lme_area_mat;
lme_mrec = lme_srec ./ lme_area_mat;

%%
save([dpath 'LME_fosi_fished_',mod,cfile '.mat'],...
    'lme_mgam','lme_mrec','-append');
