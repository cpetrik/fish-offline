% Calc LME biomass of FEISTY
% CESM FOSI

clear 
close all

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
mod = 'v15_obsfish2015_';

pp = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/';
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end
load([dpath 'Space_Means_FOSI_' mod cfile '.mat'],'sf_sbio','sp_sbio','sd_sbio',...
    'mf_sbio','mp_sbio','md_sbio',...
    'lp_sbio','ld_sbio','b_sbio');

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
Cb =NaN*ones(ni,nj);
Csf(ID)=sf_sbio;
Csp(ID)=sp_sbio;
Csd(ID)=sd_sbio;
Cmf(ID)=mf_sbio;
Cmp(ID)=mp_sbio;
Cmd(ID)=md_sbio;
Clp(ID)=lp_sbio;
Cld(ID)=ld_sbio;
Cb(ID) = b_sbio;

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
Ab_mean  = Cb  .* AREA_OCN;

AF = Asf_mean + Amf_mean;
AP = Asp_mean + Amp_mean + Alp_mean;
AD = Asd_mean + Amd_mean + Ald_mean;
All = AF + AP + AD;

%% Calc LMEs
tlme = double(lme_mask);
tlme(tlme<0) = nan;

%lme_mbio = NaN*ones(66,9);
lme_sbio = NaN*ones(66,9);
lme_area = NaN*ones(66,1);

lme_type = NaN*ones(66,4);

for L=1:66
    lid = find(tlme==L);
    %mean biomass - Change to area-weighted means
%     lme_mbio(L,1) = nanmean(Asf_mean(lid));
%     lme_mbio(L,2) = nanmean(Asp_mean(lid));
%     lme_mbio(L,3) = nanmean(Asd_mean(lid));
%     lme_mbio(L,4) = nanmean(Amf_mean(lid));
%     lme_mbio(L,5) = nanmean(Amp_mean(lid));
%     lme_mbio(L,6) = nanmean(Amd_mean(lid));
%     lme_mbio(L,7) = nanmean(Alp_mean(lid));
%     lme_mbio(L,8) = nanmean(Ald_mean(lid));
%     lme_mbio(L,9) = nanmean(Ab_mean(lid));
    
    %total biomass
    lme_sbio(L,1) = sum(Asf_mean(lid),'omitnan');
    lme_sbio(L,2) = sum(Asp_mean(lid),'omitnan');
    lme_sbio(L,3) = sum(Asd_mean(lid),'omitnan');
    lme_sbio(L,4) = sum(Amf_mean(lid),'omitnan');
    lme_sbio(L,5) = sum(Amp_mean(lid),'omitnan');
    lme_sbio(L,6) = sum(Amd_mean(lid),'omitnan');
    lme_sbio(L,7) = sum(Alp_mean(lid),'omitnan');
    lme_sbio(L,8) = sum(Ald_mean(lid),'omitnan');
    lme_sbio(L,9) = sum(Ab_mean(lid),'omitnan');

    lme_type(L,1) = sum(AF(lid),'omitnan');
    lme_type(L,2) = sum(AP(lid),'omitnan');
    lme_type(L,3) = sum(AD(lid),'omitnan');
    lme_type(L,4) = sum(All(lid),'omitnan');

    %LME area
    lme_area(L,1) = sum(AREA_OCN(lid),'omitnan');
end

%% Change to area-weighted means
lme_area_mat = repmat(lme_area,1,9);
lme_mbio = lme_sbio ./ lme_area_mat;

lme_area_mat2 = repmat(lme_area,1,4);
lme_mtype = lme_type ./ lme_area_mat2;

%%
save([dpath 'LME_fosi_fished_',mod,cfile '.mat'],...
    'lme_mbio','lme_sbio','lme_mtype','-append');
