% Calc LME biomass of FEISTY
% CESM FOSI

clear all
close all

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
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

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_noCC_RE00100';

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
dpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

sims = {'v13_All_fish03_';'v13_climatol_';'v13_varFood_';'v13_varTemp_'};
    
for n=1:length(sims)
    close all
    mod = sims{n};
    
load([dpath 'Annual_Means_FOSI_' mod cfile '.mat'],'sf_abio','sp_abio','sd_abio',...
    'mf_abio','mp_abio','md_abio',...
    'lp_abio','ld_abio','b_abio');

[nid,nt] = size(ld_abio);

%% g/m2 --> total g
AREA = TAREA(ID) * 1e-4;
AREA_OCN = repmat(AREA,1,nt);
Asf_mean = sf_abio .* AREA_OCN;
Asp_mean = sp_abio .* AREA_OCN;
Asd_mean = sd_abio .* AREA_OCN;
Amf_mean = mf_abio .* AREA_OCN;
Amp_mean = mp_abio .* AREA_OCN;
Amd_mean = md_abio .* AREA_OCN;
Alp_mean = lp_abio .* AREA_OCN;
Ald_mean = ld_abio .* AREA_OCN;
Ab_mean  = b_abio  .* AREA_OCN;

%% Calc LMEs
glme = double(lme_mask);
glme(glme<0) = nan;
tlme = glme(ID);

lme_msfb = NaN*ones(66,nt);
lme_ssfb = NaN*ones(66,nt);
lme_mspb = lme_msfb;
lme_msdb = lme_msfb;
lme_mmpb = lme_msfb;
lme_mmdb = lme_msfb;
lme_mlpb = lme_msfb;
lme_mldb = lme_msfb;
lme_mmfb = lme_msfb;
lme_mbb = lme_msfb;
lme_sspb = lme_msfb;
lme_ssdb = lme_msfb;
lme_smpb = lme_msfb;
lme_smdb = lme_msfb;
lme_slpb = lme_msfb;
lme_sldb = lme_msfb;
lme_smfb = lme_msfb;
lme_sbb = lme_msfb;

for L=1:66
    lid = find(tlme==L);
    if (~isempty(lid))
        %mean biomass
        lme_msfb(L,:) = nanmean(Asf_mean(lid,:));
        lme_mspb(L,:) = nanmean(Asp_mean(lid,:));
        lme_msdb(L,:) = nanmean(Asd_mean(lid,:));
        lme_mmfb(L,:) = nanmean(Amf_mean(lid,:));
        lme_mmpb(L,:) = nanmean(Amp_mean(lid,:));
        lme_mmdb(L,:) = nanmean(Amd_mean(lid,:));
        lme_mlpb(L,:) = nanmean(Alp_mean(lid,:));
        lme_mldb(L,:) = nanmean(Ald_mean(lid,:));
        lme_mbb(L,:) = nanmean(Ab_mean(lid,:));
        %total biomass
        lme_ssfb(L,:) = nansum(Asf_mean(lid,:));
        lme_sspb(L,:) = nansum(Asp_mean(lid,:));
        lme_ssdb(L,:) = nansum(Asd_mean(lid,:));
        lme_smfb(L,:) = nansum(Amf_mean(lid,:));
        lme_smpb(L,:) = nansum(Amp_mean(lid,:));
        lme_smdb(L,:) = nansum(Amd_mean(lid,:));
        lme_slpb(L,:) = nansum(Alp_mean(lid,:));
        lme_sldb(L,:) = nansum(Ald_mean(lid,:));
        lme_sbb(L,:) = nansum(Ab_mean(lid,:));
    end
end

%%
save([dpath 'LME_fosi_fished_',mod,cfile '.mat'],...
    'lme_msfb','lme_mspb','lme_msdb','lme_mmfb','lme_mmpb','lme_mmdb',...
    'lme_mlpb','lme_mldb','lme_mbb',...
    'lme_ssfb','lme_sspb','lme_ssdb','lme_smfb','lme_smpb','lme_smdb',...
    'lme_slpb','lme_sldb','lme_sbb');

end


