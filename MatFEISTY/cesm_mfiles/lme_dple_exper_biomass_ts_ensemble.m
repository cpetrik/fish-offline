% Calc LME biomass of FEISTY
% CESM DPLE

clear all
close all

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

AREA = TAREA(ID) * 1e-4;

glme = double(lme_mask);
glme(glme<0) = nan;
tlme = glme(ID);

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';
pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/DPLE/';
fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
ppath = [pp cfile '/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

%pick year
StartYr = 2015;
%loop over members
submem = 1:40;

lme_msfb = NaN*ones(66,10,length(submem));
lme_ssfb = NaN*ones(66,10,length(submem));
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


for mem=1:length(submem) %will loop over
    Member = submem(mem);
    exper = ['v14_Y' num2str(StartYr) '_M' num2str(Member) '_All_fish03_' ];
    
    load([fpath 'Annual_Means_DPLE_' exper cfile '.mat'],...
        'sf_abio','sp_abio','sd_abio',...
        'mf_abio','mp_abio','md_abio',...
        'lp_abio','ld_abio','b_abio');
    
    [nid,nt] = size(ld_abio);
    
    %% g/m2 --> total g
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
    
    for L=1:66
        lid = find(tlme==L);
        if (~isempty(lid))
            %mean biomass
            lme_msfb(L,:,mem) = nanmean(Asf_mean(lid,:));
            lme_mspb(L,:,mem) = nanmean(Asp_mean(lid,:));
            lme_msdb(L,:,mem) = nanmean(Asd_mean(lid,:));
            lme_mmfb(L,:,mem) = nanmean(Amf_mean(lid,:));
            lme_mmpb(L,:,mem) = nanmean(Amp_mean(lid,:));
            lme_mmdb(L,:,mem) = nanmean(Amd_mean(lid,:));
            lme_mlpb(L,:,mem) = nanmean(Alp_mean(lid,:));
            lme_mldb(L,:,mem) = nanmean(Ald_mean(lid,:));
            lme_mbb(L,:,mem) = nanmean(Ab_mean(lid,:));
            %total biomass
            lme_ssfb(L,:,mem) = nansum(Asf_mean(lid,:));
            lme_sspb(L,:,mem) = nansum(Asp_mean(lid,:));
            lme_ssdb(L,:,mem) = nansum(Asd_mean(lid,:));
            lme_smfb(L,:,mem) = nansum(Amf_mean(lid,:));
            lme_smpb(L,:,mem) = nansum(Amp_mean(lid,:));
            lme_smdb(L,:,mem) = nansum(Amd_mean(lid,:));
            lme_slpb(L,:,mem) = nansum(Alp_mean(lid,:));
            lme_sldb(L,:,mem) = nansum(Ald_mean(lid,:));
            lme_sbb(L,:,mem) = nansum(Ab_mean(lid,:));
        end
    end
    
    %%
    
end
save([fpath 'LME_DPLE_fished_ensemble_',cfile '.mat'],...
    'lme_msfb','lme_mspb','lme_msdb','lme_mmfb','lme_mmpb','lme_mmdb',...
    'lme_mlpb','lme_mldb','lme_mbb',...
    'lme_ssfb','lme_sspb','lme_ssdb','lme_smfb','lme_smpb','lme_smdb',...
    'lme_slpb','lme_sldb','lme_sbb');


