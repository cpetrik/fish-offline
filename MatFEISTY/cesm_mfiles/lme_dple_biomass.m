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

for mem=1:length(submem) %will loop over
    Member = submem(mem);
    exper = ['v14_Y' num2str(StartYr) '_M' num2str(Member) '_All_fish03_' ];
    
    load([fpath 'Space_Means_DPLE_' exper cfile '.mat'],...
        'sf_sbio','sp_sbio','sd_sbio',...
        'mf_sbio','mp_sbio','md_sbio',...
        'lp_sbio','ld_sbio','b_sbio');
    
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
    
    %% Calc LMEs
    tlme = double(lme_mask);
    tlme(tlme<0) = nan;
    
    lme_mbio = NaN*ones(66,9);
    lme_sbio = NaN*ones(66,9);
    
    for L=1:66
        lid = find(tlme==L);
        %mean biomass
        lme_mbio(L,1) = nanmean(Asf_mean(lid));
        lme_mbio(L,2) = nanmean(Asp_mean(lid));
        lme_mbio(L,3) = nanmean(Asd_mean(lid));
        lme_mbio(L,4) = nanmean(Amf_mean(lid));
        lme_mbio(L,5) = nanmean(Amp_mean(lid));
        lme_mbio(L,6) = nanmean(Amd_mean(lid));
        lme_mbio(L,7) = nanmean(Alp_mean(lid));
        lme_mbio(L,8) = nanmean(Ald_mean(lid));
        lme_mbio(L,9) = nanmean(Ab_mean(lid));
        %total biomass
        lme_sbio(L,1) = nansum(Asf_mean(lid));
        lme_sbio(L,2) = nansum(Asp_mean(lid));
        lme_sbio(L,3) = nansum(Asd_mean(lid));
        lme_sbio(L,4) = nansum(Amf_mean(lid));
        lme_sbio(L,5) = nansum(Amp_mean(lid));
        lme_sbio(L,6) = nansum(Amd_mean(lid));
        lme_sbio(L,7) = nansum(Alp_mean(lid));
        lme_sbio(L,8) = nansum(Ald_mean(lid));
        lme_sbio(L,9) = nansum(Ab_mean(lid));
    end
    
    %%
    save([fpath 'LME_DPLE_fished_',exper,cfile '.mat'],...
        'lme_mbio','lme_sbio','-append');
    
end


