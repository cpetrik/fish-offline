% FEISTY CESM DPLE members by start year
% LME catches

clear all
close all

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);

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
    
    load([fpath 'Space_Means_DPLE_' exper cfile '.mat'],'units_yield','units_catch',...
        'mf_smy','mp_smy','md_smy','lp_smy','ld_smy',...
        'mf_smc','mp_smc','md_smc','lp_smc','ld_smc',...
        'mf_sty','mp_sty','md_sty','lp_sty','ld_sty',...
        'mf_stc','mp_stc','md_stc','lp_stc','ld_stc');
    
    %% Plots in space
    Cmf=NaN*ones(ni,nj);
    Cmp=NaN*ones(ni,nj);
    Cmd=NaN*ones(ni,nj);
    Clp=NaN*ones(ni,nj);
    Cld=NaN*ones(ni,nj);
    
    Cmf(GRD.ID)=mf_smy;
    Cmp(GRD.ID)=mp_smy;
    Cmd(GRD.ID)=md_smy;
    Clp(GRD.ID)=lp_smy;
    Cld(GRD.ID)=ld_smy;
    
    AllF = Cmf;
    AllP = Cmp+Clp;
    AllD = Cmd+Cld;
    AllM = Cmp+Cmf+Cmd;
    AllL = Clp+Cld;
    FracPD = AllP ./ (AllP+AllD);
    
    Tmf=NaN*ones(ni,nj);
    Tmp=NaN*ones(ni,nj);
    Tmd=NaN*ones(ni,nj);
    Tlp=NaN*ones(ni,nj);
    Tld=NaN*ones(ni,nj);
    
    Tmf(GRD.ID)=mf_sty;
    Tmp(GRD.ID)=mp_sty;
    Tmd(GRD.ID)=md_sty;
    Tlp(GRD.ID)=lp_sty;
    Tld(GRD.ID)=ld_sty;
    
    % tAllF = Tmf;
    % tAllP = Tmp+Tlp;
    % tAllD = Tmd+Tld;
    % tAllM = Tmp+Tmf+Tmd;
    % tAllL = Tlp+Tld;
    % tFracPD = tAllP ./ (tAllP+tAllD);
    
    %% Catches
    %TAREA units 'cm^2' to m2
    AREA_OCN = TAREA * 1e-4;
    
    % g/m2/d --> total g
    Amf_mcatch = Cmf .* AREA_OCN * 365; %mean fish catch per yr
    Amp_mcatch = Cmp .* AREA_OCN * 365;
    Amd_mcatch = Cmd .* AREA_OCN * 365;
    Alp_mcatch = Clp .* AREA_OCN * 365;
    Ald_mcatch = Cld .* AREA_OCN * 365;
    
    Amf_tot = Tmf .* AREA_OCN;
    Amp_tot = Tmp .* AREA_OCN;
    Amd_tot = Tmd .* AREA_OCN;
    Alp_tot = Tlp .* AREA_OCN;
    Ald_tot = Tld .* AREA_OCN;
    
    %% Calc LMEs
    tlme = double(lme_mask);
    tlme(tlme<0) = nan;
    
    lme_mcatch = NaN*ones(66,5);
    lme_tcatch = NaN*ones(66,9);
    lme_area = NaN*ones(66,1);
    
    for L=1:66
        lid = find(tlme==L);
        %total catch g
        lme_mcatch(L,1) = nansum(Amf_mcatch(lid));
        lme_mcatch(L,2) = nansum(Amp_mcatch(lid));
        lme_mcatch(L,3) = nansum(Amd_mcatch(lid));
        lme_mcatch(L,4) = nansum(Alp_mcatch(lid));
        lme_mcatch(L,5) = nansum(Ald_mcatch(lid));
        %total catch g
        lme_tcatch(L,1) = nansum(Amf_tot(lid));
        lme_tcatch(L,2) = nansum(Amp_tot(lid));
        lme_tcatch(L,3) = nansum(Amd_tot(lid));
        lme_tcatch(L,4) = nansum(Alp_tot(lid));
        lme_tcatch(L,5) = nansum(Ald_tot(lid));
        %total area of LME
        lme_area(L,1) = nansum(AREA_OCN(lid));
    end
    
    %%
    save([fpath 'LME_DPLE_fished_',exper,cfile '.mat'],...
        'lme_mcatch','lme_tcatch','lme_area');
end



