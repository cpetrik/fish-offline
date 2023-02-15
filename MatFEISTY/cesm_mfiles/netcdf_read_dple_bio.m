% FEISTY output at all locations

clear 
close all

%%
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/DPLE/'];

%pick year
StartYr = 1954;
%loop over members
submem = 1:40;

%%
for mem=1:length(submem) %will loop over
    Member = submem(mem);
    harv = ['v15_Y' num2str(StartYr) '_M' num2str(Member) '_All_fish03_' ];
    
    %% SP
    ncid = netcdf.open([fpath 'DPLE_' harv 'sml_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);
    
    %%
    [ni,nt] = size(biomass);
    
    SP.bio = biomass;
    SP.prod = prod;
    clear biomass prod
    
    %% SF
    ncid = netcdf.open([fpath 'DPLE_' harv 'sml_f.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);
    
    SF.bio = biomass(:,1:nt);
    SF.prod = prod;
    clear biomass prod
    
    % SD
    ncid = netcdf.open([fpath 'DPLE_' harv 'sml_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);
    
    SD.bio = biomass;
    SD.prod = prod;
    clear biomass prod
    
    % MP
    ncid = netcdf.open([fpath 'DPLE_' harv 'med_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);
    
    MP.bio = biomass;
    MP.prod = prod;
    MP.yield = yield;
    clear yield biomass prod
    
    % MF
    ncid = netcdf.open([fpath 'DPLE_' harv 'med_f.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);
    
    MF.bio = biomass;
    MF.prod = prod;
    MF.yield = yield;
    clear yield biomass prod
    
    % MD
    ncid = netcdf.open([fpath 'DPLE_' harv 'med_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);
    
    MD.bio = biomass;
    MD.prod = prod;
    MD.yield = yield;
    clear yield biomass prod
    
    % LP
    ncid = netcdf.open([fpath 'DPLE_' harv 'lrg_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);
    
    LP.bio = biomass;
    LP.prod = prod;
    LP.yield = yield;
    clear yield biomass prod
    
    % LD
    ncid = netcdf.open([fpath 'DPLE_' harv 'lrg_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);
    
    LD.bio = biomass;
    LD.prod = prod;
    LD.yield = yield;
    clear yield biomass prod
    
    % Benthic material
    ncid = netcdf.open([fpath 'DPLE_' harv 'bent.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);
    
    Bent.bio = biomass;
    clear biomass
    
    %% MZ loss
    ncid = netcdf.open([fpath 'DPLE_' harv 'mzoo.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);
    
    MZ.frac = fraction;
    clear fraction
    
    %% Catch
    % Totals only in lmes
    cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
    load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
    load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
    load([cpath 'LME-mask-POP_gx1v6.mat']);
    ID = GRD.ID;
    
    tlme = double(lme_mask);
    tlme(tlme<0) = nan;
    tlme(~isnan(tlme)) = 1;
    lme_grid = tlme(ID);
    
    %TAREA units 'cm^2'
    AREA_OCN = TAREA * 1e-4;
    area = AREA_OCN(ID);
    area_km2 = area * 1e-6;
    
    MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
    nyr = nt/12;
    
    % initialized in Nov
    mos = nan*ones(ni,nt);
    mos(:,1:2) = repmat(MNTH(11:12),ni,1);
    mos(:,3:end) = repmat(MNTH,ni,floor(nyr));
    
    mns = repmat(MNTH,ni,1);
    
    area_mat = repmat(area_km2,1,nt);
    lme_mat = repmat(lme_grid,1,nt);
    
    % Units
    units_yield = 'g_m2_day';
    units_catch = 'g_km2_mo';
    
    MF.catch = MF.yield .*mos .*area_mat .*lme_mat;
    MP.catch = MP.yield .*mos .*area_mat .*lme_mat;
    MD.catch = MD.yield .*mos .*area_mat .*lme_mat;
    LP.catch = LP.yield .*mos .*area_mat .*lme_mat;
    LD.catch = LD.yield .*mos .*area_mat .*lme_mat;
    
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
    mz_tmfrac =nanmean(MZ.frac,1);
    
    sp_tprod = nanmean(SP.prod,1);
    sf_tprod = nanmean(SF.prod,1);
    sd_tprod = nanmean(SD.prod,1);
    mp_tprod = nanmean(MP.prod,1);
    mf_tprod = nanmean(MF.prod,1);
    md_tprod = nanmean(MD.prod,1);
    lp_tprod = nanmean(LP.prod,1);
    ld_tprod = nanmean(LD.prod,1);
    
    %mean yield per mo
    mf_tmy=nanmean(MF.yield,1);
    mp_tmy=nanmean(MP.yield,1);
    md_tmy=nanmean(MD.yield,1);
    lp_tmy=nanmean(LP.yield,1);
    ld_tmy=nanmean(LD.yield,1);
    %mean catch per mo
    mf_tmc=nanmean(MF.catch,1);
    mp_tmc=nanmean(MP.catch,1);
    md_tmc=nanmean(MD.catch,1);
    lp_tmc=nanmean(LP.catch,1);
    ld_tmc=nanmean(LD.catch,1);
    
    %total yield per mo
    mf_tty=nansum(MF.yield,1);
    mp_tty=nansum(MP.yield,1);
    md_tty=nansum(MD.yield,1);
    lp_tty=nansum(LP.yield,1);
    ld_tty=nansum(LD.yield,1);
    %total catch per mo
    mf_ttc=nansum(MF.catch,1);
    mp_ttc=nansum(MP.catch,1);
    md_ttc=nansum(MD.catch,1);
    lp_ttc=nansum(LP.catch,1);
    ld_ttc=nansum(LD.catch,1);
    
    %% Space
    %exclude weird jump at yr 43 (mo=517)
    % tid = 1:516;
    % sp_sbio = nanmean(SP.bio(:,tid),2);
    % sf_sbio = nanmean(SF.bio(:,tid),2);
    % sd_sbio = nanmean(SD.bio(:,tid),2);
    % mp_sbio = nanmean(MP.bio(:,tid),2);
    % mf_sbio = nanmean(MF.bio(:,tid),2);
    % md_sbio = nanmean(MD.bio(:,tid),2);
    % lp_sbio = nanmean(LP.bio(:,tid),2);
    % ld_sbio = nanmean(LD.bio(:,tid),2);
    % b_sbio  = nanmean(Bent.bio(:,tid),2);
    % mz_smfrac= nanmean(MZ.frac(:,tid),2);
    
    sp_sbio = nanmean(SP.bio,2);
    sf_sbio = nanmean(SF.bio,2);
    sd_sbio = nanmean(SD.bio,2);
    mp_sbio = nanmean(MP.bio,2);
    mf_sbio = nanmean(MF.bio,2);
    md_sbio = nanmean(MD.bio,2);
    lp_sbio = nanmean(LP.bio,2);
    ld_sbio = nanmean(LD.bio,2);
    b_sbio  = nanmean(Bent.bio,2);
    mz_smfrac= nanmean(MZ.frac,2);
    
    sp_sprod = nanmean(SP.prod,2);
    sf_sprod = nanmean(SF.prod,2);
    sd_sprod = nanmean(SD.prod,2);
    mp_sprod = nanmean(MP.prod,2);
    mf_sprod = nanmean(MF.prod,2);
    md_sprod = nanmean(MD.prod,2);
    lp_sprod = nanmean(LP.prod,2);
    ld_sprod = nanmean(LD.prod,2);
    
    %mean yield per mo
    mf_smy=nanmean(MF.yield,2);
    mp_smy=nanmean(MP.yield,2);
    md_smy=nanmean(MD.yield,2);
    lp_smy=nanmean(LP.yield,2);
    ld_smy=nanmean(LD.yield,2);
    %mean catch per mo
    mf_smc=nanmean(MF.catch,2);
    mp_smc=nanmean(MP.catch,2);
    md_smc=nanmean(MD.catch,2);
    lp_smc=nanmean(LP.catch,2);
    ld_smc=nanmean(LD.catch,2);
    
    %total yield per mo
    mf_sty=nansum(MF.yield,2);
    mp_sty=nansum(MP.yield,2);
    md_sty=nansum(MD.yield,2);
    lp_sty=nansum(LP.yield,2);
    ld_sty=nansum(LD.yield,2);
    %total catch per mo
    mf_stc=nansum(MF.catch,2);
    mp_stc=nansum(MP.catch,2);
    md_stc=nansum(MD.catch,2);
    lp_stc=nansum(LP.catch,2);
    ld_stc=nansum(LD.catch,2);
    
    %% Total times overcon happens
    MZ.over = nan*ones(size(MZ.frac));
    MZ.over(MZ.frac > 1) = ones;
    MZ.over(MZ.frac <= 1) = zeros;
    % Time
    mz_ttf=nansum(MZ.over,1);
    % Space
    mz_stf=nansum(MZ.over,2);
    
    %% Annual means
    nyr = floor(nt/12);
    % initialized in Nov
    st=3:12:length(time);
    en=14:12:length(time);
%     mz_mtf = nan*ones(ni,nyr);
%     mf_tac = nan*ones(ni,nyr);
%     mp_tac = nan*ones(ni,nyr);
%     md_tac = nan*ones(ni,nyr);
%     lp_tac = nan*ones(ni,nyr);
%     ld_tac = nan*ones(ni,nyr);
%     
%     for n=1:length(st)
%         % total overcon
%         mz_mtf(:,n)=nansum(MZ.over(:,st(n):en(n)),2);
%         
%         % mean biomass
%         sp_abio(:,n)=nanmean(SP.bio(:,st(n):en(n)),2);
%         sf_abio(:,n)=nanmean(SF.bio(:,st(n):en(n)),2);
%         sd_abio(:,n)=nanmean(SD.bio(:,st(n):en(n)),2);
%         mp_abio(:,n)=nanmean(MP.bio(:,st(n):en(n)),2);
%         mf_abio(:,n)=nanmean(MF.bio(:,st(n):en(n)),2);
%         md_abio(:,n)=nanmean(MD.bio(:,st(n):en(n)),2);
%         lp_abio(:,n)=nanmean(LP.bio(:,st(n):en(n)),2);
%         ld_abio(:,n)=nanmean(LD.bio(:,st(n):en(n)),2);
%         b_abio(:,n)=nanmean(Bent.bio(:,st(n):en(n)),2);
%         
%         % mean prod
%         sp_aprod(:,n)=nanmean(SP.prod(:,st(n):en(n)),2);
%         sf_aprod(:,n)=nanmean(SF.prod(:,st(n):en(n)),2);
%         sd_aprod(:,n)=nanmean(SD.prod(:,st(n):en(n)),2);
%         mp_aprod(:,n)=nanmean(MP.prod(:,st(n):en(n)),2);
%         mf_aprod(:,n)=nanmean(MF.prod(:,st(n):en(n)),2);
%         md_aprod(:,n)=nanmean(MD.prod(:,st(n):en(n)),2);
%         lp_aprod(:,n)=nanmean(LP.prod(:,st(n):en(n)),2);
%         ld_aprod(:,n)=nanmean(LD.prod(:,st(n):en(n)),2);
%         
%         % catch
%         mp_tac(:,n)=nansum(MP.catch(:,st(n):en(n)),2);
%         mf_tac(:,n)=nansum(MF.catch(:,st(n):en(n)),2);
%         md_tac(:,n)=nansum(MD.catch(:,st(n):en(n)),2);
%         lp_tac(:,n)=nansum(LP.catch(:,st(n):en(n)),2);
%         ld_tac(:,n)=nansum(LD.catch(:,st(n):en(n)),2);
%         
%     end
%     
%     tmn = mf_tac + mp_tac + md_tac + lp_tac + ld_tac;
%     stmn = sum(tmn);
%     
%     mp_tsac = nansum(mp_tac);
%     mf_tsac = nansum(mf_tac);
%     md_tsac = nansum(md_tac);
%     lp_tsac = nansum(lp_tac);
%     ld_tsac = nansum(ld_tac);
    
    %%
    fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/DPLE/'];
    save([fpath 'Time_Means_DPLE_' harv cfile '.mat'],'time',...
        'sf_tmean','sp_tmean','sd_tmean',...
        'mf_tmean','mp_tmean','md_tmean',...
        'lp_tmean','ld_tmean','b_tmean',...
        'sf_tprod','sp_tprod','sd_tprod',...
        'mf_tprod','mp_tprod','md_tprod',...
        'lp_tprod','ld_tprod','b_tmean',...
        'mz_tmfrac','mz_ttf','units_yield','units_catch',...
        'mf_tmy','mp_tmy','md_tmy','lp_tmy','ld_tmy',...
        'mf_tmc','mp_tmc','md_tmc','lp_tmc','ld_tmc',...
        'mf_tty','mp_tty','md_tty','lp_tty','ld_tty',...
        'mf_ttc','mp_ttc','md_ttc','lp_ttc','ld_ttc')
    
    save([fpath 'Space_Means_DPLE_' harv cfile '.mat'],'time',...
        'sf_sbio','sp_sbio','sd_sbio',...
        'mf_sbio','mp_sbio','md_sbio',...
        'lp_sbio','ld_sbio','b_sbio',...
        'sf_sprod','sp_sprod','sd_sprod',...
        'mf_sprod','mp_sprod','md_sprod',...
        'lp_sprod','ld_sprod',...
        'mz_smfrac','mz_stf','units_yield','units_catch',...
        'mf_smy','mp_smy','md_smy','lp_smy','ld_smy',...
        'mf_smc','mp_smc','md_smc','lp_smc','ld_smc',...
        'mf_sty','mp_sty','md_sty','lp_sty','ld_sty',...
        'mf_stc','mp_stc','md_stc','lp_stc','ld_stc')
    
%     save([fpath 'Annual_Means_DPLE_' harv cfile '.mat'],'time',...
%         'sf_abio','sp_abio','sd_abio',...
%         'mf_abio','mp_abio','md_abio',...
%         'lp_abio','ld_abio','b_abio',...
%         'sf_aprod','sp_aprod','sd_aprod',...
%         'mf_aprod','mp_aprod','md_aprod',...
%         'lp_aprod','ld_aprod',...
%         'mz_mtf','mf_tac','mp_tac','md_tac','lp_tac','ld_tac',...
%         'mf_tsac','mp_tsac','md_tsac','lp_tsac','ld_tsac',...
%         'units_yield','units_catch')
    
end
