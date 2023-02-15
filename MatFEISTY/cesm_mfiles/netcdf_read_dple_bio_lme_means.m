% FEISTY output at all locations

clear
close all

%%
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/DPLE/'];

%pick year
StartYr = 1954;
%loop over members
submem = 1:40;

%%
for mem=1:length(submem) %will loop over
    Member = submem(mem);
    harv = ['v15_Y' num2str(StartYr) '_M' num2str(Member) '_All_fish03_' ];

    %% SP
    ncid = netcdf.open([spath 'DPLE_' harv 'sml_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    [ni,nt] = size(biomass);

    SP.bio = biomass;
    %     SP.prod = prod;
    clear biomass prod

    %% SF
    ncid = netcdf.open([spath 'DPLE_' harv 'sml_f.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    SF.bio = biomass(:,1:nt);
    %     SF.prod = prod;
    clear biomass prod

    % SD
    ncid = netcdf.open([spath 'DPLE_' harv 'sml_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    SD.bio = biomass;
    %     SD.prod = prod;
    clear biomass prod

    % MP
    ncid = netcdf.open([spath 'DPLE_' harv 'med_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MP.bio = biomass;
    %     MP.prod = prod;
    %     MP.yield = yield;
    clear yield biomass prod

    % MF
    ncid = netcdf.open([spath 'DPLE_' harv 'med_f.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MF.bio = biomass;
    %     MF.prod = prod;
    %     MF.yield = yield;
    clear yield biomass prod

    % MD
    ncid = netcdf.open([spath 'DPLE_' harv 'med_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    MD.bio = biomass;
    %     MD.prod = prod;
    %     MD.yield = yield;
    clear yield biomass prod

    % LP
    ncid = netcdf.open([spath 'DPLE_' harv 'lrg_p.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LP.bio = biomass;
    %     LP.prod = prod;
    %     LP.yield = yield;
    clear yield biomass prod

    % LD
    ncid = netcdf.open([spath 'DPLE_' harv 'lrg_d.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    LD.bio = biomass;
    %     LD.prod = prod;
    %     LD.yield = yield;
    clear yield biomass prod

    % Benthic material
    ncid = netcdf.open([spath 'DPLE_' harv 'bent.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:nvars
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 99999) = NaN;']);
    end
    netcdf.close(ncid);

    Bent.bio = biomass;
    clear biomass

    %% LME
    cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
    load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
    load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
    load([cpath 'LME-mask-POP_gx1v6.mat']);
    ID = GRD.ID;

    tlme = double(lme_mask);
    tlme(tlme<0) = nan;
    lme_grid = tlme(ID);


    %%
    lme_mSP = NaN*ones(66,nt,40);
    lme_mSF = lme_mSP; lme_mSD = lme_mSP;
    lme_mMF = lme_mSP; lme_mMP = lme_mSP; lme_mMD = lme_mSP;
    lme_mLP = lme_mSP; lme_mLD = lme_mSP; lme_mB = lme_mSP;
    lme_tSP = NaN*ones(66,nt,40);
    lme_tSF = lme_tSP; lme_tSD = lme_tSP;
    lme_tMF = lme_tSP; lme_tMP = lme_tSP; lme_tMD = lme_tSP;
    lme_tLP = lme_tSP; lme_tLD = lme_tSP; lme_tB = lme_tSP;

    for L=1:66
        lid = find(lme_grid==L);
        if (~isempty(lid))
        %mean biomass
        lme_mSF(L,:,mem) = nanmean(SF.bio(lid,:));
        lme_mSP(L,:,mem) = nanmean(SP.bio(lid,:));
        lme_mSD(L,:,mem) = nanmean(SD.bio(lid,:));
        lme_mMF(L,:,mem) = nanmean(MF.bio(lid,:));
        lme_mMP(L,:,mem) = nanmean(MP.bio(lid,:));
        lme_mMD(L,:,mem) = nanmean(MD.bio(lid,:));
        lme_mLP(L,:,mem) = nanmean(LP.bio(lid,:));
        lme_mLD(L,:,mem) = nanmean(LD.bio(lid,:));
        lme_mB(L,:,mem)  = nanmean(Bent.bio(lid,:));
        %total biomass
        lme_tSF(L,:,mem) = nansum(SF.bio(lid,:));
        lme_tSP(L,:,mem) = nansum(SP.bio(lid,:));
        lme_tSD(L,:,mem) = nansum(SD.bio(lid,:));
        lme_tMF(L,:,mem) = nansum(MF.bio(lid,:));
        lme_tMP(L,:,mem) = nansum(MP.bio(lid,:));
        lme_tMD(L,:,mem) = nansum(MD.bio(lid,:));
        lme_tLP(L,:,mem) = nansum(LP.bio(lid,:));
        lme_tLD(L,:,mem) = nansum(LD.bio(lid,:));
        lme_tB(L,:,mem)  = nansum(Bent.bio(lid,:));
        end
    end

end
%%
mod = ['v15_Y' num2str(StartYr) '_All_fish03' ];

dpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/DPLE/'];
save([dpath 'Time_Means_DPLE_LME_' mod '.mat'],'time',...
    'lme_mSP','lme_mSF','lme_mSD','lme_mMF','lme_mMP','lme_mMD',...
    'lme_mLP','lme_mLD','lme_mB',...
    'lme_tSP','lme_tSF','lme_tSD','lme_tMF','lme_tMP','lme_tMD',...
    'lme_tLP','lme_tLD','lme_tB')

