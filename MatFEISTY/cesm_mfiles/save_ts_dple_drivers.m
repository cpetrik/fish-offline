% Visualize input forcing of FEISTY
% CESM DPLE
% Time series plots and maps

clear
close all

%% FOSI
cpath = '/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

ID = GRD.ID;

tlme = double(lme_mask);
tlme(tlme<0) = nan;
lme_grid = tlme(ID);

[clim_Tp,clim_Tb,clim_POC,clim_zooC,clim_loss] = fosi_clim_for_dple(cpath);

%% DPLE
fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/DPLE/';

% GET DPLE time info ----------------------------
ncdisp([fpath 'DPLE-FIESTY-forcing_zooC_150m.nc'])

% Pelagic temperature
ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_TEMP_150m.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end

%%
[ni,nj] = size(TLONG);
mos = length(L);

% FOSI is 1948-2015
% DPLE is 1954-2017
% leadtimes = 1:24;
firstyear = 1954; %of DPLE initializations
lastyear  = 2015; %of DPLE initializations
startmonth = 11;  %of DPLE initializations

%% array for members
nt = mos;

tTP = NaN*ones(40,nt);
tTB = tTP; tDet = tTP;
tMZ = tTP; tMZl = tTP;

lme_mTP = NaN*ones(66,nt,40);
lme_mTB = lme_mTP; lme_mDet = lme_mTP;
lme_mMZ = lme_mTP; lme_mMZl = lme_mTP; 

%% Select year and member
yr = 1; %1:length(Y)
for mem = 1:length(M)
    iy = yr;
    im = M(mem);
    
    %% READ IN ENSEMBLE MEMBER
    % Pelagic temperature
    ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_TEMP_150m.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for n = nvars
        varname = netcdf.inqVar(ncid, n-1);
        eval([ varname ' = netcdf.getVar(ncid,n-1, [0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    netcdf.close(ncid);
    
    % Bottom temperature
    ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_TEMP_bottom.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:(nvars-1)
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    [ni,nj] = size(TLONG);
    for n = nvars
        varname = netcdf.inqVar(ncid, n-1);
        eval([ varname ' = netcdf.getVar(ncid,n-1, [0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    netcdf.close(ncid);
    
    %Bottom detritus
    ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_POC_FLUX_IN_bottom.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:(nvars-1)
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    [ni,nj] = size(TLONG);
    for n = nvars
        varname = netcdf.inqVar(ncid, n-1);
        eval([ varname ' = netcdf.getVar(ncid,n-1, [0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    netcdf.close(ncid);
    
    %spC
    ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_spC_150m.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:(nvars-1)
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    [ni,nj] = size(TLONG);
    for n = nvars
        varname = netcdf.inqVar(ncid, n-1);
        eval([ varname ' = netcdf.getVar(ncid,n-1, [0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    netcdf.close(ncid);
    
    %diatC
    ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_diatC_150m.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:(nvars-1)
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    for n = nvars
        varname = netcdf.inqVar(ncid, n-1);
        eval([ varname ' = netcdf.getVar(ncid,n-1,[0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    netcdf.close(ncid);
    
    %diazC
    ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_diazC_150m.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:(nvars-1)
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    for n = nvars
        varname = netcdf.inqVar(ncid, n-1);
        eval([ varname ' = netcdf.getVar(ncid,n-1,[0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    netcdf.close(ncid);
    
    % Zooplankton
    ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_zooC_150m.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:(nvars-1)
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    for n = nvars
        varname = netcdf.inqVar(ncid, n-1);
        %eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        %zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs(1)-1],[360 180 length(z200) length(runs)]);
        eval([ varname ' = netcdf.getVar(ncid,n-1,[0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    netcdf.close(ncid);
    
    % ZooLoss
    ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_zoo_loss_150m.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = 1:(nvars-1)
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    for n = nvars
        varname = netcdf.inqVar(ncid, n-1);
        %eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        %zmeso = netcdf.getVar(ncid,n-1, [0,0,0,runs(1)-1],[360 180 length(z200) length(runs)]);
        eval([ varname ' = netcdf.getVar(ncid,n-1,[0,0,0,im-1,yr-1],[ni,nj,length(L),1,1]);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    netcdf.close(ncid);
    
    %% Doubles and nans
    TEMP_150m = double(TEMP_150m);
    TEMP_bottom = double(TEMP_bottom);
    POC_FLUX_IN_bottom = double(POC_FLUX_IN_bottom);
    zooC_150m = double(zooC_150m);
    zoo_loss_150m = double(zoo_loss_150m);
    diatC_150m = double(diatC_150m);
    diazC_150m = double(diazC_150m);
    spC_150m = double(spC_150m);
    
    POC_FLUX_IN_bottom(POC_FLUX_IN_bottom >= 9.9e+36) = nan;
    zooC_150m(zooC_150m >= 9.9e+36) = nan;
    zoo_loss_150m(zoo_loss_150m >= 9.9e+36) = nan;
    
    %% CALC LARGE FRACTION OF ZOOP FROM ALL PHYTO
    fracL = diatC_150m ./ (diatC_150m + spC_150m + diazC_150m);
    LzooC_150m = fracL .* zooC_150m;
    Lzoo_loss_150m = fracL .* zoo_loss_150m;
    
    clear diatC_150m spC_150m zooC_150m zoo_loss_150m diazC_150m
    
    %% ADD TO FOSI CLIMATOL
    
    % ADD DPLE DRIFT-CORR ANOMALIES TO FOSI VALUES
    %fosi_Tp fosi_Tb fosi_POC fosi_zooC fosi_loss
    TEMP_150m = TEMP_150m + clim_Tp;
    TEMP_bottom = TEMP_bottom + clim_Tb;
    POC_FLUX_IN_bottom = POC_FLUX_IN_bottom + clim_POC;
    LzooC_150m = LzooC_150m + clim_zooC;
    Lzoo_loss_150m = Lzoo_loss_150m + clim_loss;
    
    POC_FLUX_IN_bottom(POC_FLUX_IN_bottom<0) = 0.0;
    LzooC_150m(LzooC_150m<0) = 0.0;
    Lzoo_loss_150m(Lzoo_loss_150m<0) = 0.0;
    
    %% Vectorize ocean grid cells
    TEMP_150m = reshape(TEMP_150m,ni*nj,nt);
    TEMP_bottom = reshape(TEMP_bottom,ni*nj,nt);
    POC_FLUX_IN_bottom = reshape(POC_FLUX_IN_bottom,ni*nj,nt);
    LzooC_150m = reshape(LzooC_150m,ni*nj,nt);
    Lzoo_loss_150m = reshape(Lzoo_loss_150m,ni*nj,nt);
    
    TEMP_150m = TEMP_150m(ID,:);
    TEMP_bottom = TEMP_bottom(ID,:);
    POC_FLUX_IN_bottom = POC_FLUX_IN_bottom(ID,:);
    LzooC_150m = LzooC_150m(ID,:);
    Lzoo_loss_150m = Lzoo_loss_150m(ID,:);
    
    %% Global time means
    tTP(mem,:) = mean(TEMP_150m,'omitnan');
    tTB(mem,:) = mean(TEMP_bottom,'omitnan');
    tDet(mem,:) = mean(POC_FLUX_IN_bottom,'omitnan');
    tMZ(mem,:) = mean(LzooC_150m,'omitnan');
    tMZl(mem,:) = mean(Lzoo_loss_150m,'omitnan');
    
    %% LME time means
    for L=1:66
        lid = find(lme_grid==L);
        if (~isempty(lid))
            %mean biomass
            lme_mTP(L,:,mem) = mean(TEMP_150m(lid,:),'omitnan');
            lme_mTB(L,:,mem) = mean(TEMP_bottom(lid,:),'omitnan');
            lme_mDet(L,:,mem) = mean(POC_FLUX_IN_bottom(lid,:),'omitnan');
            lme_mMZ(L,:,mem) = mean(LzooC_150m(lid,:),'omitnan');
            lme_mMZl(L,:,mem) =mean(Lzoo_loss_150m(lid,:),'omitnan');
        end
    end
    
    %%
    clear TEMP_150m TEMP_bottom POC_FLUX_IN_bottom time LzooC_150m Lzoo_loss_150m

end

%%
save([fpath 'Time_Means_DPLE_LME_drivers_Y',num2str(Y(yr)),'.mat'],...
    'lme_mTP','lme_mTB','lme_mDet','lme_mMZ','lme_mMZl',...
    'tTP','tTB','tDet','tMZ','tMZl','Y','M','L','iy','yr');
