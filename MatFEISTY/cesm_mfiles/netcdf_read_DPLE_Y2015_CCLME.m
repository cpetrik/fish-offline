% FEISTY output at CCLME locations

clear all
close all

%% Map data
cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;
glme = double(lme_mask);
glme(glme<0) = nan;
tlme = glme(ID);

%%
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];

%pick year
StartYr = 2015;
%loop over members
submem = 1:40;
%%
for mem=1:length(submem) %will loop over
    Member = submem(mem);
    harv = ['v14_Y' num2str(StartYr) '_M' num2str(Member) '_All_fish03_' ];
    
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
    [nid,nt] = size(biomass);
    
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
    
    %% Netcdf OUTPUTS =================================================
    t=time;
    mo=t/12;
    mo=mo+StartYr;
    
    allF = SF.bio + MF.bio;
    allP = SP.bio + MP.bio + LP.bio;
    allD = SD.bio + MD.bio + LD.bio;
    allB = Bent.bio;
    allS = SF.bio + SP.bio + SD.bio;
    allM = MF.bio + MP.bio + MD.bio;
    allL = LP.bio + LD.bio;
    
    %% Reshape to lat,lon,yr
    % Extract CC LME lat-lon only
    vid = find(tlme==3);
    
    AllF = NaN*ones(ni,nj,nt);
    AllP = NaN*ones(ni,nj,nt);
    AllD = NaN*ones(ni,nj,nt);
    AllB = NaN*ones(ni,nj,nt);
    AllS = NaN*ones(ni,nj,nt);
    AllM = NaN*ones(ni,nj,nt);
    AllL = NaN*ones(ni,nj,nt);
    
    for z=1:nt
        Zf=NaN*ones(ni,nj);
        Zp=NaN*ones(ni,nj);
        Zd=NaN*ones(ni,nj);
        Zs=NaN*ones(ni,nj);
        Zm=NaN*ones(ni,nj);
        Zl=NaN*ones(ni,nj);
        Zb=NaN*ones(ni,nj);
        
        Zf(GRD.ID(vid))=allF(vid,z);
        Zp(GRD.ID(vid))=allP(vid,z);
        Zd(GRD.ID(vid))=allD(vid,z);
        Zs(GRD.ID(vid))=allS(vid,z);
        Zm(GRD.ID(vid))=allM(vid,z);
        Zl(GRD.ID(vid))=allL(vid,z);
        Zb(GRD.ID(vid))=allB(vid,z);
        
        AllF(:,:,z) = Zf;
        AllP(:,:,z) = Zp;
        AllD(:,:,z) = Zd;
        AllB(:,:,z) = Zb;
        AllS(:,:,z) = Zs;
        AllM(:,:,z) = Zm;
        AllL(:,:,z) = Zl;
    end
    
    %% Extract CC LME lat-lon only to reduce size
    test = AllF(:,:,1);
    cccol = find(~isnan(nansum(test)));
    ccrow = find(~isnan(nansum(test,2)));
    
    AllF = AllF(ccrow,cccol,:);
    AllP = AllP(ccrow,cccol,:);
    AllD = AllD(ccrow,cccol,:);
    AllB = AllB(ccrow,cccol,:);
    AllS = AllS(ccrow,cccol,:);
    AllM = AllM(ccrow,cccol,:);
    AllL = AllL(ccrow,cccol,:);
    
    lat2 = TLAT(ccrow,cccol);
    lon2 = TLONG(ccrow,cccol);
    
    %%
    
    fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
    save([fpath 'CESM_DPLE_CCLME_outputs_monthly_' harv cfile '.mat'],...
        'time','mo','lat2','lon2','cccol','ccrow',...
        'AllF','AllD','AllP','AllB',...
        'AllS','AllM','AllL');
    
end

%% test
plotminlat=21; %Set these bounds for your data
plotmaxlat=50;
plotminlon=-130;
plotmaxlon=-110;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%%
test = AllD(:,:,1);

figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(lat2,lon2,test)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
