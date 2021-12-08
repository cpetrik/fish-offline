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
for mem=1%:length(submem) %will loop over
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
    
    %% Extract CC LME lat-lon only
    gid = find(glme==3);
    vid = find(tlme==3);
    
    Zf=NaN*ones(ni,nj,nt);
    Zp=NaN*ones(ni,nj,nt);
    Zd=NaN*ones(ni,nj,nt);
    Zs=NaN*ones(ni,nj,nt);
    Zm=NaN*ones(ni,nj,nt);
    Zl=NaN*ones(ni,nj,nt);
    Zb=NaN*ones(ni,nj,nt);
    
    Zf(GRD.ID(vid))=allF(vid,:);
    Zp(GRD.ID)=sp_sbio;
    Zd(GRD.ID)=sd_sbio;
    Zs(GRD.ID)=mf_sbio;
    Zm(GRD.ID)=mp_sbio;
    Zl(GRD.ID)=md_sbio;
    Zld(GRD.ID)=ld_sbio;
    Zb(GRD.ID)=b_sbio;
    %%
    
    fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
    save([fpath 'CESM_DPLE_outputs_monthly_' harv cfile '.mat'],'time','mo',...
        'allF','allD','allP','allB','yr20');
    
end


%% test
plotminlat=22; %Set these bounds for your data
plotmaxlat=50;
plotminlon=-130;
plotmaxlon=-110;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;
%%
ZF=NaN*ones(ni,nj);
ZF(GRD.ID)=allF(:,1);

figure(3)
axesm ('Robinson','MapLatLimit',latlim,'MapLonLimit',lonlim,'frame','on',...
    'Grid','off','FLineWidth',1)
surfm(TLAT,TLONG,glme)
h=patchm(coastlat+0.5,coastlon+0.5,'w','FaceColor',[0.75 0.75 0.75]);
caxis([1 10]);
hcb = colorbar('h');
set(gcf,'renderer','painters')
