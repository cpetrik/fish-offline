% Read CESM 4P2Z netcdf

clear
close all

fpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/4P2Z/Hist/';

%% All inputs
ncdisp([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.mesozooC_zint_150m.18500101-19000101.nc'])
ncdisp([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.mesozoo_loss_zint_150m.18500101-19000101.nc'])
ncdisp([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.pocToFloor_2.18500101-19000101.nc'])
ncdisp([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.TEMP_BOTTOM_2.18500101-19000101.nc'])
ncdisp([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.TEMP_mean_150m.18500101-19000101.nc'])

%%
% mesozooC_zint_150m
% Size:       320x384x18250
% Dimensions: nlon,nlat,time
% Datatype:   single
FillValue    = 9.969209968386869e+36;
mesozooC_long_name     = 'Mesozooplankton Carbon 0-150m Vertical Integral';
mesozooC_units         = 'mmol/m^3 cm';
missing_value = 9.969209968386869e+36;

% mesozoo_loss_zint_150m
mesozoo_loss_long_name     = 'Mesozooplankton Loss Vertical Integral, 0-150m';
mesozoo_loss_units         = 'mmol/m^3 cm/s';
 
% pocToFloor_2
poc_long_name     = 'POC Flux Hitting Sea Floor';
poc_units         = 'nmol/cm^2/s';

% TEMP_BOTTOM_2
TEMP_BOTTOM_long_name     = 'Potential temperature Value at Sea Floor';
TEMP_BOTTOM_units         = 'degC';

% TEMP_mean_150m
TEMP_long_name     = 'Potential temperature 0-150m Vertical Mean';
TEMP_units         = 'degC cm';

time_units = 'days since 0000-01-01 00:00:00';
calendar   = 'noleap';


%%
syr = [18500101;19000102;19500102;20000102];
eyr = [19000101;19500101;20000101;20141231];

ndays = [50*365; (50*365); (50*365); (15*365)];
nde = cumsum(ndays);
nds = [1;nde(1:3)+1];

%% Water cells only

ncid = netcdf.open([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.TEMP_mean_150m.',num2str(syr(j)),'-',num2str(eyr(j)),'.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for k=18
    varname = netcdf.inqVar(ncid, k-1);
    eval([ varname ' = netcdf.getVar(ncid,k-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
netcdf.close(ncid);

[ni,nj] = size(HT);

WID = find(HT(:)>0);

%%
poc = nan*ones(length(WID),nde(end));

%%

for j=1%:length(syrs)
    ncid = netcdf.open([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.TEMP_mean_150m.',num2str(syr(j)),'-',num2str(eyr(j)),'.nc'],'NC_NOWRITE');
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
    for i = nvars %1:(nvars)
        varname = netcdf.inqVar(ncid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    netcdf.close(ncid);


    bcid = netcdf.open([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.TEMP_BOTTOM_2.',num2str(syr(j)),'-',num2str(eyr(j)),'.nc'],'NC_NOWRITE');
    [~,bvars,~,~] = netcdf.inq(bcid);
    for i = bvars 
        varname = netcdf.inqVar(bcid, i-1);
        eval([ varname ' = netcdf.getVar(bcid,i-1);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    netcdf.close(bcid);


    zcid = netcdf.open([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.mesozooC_zint_150m.',num2str(syr(j)),'-',num2str(eyr(j)),'.nc'],'NC_NOWRITE');
    [~,zvars,~,~] = netcdf.inq(zcid);
    for i = zvars 
        varname = netcdf.inqVar(zcid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    netcdf.close(zcid);


    lcid = netcdf.open([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.mesozoo_loss_zint_150m.',num2str(syr(j)),'-',num2str(eyr(j)),'.nc'],'NC_NOWRITE');
    [~,lvars,~,~] = netcdf.inq(lcid);
    for i = lvars 
        varname = netcdf.inqVar(lcid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    netcdf.close(lcid);


    pcid = netcdf.open([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.pocToFloor_2.',num2str(syr(j)),'-',num2str(eyr(j)),'.nc'],'NC_NOWRITE');
    [~,pvars,~,~] = netcdf.inq(pcid);
    for i = 1:(pvars)
        varname = netcdf.inqVar(pcid, i-1);
        eval([ varname ' = netcdf.getVar(ncid,i-1);']);
        eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
    end
    netcdf.close(pcid);



    poc(:,:,nds(j):nde(j)) = pocToFloor_2;

end

%% Time
%yr = (time-time(1)+1)/365;
%yr = (time-time(1))/365;
yr = (time)/365;

%%
save([fpath 'g.e22a06.G1850ECOIAF_JRA_PHYS_DEV.TL319_g17.4p4z.004.FIESTY-forcing.mat']);





