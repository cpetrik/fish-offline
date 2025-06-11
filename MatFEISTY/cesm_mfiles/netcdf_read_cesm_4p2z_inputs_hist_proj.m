% Read CESM 4P2Z netcdf

clear
close all

fpath = '/project/Feisty/GCM_Data/CESM/4P2Z/Hist/';

%% All inputs
% ncdisp([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.mesozooC_zint_150m.18500101-19000101.nc'])
% ncdisp([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.mesozoo_loss_zint_150m.18500101-19000101.nc'])
% ncdisp([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.pocToFloor_2.18500101-19000101.nc'])
% ncdisp([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.TEMP_BOTTOM_2.18500101-19000101.nc'])
% ncdisp([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.TEMP_mean_150m.18500101-19000101.nc'])

%%
% mesozooC_zint_150m
% Size:       320x384x18250
% Dimensions: nlon,nlat,time
% Datatype:   single
FillValue    = 9.969209968386869e+36;
mesozooC_long_name     = 'Mesozooplankton Carbon 0-150m Vertical Integral';
mesozooC_units         = 'mmol/m^3 cm  == nmolC cm-2';
missing_value = 9.969209968386869e+36;

% mesozoo_loss_zint_150m
mesozoo_loss_long_name     = 'Mesozooplankton Loss Vertical Integral, 0-150m';
mesozoo_loss_units         = 'mmol/m^3 cm/s  == nmolC cm-2 s-1';

% pocToFloor_2
poc_long_name     = 'POC Flux Hitting Sea Floor';
poc_units         = 'nmol/cm^2/s  == nmolC cm-2 s-1';

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

%% Water cells only

ncid = netcdf.open([fpath 'b.e21p4.BHIST.f09_g17.4p2z.002.pop.h.ecosys.nday1.TEMP_mean_150m.',num2str(syr(1)),'-',num2str(eyr(1)),'.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for k=18
    varname = netcdf.inqVar(ncid, k-1);
    eval([ varname ' = netcdf.getVar(ncid,k-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
netcdf.close(ncid);

WID = find(HT(:)>0);

%%

for j=1:length(syr)
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

    %% reshape, double, nans, zeros, water cells only
    [ni,nj,nt] = size(pocToFloor_2);

    TEMP_mean_150m = reshape(TEMP_mean_150m,ni*nj,nt);
    TEMP_BOTTOM_2 = reshape(TEMP_BOTTOM_2,ni*nj,nt);
    pocToFloor_2 = reshape(pocToFloor_2,ni*nj,nt);
    mesozooC_zint_150m = reshape(mesozooC_zint_150m,ni*nj,nt);
    mesozoo_loss_zint_150m = reshape(mesozoo_loss_zint_150m,ni*nj,nt);

    Tp = double(TEMP_mean_150m(WID,:));
    Tb = double(TEMP_BOTTOM_2(WID,:));
    Det = double(pocToFloor_2(WID,:));
    Zmeso = double(mesozooC_zint_150m(WID,:));
    ZmLoss = double(mesozoo_loss_zint_150m(WID,:));

    Tp(Tp >= 9.9e+36) = nan;
    Tb(Tb >= 9.9e+36) = nan;
    Det(Det >= 9.9e+36) = nan;
    Zmeso(Zmeso >= 9.9e+36) = nan;
    ZmLoss(ZmLoss >= 9.9e+36) = nan;

    Zmeso(Zmeso<0) = 0.0;
    ZmLoss(ZmLoss<0) = 0.0;
    Det(Det<0) = 0.0;

    %% units

    % meso zoo: nmolC cm-2 to g(WW) m-2
    % 1e9 nmol in 1 mol C
    % 1e4 cm2 in 1 m2
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    Zmeso = Zmeso * 1e-9 * 1e4 * 12.01 * 9.0;

    % medium zoo mortality: nmolC cm-2 s-1 to g(WW) m-2 d-1
    % 1e9 nmol in 1 mol C
    % 1e4 cm2 in 1 m2
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    ZmLoss = ZmLoss * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

    % detrital flux to benthos: nmolC cm-2 s-1 to g(WW) m-2 d-1
    % 1e9 nmol in 1 mol C
    % 1e4 cm2 in 1 m2
    % 12.01 g C in 1 mol C
    % 1 g dry W in 9 g wet W (Pauly & Christiansen)
    Det = Det * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;

    %% Time
    tyr = (time)/365;
    nyrs = nt/365;

    dst = 1:365:ndays(j);
    den = 365:365:ndays(j);

    for y = 1:nyrs
        yr = floor(tyr(dst(y)));

        range = dst(y):den(y);

        ESM.Tp = Tp(:,range);
        ESM.Tb = Tb(:,range);
        ESM.Zm = Zmeso(:,range);
        ESM.dZm = ZmLoss(:,range);
        ESM.det = Det(:,range);

        % save
        save([fpath 'Data_cesm_4p2z_daily_hist_',num2str(yr),'.mat'], 'ESM');
    end

end







