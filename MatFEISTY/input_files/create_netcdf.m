function create_netcdf()
    mat_file = './feisty_input_climatol_daily_locs3.mat';
    nc_file = './feisty_input_climatol_daily_locs3.nc';
    delete(sprintf('%s',nc_file))

    % 1. Create netcdf file; dims are time = 365, loc = 3
    time=365; loc=3;
    %    (a) time variable
    nccreate(nc_file, 'time', 'Dimensions', {'time', time})

    %    (b) variables with just loc dimensions
    varnames = {'loc', 'lat', 'dep'};
    for id = 1:length(varnames)
        nccreate(nc_file, varnames{id},...
                'Dimensions', {'loc', loc})
    end

    %    (c) variables with time and loc dimensions
    varnames = {'Tp', 'Tb', 'POC', 'zoo', 'zoo_mort'};
    for id = 1:length(varnames)
        nccreate(nc_file, varnames{id},...
                'Dimensions', {'loc', loc, 'time', time})
    end

    % 2. Read .mat data, replace some negative values with 0
    load(mat_file,'GRD','ESM');

	% double time(time) ;
	% 	time:long_name = "time" ;
	% 	time:standard_name = "time" ;
	% 	time:units = "year" ;
	% 	time:calendar = "365_day" ;
	% 	time:axis = "T" ;
    ncwriteatt(nc_file, 'time', 'long_name', 'time')
    ncwriteatt(nc_file, 'time', 'standard_name', 'time')
    ncwriteatt(nc_file, 'time', 'units', 'day')
    ncwriteatt(nc_file, 'time', 'calendar', 'noleap')
    ncwriteatt(nc_file, 'time', 'axis', 'T')
    ncwrite(nc_file, 'time', 0:364)

	% double loc(loc) ;
	% 	loc:long_name = "longitude" ;
	% 	loc:standard_name = "longitude" ;
	% 	loc:units = "degrees_east" ;
	% 	loc:axis = "X" ;
    ncwriteatt(nc_file, 'loc', 'long_name', 'longitude')
    ncwriteatt(nc_file, 'loc', 'standard_name', 'longitude')
    ncwriteatt(nc_file, 'loc', 'units', 'degrees_east')
    ncwriteatt(nc_file, 'loc', 'axis', 'X')
    ncwrite(nc_file, 'loc', GRD.LON)

	% double lat(loc) ;
	% 	lat:long_name = "latitude" ;
	% 	lat:standard_name = "latitude" ;
	% 	lat:units = "degrees_north" ;
	% 	lat:axis = "Y" ;
    ncwriteatt(nc_file, 'lat', 'long_name', 'latitude')
    ncwriteatt(nc_file, 'lat', 'standard_name', 'latitude')
    ncwriteatt(nc_file, 'lat', 'units', 'degrees_north')
    ncwriteatt(nc_file, 'lat', 'axis', 'Y')
    ncwrite(nc_file, 'lat', GRD.LAT)

	% double dep(loc) ;
	% 	dep:long_name = "bottom depth" ;
	% 	dep:standard_name = "bottom depth" ;
	% 	dep:units = "m" ;
	% 	dep:axis = "Z" ;
    ncwriteatt(nc_file, 'dep', 'long_name', 'bottom depth')
    ncwriteatt(nc_file, 'dep', 'standard_name', 'bottom depth')
    ncwriteatt(nc_file, 'dep', 'units', 'm')
    ncwriteatt(nc_file, 'dep', 'axis', 'Z')
    ncwrite(nc_file, 'dep', GRD.Z)

	% float Tp(time, loc) ;
	% 	Tp:long_name = "Pelagic mean temperature" ;
	% 	Tp:units = "degrees C" ;
    ncwriteatt(nc_file, 'Tp', 'long_name', 'Pelagic mean temperature')
    ncwriteatt(nc_file, 'Tp', 'units', 'degrees C')
    ncwrite(nc_file, 'Tp', ESM.Tp)

	% float Tb(time, loc) ;
	% 	Tb:long_name = "Bottom temperature" ;
	% 	Tb:units = "degrees C" ;
    ncwriteatt(nc_file, 'Tb', 'long_name', 'Bottom temperature')
    ncwriteatt(nc_file, 'Tb', 'units', 'degrees C')
    ncwrite(nc_file, 'Tb', ESM.Tb)

	% float POC(time, loc) ;
	% 	POC:long_name = "Particulate organic matter flux to seafloor" ;
	% 	POC:units = "g m-2 d-1" ;
    ncwriteatt(nc_file, 'POC', 'long_name', 'Particulate organic matter flux to seafloor')
    ncwriteatt(nc_file, 'POC', 'units', 'g m-2 d-1')
    ncwrite(nc_file, 'POC', ESM.det)

	% float zoo(time, loc) ;
	% 	zoo:long_name = "Biomass of mesozooplankton" ;
	% 	zoo:units = "g m-2" ;
    ncwriteatt(nc_file, 'zoo', 'long_name', 'Biomass of mesozooplankton')
    ncwriteatt(nc_file, 'zoo', 'units', 'g m-2')
    ncwrite(nc_file, 'zoo', ESM.Zm)

	% float zoo_mort(time, loc) ;
	% 	zoo_mort:long_name = "Mortality loss of mesozooplankton" ;
	% 	zoo_mort:units = "g m-2 d-1" ;
    ncwriteatt(nc_file, 'zoo_mort', 'long_name', 'Mortality loss of mesozooplankton')
    ncwriteatt(nc_file, 'zoo_mort', 'units', 'g m-2 d-1')
    ncwrite(nc_file, 'zoo_mort', ESM.dZm)
end