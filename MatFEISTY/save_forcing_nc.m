function save_forcing_nc(fname, GRD, ESM)
    % delete output file if it already exists
    if isfile(fname)
        delete(sprintf('%s',fname));
    end

    % 1. Create netcdf file; dims come from ESM
    nx = size(ESM.Tp, 1);
    nt = size(ESM.Tp, 2);
    %    (a) time coordinate
    nccreate(fname, 'time', 'Dimensions', {'time', nt})

    %    (b) variables with just X dimensions (including X coordinate)
    if isfield(GRD, 'LAT')
        varnames = {'X', 'lat', 'dep'};
    else
        varnames = {'X', 'dep'};
    end
    for id = 1:length(varnames)
        nccreate(fname, varnames{id},...
                'Dimensions', {'X', nx})
    end

    %    (c) zooplankton coordinate
    nccreate(fname, 'zooplankton', 'Dimension', {'nchar', 3, 'zooplankton', 1},...
    'Datatype', 'char')

    %    (d) variables with time and X dimensions
    varnames = {'T_pelagic', 'T_bottom', 'poc_flux_bottom'};
    for id = 1:length(varnames)
        nccreate(fname, varnames{id},...
                'Dimensions', {'X', nx, 'time', nt, 'zooplankton', 1})
    end

    %    (e) variables with zooplankton, time, and X dimensions
    varnames = {'zooC', 'zoo_mort'};
    for id = 1:length(varnames)
        nccreate(fname, varnames{id},...
                'Dimensions', {'X', nx, 'time', nt, 'zooplankton', 1})
    end

	% double time(time) ;
	% 	time:long_name = "time" ;
	% 	time:standard_name = "time" ;
	% 	time:units = "year" ;
	% 	time:calendar = "365_day" ;
	% 	time:axis = "T" ;
    ncwriteatt(fname, 'time', 'long_name', 'time')
    ncwriteatt(fname, 'time', 'standard_name', 'time')
    ncwriteatt(fname, 'time', 'units', 'day')
    ncwriteatt(fname, 'time', 'calendar', 'noleap')
    ncwriteatt(fname, 'time', 'axis', 'T')
    ncwrite(fname, 'time', 0:364)

	% double X(X) ;
	% 	X:long_name = "longitude" ;
	% 	X:standard_name = "longitude" ;
	% 	X:units = "degrees_east" ;
	% 	X:axis = "X" ;
    ncwriteatt(fname, 'X', 'long_name', 'longitude')
    ncwriteatt(fname, 'X', 'standard_name', 'longitude')
    ncwriteatt(fname, 'X', 'units', 'degrees_east')
    ncwriteatt(fname, 'X', 'axis', 'X')
    if isfield(GRD, 'LON')
        ncwrite(fname, 'X', GRD.LON)
    else
        ncwrite(fname, 'X', 1:nx)
    end

    % zooplankton
    ncwrite(fname, 'zooplankton', transpose(['Zoo']))

	% double lat(X) ;
	% 	lat:long_name = "latitude" ;
	% 	lat:standard_name = "latitude" ;
	% 	lat:units = "degrees_north" ;
	% 	lat:axis = "Y" ;
    if isfield(GRD, 'LAT')
        ncwriteatt(fname, 'lat', 'long_name', 'latitude')
        ncwriteatt(fname, 'lat', 'standard_name', 'latitude')
        ncwriteatt(fname, 'lat', 'units', 'degrees_north')
        ncwriteatt(fname, 'lat', 'axis', 'Y')
        ncwrite(fname, 'lat', GRD.LAT)
    end

	% double dep(X) ;
	% 	dep:long_name = "bottom depth" ;
	% 	dep:standard_name = "bottom depth" ;
	% 	dep:units = "m" ;
	% 	dep:axis = "Z" ;
    ncwriteatt(fname, 'dep', 'long_name', 'bottom depth')
    ncwriteatt(fname, 'dep', 'standard_name', 'bottom depth')
    ncwriteatt(fname, 'dep', 'units', 'm')
    ncwriteatt(fname, 'dep', 'axis', 'Z')
    ncwrite(fname, 'dep', GRD.Z)

	% float T_pelagic(time, X) ;
	% 	T_pelagic:long_name = "Pelagic mean temperature" ;
	% 	T_pelagic:units = "degrees C" ;
    ncwriteatt(fname, 'T_pelagic', 'long_name', 'Pelagic mean temperature')
    ncwriteatt(fname, 'T_pelagic', 'units', 'degrees C')
    ncwrite(fname, 'T_pelagic', ESM.Tp)

	% float T_bottom(time, X) ;
	% 	T_bottom:long_name = "Bottom temperature" ;
	% 	T_bottom:units = "degrees C" ;
    ncwriteatt(fname, 'T_bottom', 'long_name', 'Bottom temperature')
    ncwriteatt(fname, 'T_bottom', 'units', 'degrees C')
    ncwrite(fname, 'T_bottom', ESM.Tb)

	% float poc_flux_bottom(time, X) ;
	% 	poc_flux_bottom:long_name = "Particulate organic matter flux to seafloor" ;
	% 	poc_flux_bottom:units = "g m-2 d-1" ;
    ncwriteatt(fname, 'poc_flux_bottom', 'long_name', 'Particulate organic matter flux to seafloor')
    ncwriteatt(fname, 'poc_flux_bottom', 'units', 'g m-2 d-1')
    ncwrite(fname, 'poc_flux_bottom', ESM.det)

	% float zooC(time, X) ;
	% 	zooC:long_name = "Biomass of mesozooplankton" ;
	% 	zooC:units = "g m-2" ;
    ncwriteatt(fname, 'zooC', 'long_name', 'Biomass of mesozooplankton')
    ncwriteatt(fname, 'zooC', 'units', 'g m-2')
    ncwrite(fname, 'zooC', ESM.Zm)

	% float zoo_mort(time, X) ;
	% 	zoo_mort:long_name = "Mortality loss of mesozooplankton" ;
	% 	zoo_mort:units = "g m-2 d-1" ;
    ncwriteatt(fname, 'zoo_mort', 'long_name', 'Mortality loss of mesozooplankton')
    ncwriteatt(fname, 'zoo_mort', 'units', 'g m-2 d-1')
    ncwrite(fname, 'zoo_mort', ESM.dZm)
end