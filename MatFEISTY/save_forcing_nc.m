function save_forcing_nc(fname, GRD, ESM)
    % delete output file if it already exists
    if (~isfolder('model_output/'))
        mkdir('model_output/')
    end

    full_name = ['./model_output/', fname]
    if isfile(full_name)
        delete(sprintf('%s',full_name));
    end

    % 1. Create netcdf file; dims come from ESM
    nx = size(ESM.Tp, 1);
    nt = size(ESM.Tp, 2);
    %    (a) time coordinate
    nccreate(full_name, 'time', 'Dimensions', {'time', nt})

    %    (b) variables with just X dimensions (including X coordinate)
    if isfield(GRD, 'LAT')
        varnames = {'X', 'lat', 'dep'};
    else
        varnames = {'X', 'dep'};
    end
    for id = 1:length(varnames)
        nccreate(full_name, varnames{id},...
                'Dimensions', {'X', nx})
    end

    %    (c) zooplankton coordinate
    nccreate(full_name, 'zooplankton', 'Dimension', {'nchar', 3, 'zooplankton', 1},...
    'Datatype', 'char')

    %    (d) variables with time and X dimensions
    varnames = {'T_pelagic', 'T_bottom', 'poc_flux_bottom'};
    for id = 1:length(varnames)
        nccreate(full_name, varnames{id},...
                'Dimensions', {'X', nx, 'time', nt, 'zooplankton', 1})
    end

    %    (e) variables with zooplankton, time, and X dimensions
    varnames = {'zooC', 'zoo_mort'};
    for id = 1:length(varnames)
        nccreate(full_name, varnames{id},...
                'Dimensions', {'X', nx, 'time', nt, 'zooplankton', 1})
    end

	% double time(time) ;
	% 	time:long_name = "time" ;
	% 	time:standard_name = "time" ;
	% 	time:units = "year" ;
	% 	time:calendar = "365_day" ;
	% 	time:axis = "T" ;
    ncwriteatt(full_name, 'time', 'long_name', 'time')
    ncwriteatt(full_name, 'time', 'standard_name', 'time')
    ncwriteatt(full_name, 'time', 'units', 'day')
    ncwriteatt(full_name, 'time', 'calendar', 'noleap')
    ncwriteatt(full_name, 'time', 'axis', 'T')
    ncwrite(full_name, 'time', 0:364)

	% double X(X) ;
	% 	X:long_name = "longitude" ;
	% 	X:standard_name = "longitude" ;
	% 	X:units = "degrees_east" ;
	% 	X:axis = "X" ;
    ncwriteatt(full_name, 'X', 'long_name', 'longitude')
    ncwriteatt(full_name, 'X', 'standard_name', 'longitude')
    ncwriteatt(full_name, 'X', 'units', 'degrees_east')
    ncwriteatt(full_name, 'X', 'axis', 'X')
    if isfield(GRD, 'LON')
        ncwrite(full_name, 'X', GRD.LON)
    else
        ncwrite(full_name, 'X', 1:nx)
    end

    % zooplankton
    ncwrite(full_name, 'zooplankton', transpose(['Zoo']))

	% double lat(X) ;
	% 	lat:long_name = "latitude" ;
	% 	lat:standard_name = "latitude" ;
	% 	lat:units = "degrees_north" ;
	% 	lat:axis = "Y" ;
    if isfield(GRD, 'LAT')
        ncwriteatt(full_name, 'lat', 'long_name', 'latitude')
        ncwriteatt(full_name, 'lat', 'standard_name', 'latitude')
        ncwriteatt(full_name, 'lat', 'units', 'degrees_north')
        ncwriteatt(full_name, 'lat', 'axis', 'Y')
        ncwrite(full_name, 'lat', GRD.LAT)
    end

	% double dep(X) ;
	% 	dep:long_name = "bottom depth" ;
	% 	dep:standard_name = "bottom depth" ;
	% 	dep:units = "m" ;
	% 	dep:axis = "Z" ;
    ncwriteatt(full_name, 'dep', 'long_name', 'bottom depth')
    ncwriteatt(full_name, 'dep', 'standard_name', 'bottom depth')
    ncwriteatt(full_name, 'dep', 'units', 'm')
    ncwriteatt(full_name, 'dep', 'axis', 'Z')
    ncwrite(full_name, 'dep', GRD.Z)

	% float T_pelagic(time, X) ;
	% 	T_pelagic:long_name = "Pelagic mean temperature" ;
	% 	T_pelagic:units = "degrees C" ;
    ncwriteatt(full_name, 'T_pelagic', 'long_name', 'Pelagic mean temperature')
    ncwriteatt(full_name, 'T_pelagic', 'units', 'degrees C')
    ncwrite(full_name, 'T_pelagic', ESM.Tp)

	% float T_bottom(time, X) ;
	% 	T_bottom:long_name = "Bottom temperature" ;
	% 	T_bottom:units = "degrees C" ;
    ncwriteatt(full_name, 'T_bottom', 'long_name', 'Bottom temperature')
    ncwriteatt(full_name, 'T_bottom', 'units', 'degrees C')
    ncwrite(full_name, 'T_bottom', ESM.Tb)

	% float poc_flux_bottom(time, X) ;
	% 	poc_flux_bottom:long_name = "Particulate organic matter flux to seafloor" ;
	% 	poc_flux_bottom:units = "g m-2 d-1" ;
    ncwriteatt(full_name, 'poc_flux_bottom', 'long_name', 'Particulate organic matter flux to seafloor')
    ncwriteatt(full_name, 'poc_flux_bottom', 'units', 'g m-2 d-1')
    ncwrite(full_name, 'poc_flux_bottom', ESM.det)

	% float zooC(time, X) ;
	% 	zooC:long_name = "Biomass of mesozooplankton" ;
	% 	zooC:units = "g m-2" ;
    ncwriteatt(full_name, 'zooC', 'long_name', 'Biomass of mesozooplankton')
    ncwriteatt(full_name, 'zooC', 'units', 'g m-2')
    ncwrite(full_name, 'zooC', ESM.Zm)

	% float zoo_mort(time, X) ;
	% 	zoo_mort:long_name = "Mortality loss of mesozooplankton" ;
	% 	zoo_mort:units = "g m-2 d-1" ;
    ncwriteatt(full_name, 'zoo_mort', 'long_name', 'Mortality loss of mesozooplankton')
    ncwriteatt(full_name, 'zoo_mort', 'units', 'g m-2 d-1')
    ncwrite(full_name, 'zoo_mort', ESM.dZm)
end