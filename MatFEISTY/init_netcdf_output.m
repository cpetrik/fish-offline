function init_netcdf_output(fname, nt, nx, ngroup, GRD)
    % (1) delete output file if it already exists
    if (~isfolder('model_output/'))
        mkdir('model_output/')
    end

    full_name = ['./model_output/', fname];
    if isfile(full_name)
        delete(sprintf('%s',full_name));
    end

    % (2) Create netcdf file, define dimensions and variables
    %     (a) time
    nccreate(full_name, 'time', 'Dimensions', {'time', nt})

    %     (b) variables with just X dimensions (including X coordinate)
    if isfield(GRD, 'LAT')
        varnames = {'X', 'lat', 'dep'};
    else
        varnames = {'X', 'dep'};
    end
    for id = 1:length(varnames)
        nccreate(full_name, varnames{id}, 'Dimensions', {'X', nx})
    end

    %     (c) variables with additional dimensions
    nccreate(full_name, 'zooplankton', 'Dimension', {'char3', 3, 'zooplankton', 1},...
    'Datatype', 'char')
    nccreate(full_name, 'group', 'Dimensions', {'char2', 2, 'group', ngroup},...
             'Datatype', 'char')
    nccreate(full_name, 'biomass', 'Dimensions', {'group', ngroup, 'X', nx,'time', nt})
    varnames = {'T_pelagic', 'T_bottom', 'poc_flux_bottom'};
    for id = 1:length(varnames)
        nccreate(full_name, varnames{id},...
                'Dimensions', {'X', nx, 'time', nt, 'zooplankton', 1})
    end
    varnames = {'zooC', 'zoo_mort'};
    for id = 1:length(varnames)
        nccreate(full_name, varnames{id},...
                'Dimensions', {'X', nx, 'time', nt, 'zooplankton', 1})
    end

    %     (d) define metadata
    ncwriteatt(full_name, 'time', 'long_name', 'time')
    ncwriteatt(full_name, 'time', 'standard_name', 'time')
    ncwriteatt(full_name, 'time', 'units', 'days since 0001-01-01 00:00:00')
    ncwriteatt(full_name, 'time', 'calendar', 'noleap')
    ncwriteatt(full_name, 'time', 'axis', 'T')

    ncwriteatt(full_name, 'X', 'long_name', 'longitude')
    ncwriteatt(full_name, 'X', 'standard_name', 'longitude')
    ncwriteatt(full_name, 'X', 'units', 'degrees_east')
    ncwriteatt(full_name, 'X', 'axis', 'X')

    if isfield(GRD, 'LAT')
        ncwriteatt(full_name, 'lat', 'long_name', 'latitude')
        ncwriteatt(full_name, 'lat', 'standard_name', 'latitude')
        ncwriteatt(full_name, 'lat', 'units', 'degrees_north')
        ncwriteatt(full_name, 'lat', 'axis', 'Y')
    end

    ncwriteatt(full_name, 'dep', 'long_name', 'bottom depth')
    ncwriteatt(full_name, 'dep', 'standard_name', 'bottom depth')
    ncwriteatt(full_name, 'dep', 'units', 'm')
    ncwriteatt(full_name, 'dep', 'axis', 'Z')

    ncwriteatt(full_name, 'group', 'long_name', 'prognostic classes')

    ncwriteatt(full_name, 'T_pelagic', 'long_name', 'Pelagic mean temperature')
    ncwriteatt(full_name, 'T_pelagic', 'units', 'degrees C')

    ncwriteatt(full_name, 'T_bottom', 'long_name', 'Bottom temperature')
    ncwriteatt(full_name, 'T_bottom', 'units', 'degrees C')

    ncwriteatt(full_name, 'poc_flux_bottom', 'long_name', 'Particulate organic matter flux to seafloor')
    ncwriteatt(full_name, 'poc_flux_bottom', 'units', 'g m-2 d-1')

    ncwriteatt(full_name, 'zooC', 'long_name', 'Biomass of mesozooplankton')
    ncwriteatt(full_name, 'zooC', 'units', 'g m-2')

    ncwriteatt(full_name, 'zoo_mort', 'long_name', 'Mortality loss of mesozooplankton')
    ncwriteatt(full_name, 'zoo_mort', 'units', 'g m-2 d-1')

end