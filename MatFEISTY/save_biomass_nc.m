function save_biomass_nc(fname, biomass)
    % biomass is ndays x nxs x 9
    nt      = size(biomass, 1);
    nx      = size(biomass, 2);
    ngroup  = size(biomass, 3);

    % explicitly define coordinate values
    time = 1:nt;
    X = 1:nx; % eventually pass this in as an argument
    group = ['Sf'; 'Sp'; 'Sd'; 'Mf'; 'Mp'; 'Md'; 'Lp'; 'Ld'; 'bp'];

    % delete output file if it already exists
    if (~isfolder('model_output/'))
        mkdir('model_output/')
    end

    full_name = ['./model_output/', fname]
    if isfile(full_name)
        delete(sprintf('%s',full_name));
    end

    % Create netcdf file / define dimensions / variables
    nccreate(full_name, 'time', 'Dimensions', {'time', nt})
    nccreate(full_name, 'X', 'Dimensions', {'X', nx})
    nccreate(full_name, 'group', 'Dimensions', {'nchar', 2, 'group', ngroup},...
             'Datatype', 'char')
    nccreate(full_name, 'biomass', 'Dimensions', {'group', ngroup, 'X', nx,'time', nt})

    % fill coordinate information
    ncwriteatt(full_name, 'time', 'long_name', 'time')
    ncwriteatt(full_name, 'time', 'standard_name', 'time')
    ncwriteatt(full_name, 'time', 'units', 'day')
    ncwriteatt(full_name, 'time', 'calendar', 'noleap')
    ncwriteatt(full_name, 'time', 'axis', 'T')
    ncwrite(full_name, 'time', time)

    ncwriteatt(full_name, 'X', 'long_name', 'longitude')
    ncwriteatt(full_name, 'X', 'standard_name', 'longitude')
    ncwriteatt(full_name, 'X', 'units', 'degrees_east')
    ncwriteatt(full_name, 'X', 'axis', 'X')
    ncwrite(full_name, 'X', X)

    ncwriteatt(full_name, 'group', 'long_name', 'prognostic classes')
    ncwrite(full_name, 'group', transpose(group))

    ncwrite(full_name, 'biomass', permute(biomass, [3, 2, 1]))
end