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
    if isfile(fname)
        delete(sprintf('%s',fname));
    end

    % Create netcdf file / define dimensions / variables
    nccreate(fname, 'time', 'Dimensions', {'time', nt})
    nccreate(fname, 'X', 'Dimensions', {'X', nx})
    nccreate(fname, 'group', 'Dimensions', {'nchar', 2, 'group', ngroup},...
             'Datatype', 'char')
    nccreate(fname, 'biomass', 'Dimensions', {'group', ngroup, 'X', nx,'time', nt})

    % fill coordinate information
    ncwriteatt(fname, 'time', 'long_name', 'time')
    ncwriteatt(fname, 'time', 'standard_name', 'time')
    ncwriteatt(fname, 'time', 'units', 'day')
    ncwriteatt(fname, 'time', 'calendar', 'noleap')
    ncwriteatt(fname, 'time', 'axis', 'T')
    ncwrite(fname, 'time', time)

    ncwriteatt(fname, 'X', 'long_name', 'longitude')
    ncwriteatt(fname, 'X', 'standard_name', 'longitude')
    ncwriteatt(fname, 'X', 'units', 'degrees_east')
    ncwriteatt(fname, 'X', 'axis', 'X')
    ncwrite(fname, 'X', X)

    ncwriteatt(fname, 'group', 'long_name', 'prognostic classes')
    ncwrite(fname, 'group', transpose(group))

    ncwrite(fname, 'biomass', permute(biomass, [3, 2, 1]))
end