function save_biomass_nc(fname, biomass, time, X, group, GRD, ESM)

    % time = 1:nt;
    % X = 1:nx; % eventually pass this in as an argument
    % group = ['Sf'; 'Sp'; 'Sd'; 'Mf'; 'Mp'; 'Md'; 'Lp'; 'Ld'; 'bp'];

    % Make sure file exists!
    full_name = ['./model_output/', fname]
    if ~isfile(full_name)
        error(sprintf('ERROR: can not find %s',full_name));
    end

    % Write fields to file
    ncwrite(full_name, 'time', time)
    if isfield(GRD, 'LON')
        ncwrite(full_name, 'X', GRD.LON)
    else
        ncwrite(full_name, 'X', X)
    end
    if isfield(GRD, 'LAT')
        ncwrite(full_name, 'lat', GRD.LAT)
    end
    ncwrite(full_name, 'dep', GRD.Z)
    ncwrite(full_name, 'group', transpose(group))
    ncwrite(full_name, 'zooplankton', transpose(['Zoo']))

    % Forcing fields
    ncwrite(full_name, 'T_pelagic', ESM.Tp)
    ncwrite(full_name, 'T_bottom', ESM.Tb)
    ncwrite(full_name, 'poc_flux_bottom', ESM.det)
    ncwrite(full_name, 'zooC', ESM.Zm)
    ncwrite(full_name, 'zoo_mort', ESM.dZm)

    % Biomass (add additional prognostic variables here)
    ncwrite(full_name, 'biomass', permute(biomass, [3, 2, 1]))
end