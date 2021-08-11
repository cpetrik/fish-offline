%%% Get CESM data
% Calculates quad mort loss from ESM
function ENVR = get_ESM_1meso_quad(ESM,GRD,param,DY)

    %% Get data
    ENVR.Tp(:,1)  = ESM.Tp(param.ID,DY);
    ENVR.Tb(:,1)  = ESM.Tb(param.ID,DY);
    ENVR.Zm(:,1)  = ESM.Zm(param.ID,DY);
    ENVR.det(:,1) = ESM.det(param.ID,DY);
    ENVR.fZm(:,1) = zeros(param.NX,1);
    ENVR.fB(:,1)  = zeros(param.NX,1);
    ENVR.H(:,1)   = GRD.Z(param.ID);
    ENVR.A(:,1)   = GRD.area(param.ID);
    
    %% back calc quad mort from temp, tot loss, and biomass
    Q10 = 2.0;
    T0Kelv = 273.15;
    Tref = 30.0;
    spd = 86400; %seconds per day
    dps = 1 ./ spd; %days per second
    parm_z_mort = 0.08 * dps;
    %parm_z_mort2 = 0.42 * dps;
    avg_f_loss_thres = 0.9167;
    loss_thres_zoo = 0.2;
    C_loss_thres = avg_f_loss_thres * loss_thres_zoo;
    
    % convert to orig units nmolC cm-2 and nmolC cm-2 s-1
    zooC = ENVR.Zm ./ (1e-9 * 1e4 * 12.01 * 9.0);
    zoo_loss = ESM.dZm(param.ID,DY) ./ (1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24);
    
    Tfn = Q10 .^ ( ( (ENVR.Tp + T0Kelv) - (Tref + T0Kelv) ) / 10.0 );
    Zprime = max((zooC - C_loss_thres),0.0);
    
    %zoo_loss = (zmort2 .* Zprime.^1.4) + (zmort .* Zprime);
    zmort = parm_z_mort .* Tfn;
    Lzoo_quad = zoo_loss - (zmort .* Zprime);
    
    % convert units back to gWW/m2/d
    Lzoo_quad = Lzoo_quad * 1e-9 * 1e4 * 12.01 * 9.0 * 60 * 60 * 24;
    
    %% calc quad mort from temp and biomass, don't use tot loss
%     Q10 = 2.0;
%     T0Kelv = 273.15;
%     Tref = 30.0;
%     spd = 86400; %seconds per day
%     dps = 1 ./ spd; %days per second
%     %parm_z_mort = 0.08 * dps;
%     parm_z_mort2 = 0.42 * dps;
%     avg_f_loss_thres = 0.9167;
%     loss_thres_zoo = 0.2;
%     C_loss_thres = avg_f_loss_thres * loss_thres_zoo;
%     
%     Tfn = Q10 .^ ( ( (ENVR.Tp + T0Kelv) - (Tref + T0Kelv) ) / 10.0 );
%     Zprime = max((ENVR.Zm - C_loss_thres),0.0);
%     
%     zmort2 = parm_z_mort2 .* Tfn;
%     Lzoo_quad = zmort2 .* Zprime.^1.4;
    
    ENVR.dZm(:,1) = Lzoo_quad;
    
end
