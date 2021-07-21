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
    zoo_loss = ESM.dZm(param.ID,DY);
    Tfn = exp(-4000 .* ( (1./(ENVR.Tp+273.15)) - (1./303.15) ));
    Zprime = max((ENVR.Zm - 0.01),0);
    Lzoo_quad = zoo_loss - (Tfn .* 0.1 .* Zprime);
    
    %%
    ENVR.dZm(:,1) = Lzoo_quad;
    
end
