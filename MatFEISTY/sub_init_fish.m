%============== INITIAL CONDITIONS =============%
function [Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT] = sub_init_fish(ID,Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT)

% Initialize variables not saved from spinup
% (Only biomass saved from spinup)

%===== VARIABLES =====%
    %%%! Number of spatial cells
    NX = length(ID);

    %%%! Fraction of time in pelagic
    Sf.td = ones(NX,1);
    Sp.td = ones(NX,1);
    Sd.td = ones(NX,1);
    Mf.td = ones(NX,1);
    Mp.td = ones(NX,1);
    Md.td = zeros(NX,1);
    Lp.td = ones(NX,1);
    Ld.td = zeros(NX,1);

    %%%! Rates, etc.
    nzero = {'met' 'enc_f' 'enc_p' 'enc_d' 'enc_zm' 'enc_zl' 'enc_be' ...
        'con_f' 'con_p' 'con_d' 'con_zm' 'con_zl' 'con_be' 'I' 'die' 'pred' ...
        'nmort' 'prod' 'nu' 'gamma' 'rep' 'rec' 'clev' 'caught' 'fmort' 'thab'};
    for i = 1 : length(nzero)
        Sf.(nzero{i}) = zeros(NX,1);
        Sp.(nzero{i}) = zeros(NX,1);
        Sd.(nzero{i}) = zeros(NX,1);
        Mf.(nzero{i}) = zeros(NX,1);
        Mp.(nzero{i}) = zeros(NX,1);
        Md.(nzero{i}) = zeros(NX,1);
        Lp.(nzero{i}) = zeros(NX,1);
        Ld.(nzero{i}) = zeros(NX,1);
    end
    
   
    %%%! Detritus
    BENT.pred = zeros(NX,1);

end
