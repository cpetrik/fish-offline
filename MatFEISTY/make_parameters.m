%============== Parameters of the model =============%
%============= PARAMETER TYPE ==========%
function make_parameters()
    global DT PI_be_cutoff pdc L_s L_m L_l M_s M_m M_l L_zm L_zl 
    global Z_s Z_m Z_l Lambda K_l K_j K_a h gam kt bpow
    global bent_eff rfrac D J Sm A benc bcmx amet 
    global Tu_s Tu_m Tu_l Nat_mrt MORT
    global MF_phi_MZ MF_phi_LZ MF_phi_S MP_phi_MZ MP_phi_LZ MP_phi_S MD_phi_BE 
    global LP_phi_MF LP_phi_MP LP_phi_MD LD_phi_MF LD_phi_MP LD_phi_MD LD_phi_BE  
    global MFsel MPsel MDsel LPsel LDsel Jsel
    global tstep K frate kc ke
    
    %! Integration parameters
    DT = 1.0;       % time step
    tstep = 1.0;    % time step in hours for adv-diff
    
    % define diffusivity
    K = 600.0;
    
    %! Which fishes harvested
    MFsel = 1;
    LPsel = 1;
    LDsel = 1;
    Jsel  = 0.1;
    MPsel = Jsel * LPsel;
    MDsel = Jsel * LDsel;
    
    %! Benthic-pelagic coupling cutoff (depth, m)
    PI_be_cutoff = 200;
    % 0:no coupling; 1:demersal coupled only; 2:pelagic & demersal coupled;
    pdc = 1;

    %!Individual Mass (g) = geometric mean
    M_s = 10^((log10(0.001)+log10(0.5))/2);  %0.0224
    M_m = 10^((log10(0.5)+log10(250))/2);    %11.1803
    M_l = 10^((log10(250)+log10(125000))/2); %5.5902e3
    %logspace(-3,5.0969,7) %gives end points and mid points

    %! Body lengths (mm)
    % Convert from mm to cm and use their const coeff = 0.01g/cm3
    L_s = 10.0 * (M_s/0.01)^(1/3); % small   4.6-36.8 mm, 13.1 mm
    L_m = 10.0 * (M_m/0.01)^(1/3); % medium 36.8-292 mm, 10.4 cm
    L_l = 10.0 * (M_l/0.01)^(1/3); % large  0.292-2.32 m, 0.82 m

    %! Median Zooplankton size in mm from Charlie/COBALT
    %! Median Zooplankton body mass in g wet weight
    % eq from James Watkins, Lars Rudstam and Kristen Holeck in dry weight
    % convert to wet weight with 1g dry = 9 g wet
    L_zm = 10^((log10(0.2)+log10(2))/2); % lengths (ESD)
    L_zl = 10^((log10(2)+log10(20))/2);

    %! Ratio of initial and final body sizes per size-class
    Z_s = 0.001/0.5;
    Z_m = 0.5/250;
    Z_l = 250/125000;

    %%%! Assimilation efficiency lambda (constant across everything)
    Lambda = 0.7;

    %%%! Kappa rule K as a function of body size
    % K = fraction of energy consumed diverted to somatic growth
    % subscript is larvae, juv, adult)
    K_l = 1;
    K_j = 1;
    K_a = 0.5;

    %%%! Metabolism constants (activity and basal)
    amet = 4;       % coeff on met
    h = 20;         % coeff on Cmax
    gam = 70;       % coeff on search area
    kc = 0.063;     % coeff on cmax T-dep fn 
    ke = 0.063;     % coeff on enc T-dep fn 
    kt = 0.0855;    % coeff on met T-dep fn (orig 0.063)
    bpow = 0.175;   % power on metab fn (orig 0.25)
    benc = 0.20;    % power on enc fn (orig 0.20)
    bcmx = 0.25;    % power on cmax fn (orig 0.25)
    
    %%%! Transfer efficiency of detritus to benthic prey
    bent_eff = 0.075;
    
    %%%! Reproductive efficiency
    rfrac = 0.01; 

    %! Fraction of time spent swimming (from Van Leeuwen)
    Tu_s = 1.0;
    Tu_m = 1.0; %0.5
    Tu_l = 1.0; %0.1

    %%%! Background mortality
    Nat_mrt = 0.1/365; 
    %0=none, 1=constant, 6=const wgt, T-dep, 7=const T, wgt-dep
    %2=Hartvig T-dep, 3=mizer T-dep, 4=J&C T-dep, 5=P&W T-dep
    MORT = 1;   

    %%%! Diet Preference Phi
    % The predator prey mass ratio is assumed 3 orders of mag, i.e. 1000, i.e. one step down
    % Because Medium fishes are 2 sizes bigger than Medium zoo, pref = 0.1
    % We don't have a pred-prey matrix anymore, we are instead explicit about who eats who:
    %-----
    %small forage fish eats medium zoo
    %small piscivores eats medium zoo
    %small detritivore eats medium zoo
    %medium forage fish eats medium & large zoo, all small fishes
    %medium piscivore eats medium & large zoo, all small fishes
    %medium detritivore eats detritus
    %large piscivore eats medium forage fish, medium piscivore, medium detritivore
    %large detritivore eats detritus, medium forage fish, medium piscivore, medium detrivore
    
    Sm = 0.25;  %Feeding 2 sizes down
    J = 1.0;    %Juvenile feeding reduction
    D = 0.75;   %Demersal feeding in pelagic reduction
    A = 0.5;    %Adult predation reduction %*****

    MF_phi_MZ = Sm;
    MF_phi_LZ = 1.0;
    MF_phi_S  = 1.0;
    
    MP_phi_MZ = Sm*J;
    MP_phi_LZ = J;
    MP_phi_S  = J;
    
    MD_phi_BE = 1.0;
    
    LP_phi_MF = 1.0*A;
    LP_phi_MP = 1.0;
    LP_phi_MD = 0.0;
    
    LD_phi_MF = D*A;
    LD_phi_MP = D;
    LD_phi_MD = 1.0;
    LD_phi_BE = 1.0;
    
%-----
end
