%%%% THE MODEL
%%% DEMOGRAPHIC CALCULATIONS
function [Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT,ENVR] = sub_futbio_1meso(DY,ESM,GRD,Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT,param)

%%% If biomass < individual fish mass per grid cell, set all rates to zero? %%%

%%% ESM information
ENVR = get_ESM_1meso(ESM,GRD,param,DY);
ENVR.det = sub_neg(ENVR.det);
ENVR.Zm  = sub_neg(ENVR.Zm);
ENVR.dZm  = sub_neg(ENVR.dZm);

% Update benthic biomass with new detritus avail at that time step
[BENT.mass,BENT.pred] = sub_update_be(BENT.mass,param,ENVR.det,[Md.con_be,Ld.con_be],[Md.bio,Ld.bio]);
BENT.mass = sub_check(BENT.mass);

% Pelagic-demersal coupling
%Lp: fraction of time large piscivores spends in pelagic
%Ld: fraction of time large demersals spends in pelagic
if (param.pdc == 0)
    Lp.td = ones(param.NX,1);
    Ld.td = zeros(param.NX,1);
elseif (param.pdc == 1)
    Lp.td = ones(param.NX,1);
    Ld.td = sub_tdif_dem(ENVR.H,param,Mf.bio,Mp.bio,Md.bio,BENT.mass);
elseif (param.pdc == 2)
    Lp.td = sub_tdif_pel(ENVR.H,param,Mf.bio,Mp.bio,Md.bio);
    Ld.td = sub_tdif_dem(ENVR.H,param,Mf.bio,Mp.bio,Md.bio,BENT.mass);
else
    Lp.td = ones(param.NX,1);
    Ld.td = ones(param.NX,1);
end

% Average habitat temperature
Sf.thab = sub_thab(ENVR.Tp,ENVR.Tb,Sf.td);
Sp.thab = sub_thab(ENVR.Tp,ENVR.Tb,Sp.td);
Sd.thab = sub_thab(ENVR.Tp,ENVR.Tb,Sd.td);
Mf.thab = sub_thab(ENVR.Tp,ENVR.Tb,Mf.td);
Mp.thab = sub_thab(ENVR.Tp,ENVR.Tb,Mp.td);
Md.thab = sub_thab(ENVR.Tp,ENVR.Tb,Md.td);
Lp.thab = sub_thab(ENVR.Tp,ENVR.Tb,Lp.td);
Ld.thab = sub_thab(ENVR.Tp,ENVR.Tb,Ld.td);

% Metabolism
Sf.met = sub_met(ENVR.Tp,ENVR.Tb,Sf.td,param.M_s,param);
Sp.met = sub_met(ENVR.Tp,ENVR.Tb,Sp.td,param.M_s,param);
Sd.met = sub_met(ENVR.Tp,ENVR.Tb,Sd.td,param.M_s,param);
Mf.met = sub_met(ENVR.Tp,ENVR.Tb,Mf.td,param.M_m,param);
Mp.met = sub_met(ENVR.Tp,ENVR.Tb,Mp.td,param.M_m,param);
Md.met = sub_met(ENVR.Tp,ENVR.Tb,Md.td,param.M_m,param);
Lp.met = sub_met(ENVR.Tp,ENVR.Tb,Lp.td,param.M_l,param);
Ld.met = sub_met(ENVR.Tp,ENVR.Tb,Ld.td,param.M_l,param);

% Encounter rates
%           sub_enc(params,Tp     ,Tb     ,wgt      ,prey   ,tpel ,tprey,pref)
Sf.enc_zm = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_s,ENVR.Zm,Sf.td,Sf.td,param.MZ);
Sp.enc_zm = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_s,ENVR.Zm,Sp.td,Sp.td,param.MZ);
Sd.enc_zm = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_s,ENVR.Zm,Sd.td,Sd.td,param.MZ);

Mf.enc_zm = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_m,ENVR.Zm,Mf.td,Mf.td,param.MF_phi_MZ);
Mf.enc_f  = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_m,Sf.bio,Mf.td,Mf.td,param.MF_phi_S);
Mf.enc_p  = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_m,Sp.bio,Mf.td,Mf.td,param.MF_phi_S);
Mf.enc_d  = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_m,Sd.bio,Mf.td,Mf.td,param.MF_phi_S);

Mp.enc_zm = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_m,ENVR.Zm,Mp.td,Mp.td,param.MP_phi_MZ);
Mp.enc_f  = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_m,Sf.bio,Mp.td,Mp.td,param.MP_phi_S);
Mp.enc_p  = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_m,Sp.bio,Mp.td,Mp.td,param.MP_phi_S);
Mp.enc_d  = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_m,Sd.bio,Mp.td,Mp.td,param.MP_phi_S);

Md.enc_be = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_m,BENT.mass,Md.td,1-Md.td,param.MD_phi_BE);

Lp.enc_f  = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_l,Mf.bio,Lp.td,Lp.td,param.LP_phi_MF);
Lp.enc_p  = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_l,Mp.bio,Lp.td,Lp.td,param.LP_phi_MP);
Lp.enc_d  = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_l,Md.bio,Lp.td,1-Lp.td,param.LP_phi_MD);

Ld.enc_f  = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_l,Mf.bio,Ld.td,Ld.td,param.LD_phi_MF);
Ld.enc_p  = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_l,Mp.bio,Ld.td,Ld.td,param.LD_phi_MP);
Ld.enc_d  = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_l,Md.bio,Ld.td,1-Ld.td,param.LD_phi_MD);
Ld.enc_be = sub_enc(param,ENVR.Tp,ENVR.Tb,param.M_l,BENT.mass,Ld.td,1-Ld.td,param.LD_phi_BE);

% Consumption rates
Sf.con_zm = sub_cons(param,ENVR.Tp,ENVR.Tb,Sf.td,param.M_s,Sf.enc_zm);
Sp.con_zm = sub_cons(param,ENVR.Tp,ENVR.Tb,Sp.td,param.M_s,Sp.enc_zm);
Sd.con_zm = sub_cons(param,ENVR.Tp,ENVR.Tb,Sd.td,param.M_s,Sd.enc_zm);

Mf.con_zm = sub_cons(param,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m,[Mf.enc_zm,Mf.enc_f,Mf.enc_p,Mf.enc_d]);
Mf.con_f  = sub_cons(param,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m,[Mf.enc_f,Mf.enc_zm,Mf.enc_p,Mf.enc_d]);
Mf.con_p  = sub_cons(param,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m,[Mf.enc_p,Mf.enc_zm,Mf.enc_f,Mf.enc_d]);
Mf.con_d  = sub_cons(param,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m,[Mf.enc_d,Mf.enc_zm,Mf.enc_f,Mf.enc_p]);

Mp.con_zm = sub_cons(param,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m,[Mp.enc_zm,Mp.enc_f,Mp.enc_p,Mp.enc_d]);
Mp.con_f  = sub_cons(param,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m,[Mp.enc_f,Mp.enc_zm,Mp.enc_p,Mp.enc_d]);
Mp.con_p  = sub_cons(param,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m,[Mp.enc_p,Mp.enc_zm,Mp.enc_f,Mp.enc_d]);
Mp.con_d  = sub_cons(param,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m,[Mp.enc_d,Mp.enc_zm,Mp.enc_f,Mp.enc_p]);

Md.con_be = sub_cons(param,ENVR.Tp,ENVR.Tb,Md.td,param.M_m,Md.enc_be);

Lp.con_f  = sub_cons(param,ENVR.Tp,ENVR.Tb,Lp.td,param.M_l,[Lp.enc_f,Lp.enc_p,Lp.enc_d]);
Lp.con_p  = sub_cons(param,ENVR.Tp,ENVR.Tb,Lp.td,param.M_l,[Lp.enc_p,Lp.enc_f,Lp.enc_d]);
Lp.con_d  = sub_cons(param,ENVR.Tp,ENVR.Tb,Lp.td,param.M_l,[Lp.enc_d,Lp.enc_p,Lp.enc_f]);

Ld.con_f  = sub_cons(param,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l,[Ld.enc_f,Ld.enc_p,Ld.enc_d,Ld.enc_be]);
Ld.con_p  = sub_cons(param,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l,[Ld.enc_p,Ld.enc_f,Ld.enc_d,Ld.enc_be]);
Ld.con_d  = sub_cons(param,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l,[Ld.enc_d,Ld.enc_p,Ld.enc_f,Ld.enc_be]);
Ld.con_be = sub_cons(param,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l,[Ld.enc_be,Ld.enc_f,Ld.enc_p,Ld.enc_d]);

% Offline coupling
%MZ consumption cannot exceed amount lost to higher predation
%Track how much of mort is consumed with ENVR.fZl
[Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,ENVR.fZm] = ...
    sub_offline_zm(Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,Sf.bio,Sp.bio,Sd.bio,Mf.bio,Mp.bio,ENVR.dZm);

% Total consumption rates (could factor in handling times here; g m-2 d-1)
Sf.I = Sf.con_zm;
Sp.I = Sp.con_zm;
Sd.I = Sd.con_zm;
Mf.I = Mf.con_zm + Mf.con_f + Mf.con_p + Mf.con_d;
Mp.I = Mp.con_zm + Mp.con_f + Mp.con_p + Mp.con_d;
Md.I = Md.con_be;
Lp.I = Lp.con_f + Lp.con_p + Lp.con_d;
Ld.I = Ld.con_f + Ld.con_p + Ld.con_d + Ld.con_be;

% Consumption related to Cmax
Sf.clev = sub_clev(param,Sf.I,ENVR.Tp,ENVR.Tb,Sf.td,param.M_s);
Sp.clev = sub_clev(param,Sp.I,ENVR.Tp,ENVR.Tb,Sp.td,param.M_s);
Sd.clev = sub_clev(param,Sd.I,ENVR.Tp,ENVR.Tb,Sd.td,param.M_s);
Mf.clev = sub_clev(param,Mf.I,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m);
Mp.clev = sub_clev(param,Mp.I,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m);
Md.clev = sub_clev(param,Md.I,ENVR.Tp,ENVR.Tb,Md.td,param.M_m);
Lp.clev = sub_clev(param,Lp.I,ENVR.Tp,ENVR.Tb,Lp.td,param.M_l);
Ld.clev = sub_clev(param,Ld.I,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l);

% Death rates (g m-2 d-1)
Sf.die = Mp.con_f.*Mp.bio + Mf.con_f.*Mf.bio;
Sp.die = Mp.con_p.*Mp.bio + Mf.con_p.*Mf.bio;
Sd.die = Mp.con_d.*Mp.bio + Mf.con_d.*Mf.bio;
Mf.die = Lp.con_f.*Lp.bio + Ld.con_f.*Ld.bio;
Mp.die = Lp.con_p.*Lp.bio + Ld.con_p.*Ld.bio;
Md.die = Lp.con_d.*Lp.bio + Ld.con_d.*Ld.bio;

% predation rates (m-2 d-1)
Sf.pred = Sf.die ./ Sf.bio;
Sp.pred = Sp.die ./ Sp.bio;
Sd.pred = Sd.die ./ Sd.bio;
Mf.pred = Mf.die ./ Mf.bio;
Mp.pred = Mp.die ./ Mp.bio;
Md.pred = Md.die ./ Md.bio;

% Natural mortality rates
Sf.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Sf.td,param.M_s);
Sp.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Sp.td,param.M_s);
Sd.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Sd.td,param.M_s);
Mf.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m);
Mp.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m);
Md.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Md.td,param.M_m);
Lp.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Lp.td,param.M_l);
Ld.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l);

% Energy available for somatic growth nu
[Sf.nu, Sf.prod] = sub_nu(param,Sf.I,Sf.bio,Sf.met);
[Sp.nu, Sp.prod] = sub_nu(param,Sp.I,Sp.bio,Sp.met);
[Sd.nu, Sd.prod] = sub_nu(param,Sd.I,Sd.bio,Sd.met);
[Mf.nu, Mf.prod] = sub_nu(param,Mf.I,Mf.bio,Mf.met);
[Mp.nu, Mp.prod] = sub_nu(param,Mp.I,Mp.bio,Mp.met);
[Md.nu, Md.prod] = sub_nu(param,Md.I,Md.bio,Md.met);
[Lp.nu, Lp.prod] = sub_nu(param,Lp.I,Lp.bio,Lp.met);
[Ld.nu, Ld.prod] = sub_nu(param,Ld.I,Ld.bio,Ld.met);

% Maturation (note subscript on Kappa is larvae, juv, adult)
Sf.gamma = sub_gamma(param.K_l,param.Z_s,Sf.nu,Sf.die,Sf.bio,Sf.nmort,0,0);
Sp.gamma = sub_gamma(param.K_l,param.Z_s,Sp.nu,Sp.die,Sp.bio,Sp.nmort,0,0);
Sd.gamma = sub_gamma(param.K_l,param.Z_s,Sd.nu,Sd.die,Sd.bio,Sd.nmort,0,0);
Mf.gamma = sub_gamma(param.K_a,param.Z_m,Mf.nu,Mf.die,Mf.bio,Mf.nmort,param.dfrate,param.MFsel);
Mp.gamma = sub_gamma(param.K_j,param.Z_m,Mp.nu,Mp.die,Mp.bio,Mp.nmort,param.dfrate,param.MPsel);
Md.gamma = sub_gamma(param.K_j,param.Z_m,Md.nu,Md.die,Md.bio,Md.nmort,param.dfrate,param.MDsel);
Lp.gamma = sub_gamma(param.K_a,param.Z_l,Lp.nu,Lp.die,Lp.bio,Lp.nmort,param.dfrate,param.LPsel);
Ld.gamma = sub_gamma(param.K_a,param.Z_l,Ld.nu,Ld.die,Ld.bio,Ld.nmort,param.dfrate,param.LDsel);

% Reproduction (by adults only)
[Mf.gamma,Mf.nu,Mf.rep] = sub_rep(param.NX,Mf.gamma,Mf.nu,param.K_a);
[Lp.gamma,Lp.nu,Lp.rep] = sub_rep(param.NX,Lp.gamma,Lp.nu,param.K_a);
[Ld.gamma,Ld.nu,Ld.rep] = sub_rep(param.NX,Ld.gamma,Ld.nu,param.K_a);

% Recruitment (from smaller size class)
Sf.rec = sub_rec_larv(Mf.rep,Mf.bio,param.rfrac);
Sp.rec = sub_rec_larv(Lp.rep,Lp.bio,param.rfrac);
Sd.rec = sub_rec_larv(Ld.rep,Ld.bio,param.rfrac);
Mf.rec = sub_rec(Sf.gamma,Sf.bio);
Mp.rec = sub_rec(Sp.gamma,Sp.bio);
Md.rec = sub_rec(Sd.gamma,Sd.bio);
Lp.rec = sub_rec(Mp.gamma,Mp.bio);
Ld.rec = sub_rec(Md.gamma,Md.bio);

% Fishing by rate
[Sf.bio, Sf.caught, Sf.fmort] = sub_fishing_rate(Sf.bio,param.dfrate,0);
[Sp.bio, Sp.caught, Sp.fmort] = sub_fishing_rate(Sp.bio,param.dfrate,0);
[Sd.bio, Sd.caught, Sd.fmort] = sub_fishing_rate(Sd.bio,param.dfrate,0);
[Mf.bio, Mf.caught, Mf.fmort] = sub_fishing_rate(Mf.bio,param.dfrate,param.MFsel);
[Mp.bio, Mp.caught, Mp.fmort] = sub_fishing_rate(Mp.bio,param.dfrate,param.MPsel);
[Md.bio, Md.caught, Md.fmort] = sub_fishing_rate(Md.bio,param.dfrate,param.MDsel);
[Lp.bio, Lp.caught, Lp.fmort] = sub_fishing_rate(Lp.bio,param.dfrate,param.LPsel);
[Ld.bio, Ld.caught, Ld.fmort] = sub_fishing_rate(Ld.bio,param.dfrate,param.LDsel);

% Mass balance
Sf.bio = sub_update_fi(Sf.bio,Sf.rec,Sf.nu,Sf.rep,Sf.gamma,Sf.die,Sf.nmort,Sf.fmort);
Sp.bio = sub_update_fi(Sp.bio,Sp.rec,Sp.nu,Sp.rep,Sp.gamma,Sp.die,Sp.nmort,Sp.fmort);
Sd.bio = sub_update_fi(Sd.bio,Sd.rec,Sd.nu,Sd.rep,Sd.gamma,Sd.die,Sd.nmort,Sd.fmort);

Mf.bio = sub_update_fi(Mf.bio,Mf.rec,Mf.nu,Mf.rep,Mf.gamma,Mf.die,Mf.nmort,Mf.fmort);
Mp.bio = sub_update_fi(Mp.bio,Mp.rec,Mp.nu,Mp.rep,Mp.gamma,Mp.die,Mp.nmort,Mp.fmort);
Md.bio = sub_update_fi(Md.bio,Md.rec,Md.nu,Md.rep,Md.gamma,Md.die,Md.nmort,Md.fmort);

Lp.bio = sub_update_fi(Lp.bio,Lp.rec,Lp.nu,Lp.rep,Lp.gamma,Lp.die,Lp.nmort,Lp.fmort);
Ld.bio = sub_update_fi(Ld.bio,Ld.rec,Ld.nu,Ld.rep,Ld.gamma,Ld.die,Ld.nmort,Ld.fmort);


% Forward Euler checks for demographics and movement
Sf.bio=sub_check(Sf.bio);
Sp.bio=sub_check(Sp.bio);
Sd.bio=sub_check(Sd.bio);
Mf.bio=sub_check(Mf.bio);
Mp.bio=sub_check(Mp.bio);
Md.bio=sub_check(Md.bio);
Lp.bio=sub_check(Lp.bio);
Ld.bio=sub_check(Ld.bio);

%%% MNL DEBUG INFO
print_yaml = false;
% print_yaml = true;
if (print_yaml)
    fprintf('\nforcing:\n  values:\n')
    fprintf('    T_pelagic: %.15e\n', ENVR.Tp(1));
    fprintf('    T_bottom: %.15e\n', ENVR.Tb(1));
    fprintf('    zooC: %.15e\n', ENVR.Zm(1));
    fprintf('    poc_flux_bottom: %.15e\n', ENVR.det(1));
    fprintf('    zoo_mort: %.15e\n', ENVR.dZm(1));

    fprintf('\nT_habitat:\n  dimname: fish\n  values:\n')
    fprintf('    Sf: %.15e\n', Sf.thab(1))
    fprintf('    Sp: %.15e\n', Sp.thab(1))
    fprintf('    Sd: %.15e\n', Sd.thab(1))
    fprintf('    Mf: %.15e\n', Mf.thab(1))
    fprintf('    Mp: %.15e\n', Mp.thab(1))
    fprintf('    Md: %.15e\n', Md.thab(1))
    fprintf('    Lp: %.15e\n', Lp.thab(1))
    fprintf('    Ld: %.15e\n', Ld.thab(1))

    fprintf('\nmetabolism:\n  dimname: fish\n  values:\n')
    fprintf('    Sf: %.15e\n', Sf.met(1))
    fprintf('    Sp: %.15e\n', Sp.met(1))
    fprintf('    Sd: %.15e\n', Sd.met(1))
    fprintf('    Mf: %.15e\n', Mf.met(1))
    fprintf('    Mp: %.15e\n', Mp.met(1))
    fprintf('    Md: %.15e\n', Md.met(1))
    fprintf('    Lp: %.15e\n', Lp.met(1))
    fprintf('    Ld: %.15e\n', Ld.met(1))

    fprintf('\nencounter:\n  dimname: feeding_link\n  values:\n')
    fprintf('    Sf_Zoo: %.15e\n', Sf.enc_zm(1));
    fprintf('    Sp_Zoo: %.15e\n', Sp.enc_zm(1));
    fprintf('    Sd_Zoo: %.15e\n', Sd.enc_zm(1));
    fprintf('    Mf_Zoo: %.15e\n', Mf.enc_zm(1))
    fprintf('    Mf_Sf: %.15e\n', Mf.enc_f(1))
    fprintf('    Mf_Sp: %.15e\n', Mf.enc_p(1))
    fprintf('    Mf_Sd: %.15e\n', Mf.enc_d(1))
    fprintf('    Mp_Zoo: %.15e\n', Mp.enc_zm(1))
    fprintf('    Mp_Sf: %.15e\n', Mp.enc_f(1))
    fprintf('    Mp_Sp: %.15e\n', Mp.enc_p(1))
    fprintf('    Mp_Sd: %.15e\n', Mp.enc_d(1))
    fprintf('    Md_benthic_prey: %.15e\n', Md.enc_be(1))
    fprintf('    Lp_Mf: %.15e\n', Lp.enc_f(1))
    fprintf('    Lp_Mp: %.15e\n', Lp.enc_p(1))
    fprintf('    Lp_Md: %.15e\n', Lp.enc_d(1))
    fprintf('    Ld_Mf: %.15e\n', Ld.enc_f(1))
    fprintf('    Ld_Mp: %.15e\n', Ld.enc_p(1))
    fprintf('    Ld_Md: %.15e\n', Ld.enc_d(1))
    fprintf('    Ld_benthic_prey: %.15e\n', Ld.enc_be(1))

    fprintf('\nconsumption:\n  dimname: feeding_link\n  values:\n')
    fprintf('    Sf_Zoo: %.15e\n', Sf.con_zm(1))
    fprintf('    Sp_Zoo: %.15e\n', Sp.con_zm(1))
    fprintf('    Sd_Zoo: %.15e\n', Sd.con_zm(1))
    fprintf('    Mf_Zoo: %.15e\n', Mf.con_zm(1))
    fprintf('    Mf_Sf: %.15e\n', Mf.con_f(1))
    fprintf('    Mf_Sp: %.15e\n', Mf.con_p(1))
    fprintf('    Mf_Sd: %.15e\n', Mf.con_d(1))
    fprintf('    Mp_Zoo: %.15e\n', Mp.con_zm(1))
    fprintf('    Mp_Sf: %.15e\n', Mp.con_f(1))
    fprintf('    Mp_Sp: %.15e\n', Mp.con_p(1))
    fprintf('    Mp_Sd: %.15e\n', Mp.con_d(1))
    fprintf('    Md_benthic_prey: %.15e\n', Md.con_be(1))
    fprintf('    Lp_Mf: %.15e\n', Lp.con_f(1))
    fprintf('    Lp_Mp: %.15e\n', Lp.con_p(1))
    fprintf('    Lp_Md: %.15e\n', Lp.con_d(1))
    fprintf('    Ld_Mf: %.15e\n', Ld.con_f(1))
    fprintf('    Ld_Mp: %.15e\n', Ld.con_p(1))
    fprintf('    Ld_Md: %.15e\n', Ld.con_d(1))
    fprintf('    Ld_benthic_prey: %.15e\n', Ld.con_be(1))

    fprintf('\ningestion:\n  dimname: fish\n  values:\n')
    fprintf('    Sf: %.15e\n', Sf.I(1))
    fprintf('    Sp: %.15e\n', Sp.I(1))
    fprintf('    Sd: %.15e\n', Sd.I(1))
    fprintf('    Mf: %.15e\n', Mf.I(1))
    fprintf('    Mp: %.15e\n', Mp.I(1))
    fprintf('    Md: %.15e\n', Md.I(1))
    fprintf('    Lp: %.15e\n', Lp.I(1))
    fprintf('    Ld: %.15e\n', Ld.I(1))

    fprintf('\npredation:\n  dimname: fish\n  values:\n')
    fprintf('    Sf: %.15e\n', Sf.pred(1))
    fprintf('    Sp: %.15e\n', Sp.pred(1))
    fprintf('    Sd: %.15e\n', Sd.pred(1))
    fprintf('    Mf: %.15e\n', Mf.pred(1))
    fprintf('    Mp: %.15e\n', Mp.pred(1))
    fprintf('    Md: %.15e\n', Md.pred(1))
    fprintf('    Lp: %.15e\n', Lp.pred(1))
    fprintf('    Ld: %.15e\n', Ld.pred(1))

    fprintf('\nmortality:\n  dimname: fish\n  values:\n')
    fprintf('    Sf: %.15e\n', Sf.nmort(1))
    fprintf('    Sp: %.15e\n', Sp.nmort(1))
    fprintf('    Sd: %.15e\n', Sd.nmort(1))
    fprintf('    Mf: %.15e\n', Mf.nmort(1))
    fprintf('    Mp: %.15e\n', Mp.nmort(1))
    fprintf('    Md: %.15e\n', Md.nmort(1))
    fprintf('    Lp: %.15e\n', Lp.nmort(1))
    fprintf('    Ld: %.15e\n', Ld.nmort(1))

    fprintf('\nfishing:\n  dimname: fish\n  values:\n')
    fprintf('    Sf: %.15e\n', Sf.fmort(1))
    fprintf('    Sp: %.15e\n', Sp.fmort(1))
    fprintf('    Sd: %.15e\n', Sd.fmort(1))
    fprintf('    Mf: %.15e\n', Mf.fmort(1))
    fprintf('    Mp: %.15e\n', Mp.fmort(1))
    fprintf('    Md: %.15e\n', Md.fmort(1))
    fprintf('    Lp: %.15e\n', Lp.fmort(1))
    fprintf('    Ld: %.15e\n', Ld.fmort(1))

    fprintf('\navail_energy:\n  dimname: fish\n  values:\n')
    fprintf('    Sf: %.15e\n', Sf.nu(1))
    fprintf('    Sp: %.15e\n', Sp.nu(1))
    fprintf('    Sd: %.15e\n', Sd.nu(1))
    fprintf('    Mf: %.15e\n', Mf.nu(1))
    fprintf('    Mp: %.15e\n', Mp.nu(1))
    fprintf('    Md: %.15e\n', Md.nu(1))
    fprintf('    Lp: %.15e\n', Lp.nu(1))
    fprintf('    Ld: %.15e\n', Ld.nu(1))

    fprintf('\ngrowth:\n  dimname: fish\n  values:\n')
    fprintf('    Sf: %.15e\n', Sf.gamma(1))
    fprintf('    Sp: %.15e\n', Sp.gamma(1))
    fprintf('    Sd: %.15e\n', Sd.gamma(1))
    fprintf('    Mf: %.15e\n', Mf.gamma(1))
    fprintf('    Mp: %.15e\n', Mp.gamma(1))
    fprintf('    Md: %.15e\n', Md.gamma(1))
    fprintf('    Lp: %.15e\n', Lp.gamma(1))
    fprintf('    Ld: %.15e\n', Ld.gamma(1))

    fprintf('\nreproduction:\n  dimname: fish\n  values:\n')
    fprintf('    Sf: %.15e\n', Sf.rep(1))
    fprintf('    Sp: %.15e\n', Sp.rep(1))
    fprintf('    Sd: %.15e\n', Sd.rep(1))
    fprintf('    Mf: %.15e\n', Mf.rep(1))
    fprintf('    Mp: %.15e\n', Mp.rep(1))
    fprintf('    Md: %.15e\n', Md.rep(1))
    fprintf('    Lp: %.15e\n', Lp.rep(1))
    fprintf('    Ld: %.15e\n', Ld.rep(1))

    fprintf('\nrecruitment:\n  dimname: fish\n  values:\n')
    fprintf('    Sf: %.15e\n', Sf.rec(1))
    fprintf('    Sp: %.15e\n', Sp.rec(1))
    fprintf('    Sd: %.15e\n', Sd.rec(1))
    fprintf('    Mf: %.15e\n', Mf.rec(1))
    fprintf('    Mp: %.15e\n', Mp.rec(1))
    fprintf('    Md: %.15e\n', Md.rec(1))
    fprintf('    Lp: %.15e\n', Lp.rec(1))
    fprintf('    Ld: %.15e\n', Ld.rec(1))

    fprintf('\nbiomass:\n  dimname: group\n  values:\n')
    fprintf('    Sf: %.15e\n', Sf.bio(1))
    fprintf('    Sp: %.15e\n', Sp.bio(1))
    fprintf('    Sd: %.15e\n', Sd.bio(1))
    fprintf('    Mf: %.15e\n', Mf.bio(1))
    fprintf('    Mp: %.15e\n', Mp.bio(1))
    fprintf('    Md: %.15e\n', Md.bio(1))
    fprintf('    Lp: %.15e\n', Lp.bio(1))
    fprintf('    Ld: %.15e\n', Ld.bio(1))
    fprintf('    benthic_prey: %.15e\n', BENT.mass(1))
end % end if (print_yaml)

end
