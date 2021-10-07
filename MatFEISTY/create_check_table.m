% Create check table for FEISTY functions
clear all
close all

%% ! Make core parameters/constants 
param = make_parameters_1meso(); 

%! Grid locations - just do 2, one coastal, one deep
param.NX = 2;
param.ID = 1:param.NX;
NX = param.NX;
ID = 1:param.NX;

%! How long to run the model
YEARS = 68;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%% ! Create a directory for output
exper = 'v13_';
[fname,simname] = sub_fname_cesm_fosi_exper(param,exper);

%% ! Initialize
init_sim = [exper simname];
load(['/Volumes/MIP/NC/CESM_MAPP/',simname '/Last_mo_spin_' init_sim '.mat']);
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%%
Sml_f.bio = round(nanmean(Sml_f.bio),2,'significant');
Sml_f.bio(2,1) = round(nanmean(Sml_f.bio),2,'significant');

Sml_p.bio = round(nanmean(Sml_p.bio),2,'significant');
Sml_p.bio(2,1) = round(nanmean(Sml_p.bio),2,'significant');

Sml_d.bio = round(nanmean(Sml_d.bio),2,'significant');
Sml_d.bio(2,1) = round(nanmean(Sml_d.bio),2,'significant');

Med_f.bio = round(nanmean(Med_f.bio),2,'significant');
Med_f.bio(2,1) = round(nanmean(Med_f.bio),2,'significant');

Med_p.bio = round(nanmean(Med_p.bio),2,'significant');
Med_p.bio(2,1) = round(nanmean(Med_p.bio),2,'significant');

Med_d.bio = round(nanmean(Med_d.bio),2,'significant');
Med_d.bio(2,1) = round(nanmean(Med_d.bio),2,'significant');

Lrg_p.bio = round(nanmean(Lrg_p.bio),2,'significant');
Lrg_p.bio(2,1) = round(nanmean(Lrg_p.bio),2,'significant');

Lrg_d.bio = round(nanmean(Lrg_d.bio),2,'significant');
Lrg_d.bio(2,1) = round(nanmean(Lrg_d.bio),2,'significant');

BENT.bio = round(nanmedian(BENT.mass),2,'significant');
BENT.bio(2,1) = round(nanmedian(BENT.mass),2,'significant');
BENT.mass = BENT.bio;

%%
Sf = Sml_f;
Sp = Sml_p;
Sd = Sml_d;
Mf = Med_f;
Mp = Med_p;
Md = Med_d;
Lp = Lrg_p;
Ld = Lrg_d;

%%
load('/Volumes/MIP/GCM_DATA/CESM/FOSI/Data_grid_POP_gx1v6.mat','GRD');
YR=1;
ti = num2str(YR);
load(['/Volumes/MIP/GCM_DATA/CESM/FOSI/Data_cesm_fosi_v6_daily_',ti,'.mat'],'ESM');
    
cst = GRD.Z<=200;
esm.Tp = round(nanmean(ESM.Tp(cst,1)),2,'significant');
esm.Tp(2,1) = round(nanmean(ESM.Tp(~cst,1)),2,'significant');
esm.Tb = round(nanmean(ESM.Tb(cst,1)),2,'significant');
esm.Tb(2,1) = round(nanmean(ESM.Tb(~cst,1)),2,'significant');
esm.Zm = round(nanmean(ESM.Zm(cst,1)),2,'significant');
esm.Zm(2,1) = round(nanmean(ESM.Zm(~cst,1)),2,'significant');
esm.dZm = round(nanmean(ESM.dZm(cst,1)),2,'significant');
esm.dZm(2,1) = round(nanmean(ESM.dZm(~cst,1)),2,'significant');
esm.det = round(nanmean(ESM.det(cst,1)),2,'significant');
esm.det(2,1) = round(nanmean(ESM.det(~cst,1)),2,'significant');

clear ESM
ESM = esm;

gridZ = round(nanmean(GRD.Z(cst,1)),2,'significant');
gridZ(2,1) = round(nanmean(GRD.Z(~cst,1)),2,'significant');
GRD.Z = gridZ;

%% 
DY = 1;
% [Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT,ENVR] = sub_futbio_1meso_mzpref(DY,ESM,GRD,...
%     Sf,Sp,Sd,Mf,Mp,Md,Lp,Ld,BENT,param);

ENVR = get_ESM_1meso(ESM,GRD,param,DY);

[BENT.mass,BENT.pred] = sub_update_be(BENT.mass,param,ENVR.det,[Md.con_be,Ld.con_be],[Md.bio,Ld.bio]);

Ld.td = sub_tdif_dem(ENVR.H,param,Mf.bio,Mp.bio,Md.bio,BENT.mass);

%% Metabolism
Sf.met = sub_met(ENVR.Tp,ENVR.Tb,Sf.td,param.M_s,param);
Sp.met = sub_met(ENVR.Tp,ENVR.Tb,Sp.td,param.M_s,param);
Sd.met = sub_met(ENVR.Tp,ENVR.Tb,Sd.td,param.M_s,param);
Mf.met = sub_met(ENVR.Tp,ENVR.Tb,Mf.td,param.M_m,param);
Mp.met = sub_met(ENVR.Tp,ENVR.Tb,Mp.td,param.M_m,param);
Md.met = sub_met(ENVR.Tp,ENVR.Tb,Md.td,param.M_m,param);
Lp.met = sub_met(ENVR.Tp,ENVR.Tb,Lp.td,param.M_l,param);
Ld.met = sub_met(ENVR.Tp,ENVR.Tb,Ld.td,param.M_l,param);

%% Encounter rates
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

%% Consumption rates
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

%% Offline coupling
[Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,ENVR.fZm] = ...
    sub_offline_zm(Sf.con_zm,Sp.con_zm,Sd.con_zm,Mf.con_zm,Mp.con_zm,...
    Sf.bio,Sp.bio,Sd.bio,Mf.bio,Mp.bio,ENVR.dZm);

%% Total consumption rates (could factor in handling times here; g m-2 d-1)
Sf.I = Sf.con_zm;
Sp.I = Sp.con_zm;
Sd.I = Sd.con_zm;
Mf.I = Mf.con_zm + Mf.con_f + Mf.con_p + Mf.con_d;
Mp.I = Mp.con_zm + Mp.con_f + Mp.con_p + Mp.con_d;
Md.I = Md.con_be;
Lp.I = Lp.con_f + Lp.con_p + Lp.con_d;
Ld.I = Ld.con_f + Ld.con_p + Ld.con_d + Ld.con_be;

%% Consumption related to Cmax
Sf.clev = sub_clev(param,Sf.I,ENVR.Tp,ENVR.Tb,Sf.td,param.M_s);
Sp.clev = sub_clev(param,Sp.I,ENVR.Tp,ENVR.Tb,Sp.td,param.M_s);
Sd.clev = sub_clev(param,Sd.I,ENVR.Tp,ENVR.Tb,Sd.td,param.M_s);
Mf.clev = sub_clev(param,Mf.I,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m);
Mp.clev = sub_clev(param,Mp.I,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m);
Md.clev = sub_clev(param,Md.I,ENVR.Tp,ENVR.Tb,Md.td,param.M_m);
Lp.clev = sub_clev(param,Lp.I,ENVR.Tp,ENVR.Tb,Lp.td,param.M_l);
Ld.clev = sub_clev(param,Ld.I,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l);

%temp save here

%% Death rates (g m-2 d-1)
Sf.die = Mp.con_f.*Mp.bio + Mf.con_f.*Mf.bio;
Sp.die = Mp.con_p.*Mp.bio + Mf.con_p.*Mf.bio;
Sd.die = Mp.con_d.*Mp.bio + Mf.con_d.*Mf.bio;
Mf.die = Lp.con_f.*Lp.bio + Ld.con_f.*Ld.bio;
Mp.die = Lp.con_p.*Lp.bio + Ld.con_p.*Ld.bio;
Md.die = Lp.con_d.*Lp.bio + Ld.con_d.*Ld.bio;

%% biomass-specific predation rates (m-2 d-1)
Sf.pred = Sf.die ./ Sf.bio;
Sp.pred = Sp.die ./ Sp.bio;
Sd.pred = Sd.die ./ Sd.bio;
Mf.pred = Mf.die ./ Mf.bio;
Mp.pred = Mp.die ./ Mp.bio;
Md.pred = Md.die ./ Md.bio;

%% Natural mortality rates
Sf.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Sf.td,param.M_s);
Sp.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Sp.td,param.M_s);
Sd.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Sd.td,param.M_s);
Mf.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Mf.td,param.M_m);
Mp.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Mp.td,param.M_m);
Md.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Md.td,param.M_m);
Lp.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Lp.td,param.M_l);
Ld.nmort = sub_nmort(param,ENVR.Tp,ENVR.Tb,Ld.td,param.M_l);

%% Energy available for somatic growth nu
[Sf.nu, Sf.prod] = sub_nu(param,Sf.I,Sf.bio,Sf.met);
[Sp.nu, Sp.prod] = sub_nu(param,Sp.I,Sp.bio,Sp.met);
[Sd.nu, Sd.prod] = sub_nu(param,Sd.I,Sd.bio,Sd.met);
[Mf.nu, Mf.prod] = sub_nu(param,Mf.I,Mf.bio,Mf.met);
[Mp.nu, Mp.prod] = sub_nu(param,Mp.I,Mp.bio,Mp.met);
[Md.nu, Md.prod] = sub_nu(param,Md.I,Md.bio,Md.met);
[Lp.nu, Lp.prod] = sub_nu(param,Lp.I,Lp.bio,Lp.met);
[Ld.nu, Ld.prod] = sub_nu(param,Ld.I,Ld.bio,Ld.met);

%% Maturation (note subscript on Kappa is larvae, juv, adult)
Sf.gamma = sub_gamma(param.K_l,param.Z_s,Sf.nu,Sf.die,Sf.bio,Sf.nmort,param.dfrate,0);
Sp.gamma = sub_gamma(param.K_l,param.Z_s,Sp.nu,Sp.die,Sp.bio,Sp.nmort,param.dfrate,0);
Sd.gamma = sub_gamma(param.K_l,param.Z_s,Sd.nu,Sd.die,Sd.bio,Sd.nmort,param.dfrate,0);
Mf.gamma = sub_gamma(param.K_a,param.Z_m,Mf.nu,Mf.die,Mf.bio,Mf.nmort,param.dfrate,param.MFsel);
Mp.gamma = sub_gamma(param.K_j,param.Z_m,Mp.nu,Mp.die,Mp.bio,Mp.nmort,param.dfrate,param.MPsel);
Md.gamma = sub_gamma(param.K_j,param.Z_m,Md.nu,Md.die,Md.bio,Md.nmort,param.dfrate,param.MDsel);
Lp.gamma = sub_gamma(param.K_a,param.Z_l,Lp.nu,Lp.die,Lp.bio,Lp.nmort,param.dfrate,param.LPsel);
Ld.gamma = sub_gamma(param.K_a,param.Z_l,Ld.nu,Ld.die,Ld.bio,Ld.nmort,param.dfrate,param.LDsel);

%% Reproduction (by adults only)
[Mf.gamma,Mf.nu,Mf.rep] = sub_rep(param.NX,Mf.gamma,Mf.nu,param.K_a);
[Lp.gamma,Lp.nu,Lp.rep] = sub_rep(param.NX,Lp.gamma,Lp.nu,param.K_a);
[Ld.gamma,Ld.nu,Ld.rep] = sub_rep(param.NX,Ld.gamma,Ld.nu,param.K_a);

%% Recruitment (from smaller size class)
Sf.rec = sub_rec_larv(Mf.rep,Mf.bio,param.rfrac);
Sp.rec = sub_rec_larv(Lp.rep,Lp.bio,param.rfrac);
Sd.rec = sub_rec_larv(Ld.rep,Ld.bio,param.rfrac);
Mf.rec = sub_rec(Sf.gamma,Sf.bio);
Mp.rec = sub_rec(Sp.gamma,Sp.bio);
Md.rec = sub_rec(Sd.gamma,Sd.bio);
Lp.rec = sub_rec(Mp.gamma,Mp.bio);
Ld.rec = sub_rec(Md.gamma,Md.bio);

%%
Sf.bio_init=Sf.bio;
Sp.bio_init=Sp.bio;
Sd.bio_init=Sd.bio;
Mf.bio_init=Mf.bio;
Mp.bio_init=Mp.bio;
Md.bio_init=Md.bio;
Lp.bio_init=Lp.bio;
Ld.bio_init=Ld.bio;

%% Mass balance
Sf.bio = sub_update_fi(Sf.bio,Sf.rec,Sf.nu,Sf.rep,Sf.gamma,Sf.die,Sf.nmort);
Sp.bio = sub_update_fi(Sp.bio,Sp.rec,Sp.nu,Sp.rep,Sp.gamma,Sp.die,Sp.nmort);
Sd.bio = sub_update_fi(Sd.bio,Sd.rec,Sd.nu,Sd.rep,Sd.gamma,Sd.die,Sd.nmort);

Mf.bio = sub_update_fi(Mf.bio,Mf.rec,Mf.nu,Mf.rep,Mf.gamma,Mf.die,Mf.nmort);
Mp.bio = sub_update_fi(Mp.bio,Mp.rec,Mp.nu,Mp.rep,Mp.gamma,Mp.die,Mp.nmort);
Md.bio = sub_update_fi(Md.bio,Md.rec,Md.nu,Md.rep,Md.gamma,Md.die,Md.nmort);

Lp.bio = sub_update_fi(Lp.bio,Lp.rec,Lp.nu,Lp.rep,Lp.gamma,Lp.die,Lp.nmort);
Ld.bio = sub_update_fi(Ld.bio,Ld.rec,Ld.nu,Ld.rep,Ld.gamma,Ld.die,Ld.nmort);

%% Fishing by rate
[Mf.bio, Mf.caught, Mf.fmort] = sub_fishing_rate(Mf.bio,param.dfrate,param.MFsel);
[Mp.bio, Mp.caught, Mp.fmort] = sub_fishing_rate(Mp.bio,param.dfrate,param.MPsel);
[Md.bio, Md.caught, Md.fmort] = sub_fishing_rate(Md.bio,param.dfrate,param.MDsel);
[Lp.bio, Lp.caught, Lp.fmort] = sub_fishing_rate(Lp.bio,param.dfrate,param.LPsel);
[Ld.bio, Ld.caught, Ld.fmort] = sub_fishing_rate(Ld.bio,param.dfrate,param.LDsel);

%%

%% save
Fstat = array2table(fish_stat,'VariableNames',{'r','p','RMSE','Bias'},...
    'RowNames',{'SAU All Fish','SAU F','SAU P','SAU D','SAU Frac Pelagic',...
    'vanD Frac Pelagic','Stock All Fish'});
% writetable(Fstat,[fpath 'FOSI_',mod,'obs_LME_all_ms_stats_' cfile2 '.csv'],'Delimiter',',',...
%     'WriteRowNames',true)
% save([fpath 'FOSI_',mod,'obs_LME_all_ms_stats_' cfile2 '.mat'],'fish_stat')
