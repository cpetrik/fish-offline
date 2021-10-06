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

%% save
Fstat = array2table(fish_stat,'VariableNames',{'r','p','RMSE','Bias'},...
    'RowNames',{'SAU All Fish','SAU F','SAU P','SAU D','SAU Frac Pelagic',...
    'vanD Frac Pelagic','Stock All Fish'});
writetable(Fstat,[fpath 'FOSI_',mod,'obs_LME_all_ms_stats_' cfile2 '.csv'],'Delimiter',',',...
    'WriteRowNames',true)
save([fpath 'FOSI_',mod,'obs_LME_all_ms_stats_' cfile2 '.mat'],'fish_stat')
