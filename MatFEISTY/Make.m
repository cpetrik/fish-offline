%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
spinup_cesm = true;
fosi_cesm = false;

tic

if spinup_cesm
%    Locs_CESM_4p4z_spinup()
%     Spinup_cesm()
%     Spinup_cesm_quad()
    Spinup_cesm_quad_v2()
    Spinup_cesm_quad_v3()
%     Spinup_4p4z()
%     Spinup_4p4z_comb()
end
if fosi_cesm
%     Locs_CESM_4p4z()
%    FOSI_cesm()
%    FOSI_cesm_catch()
   FOSI_cesm_quad()
%   FOSI_cesm_quad_catch()
%     CESM_4p4z()
%     CESM_4p4z_comb()
%     FOSI_cesm_climatol()
%     FOSI_cesm_varTemp()
%     FOSI_cesm_varFood()
end

toc
