%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
spinup_cesm = false;
fosi_cesm = true;

tic

if spinup_cesm
    Locs_CESM_4p4z_spinup()
%     Spinup_cesm()
%     Spinup_4p4z()
%     Spinup_4p4z_comb()
end
if fosi_cesm
%     Locs_CESM_4p4z()
    FOSI_cesm_catch()
%     CESM_4p4z()
%     CESM_4p4z_comb()
%     FOSI_cesm_climatol()
%     FOSI_cesm_varTemp()
%     FOSI_cesm_varFood()
end

toc
