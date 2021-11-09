%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
spinup_cesm = true;
fosi_cesm = false;

tic

if spinup_cesm
%     Locs_CESM_4p4z_spinup()
    Spinup_cesm()
    Spinup_FOSI_climatol()
    Spinup_FOSI_varFood()
    Spinup_FOSI_varTemp()
%     Spinup_4p4z()
%     Spinup_4p4z_comb()
end
if fosi_cesm
%     Locs_CESM_4p4z()
    FOSI_cesm()
    FOSI_cesm_catch()
    FOSI_cesm_climatol()
    FOSI_cesm_varFood()
    FOSI_cesm_varTemp()
%     FOSI_cesm_search()
%     CESM_4p4z()
%     CESM_4p4z_comb()
end

toc
