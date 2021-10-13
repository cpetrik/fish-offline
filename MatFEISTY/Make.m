%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
spinup_cesm = false;
fosi_cesm = true;

tic

if spinup_cesm
%    Locs_CESM_4p4z_spinup() 
    Spinup_cesm()
%     Spinup_4p4z()
%     Spinup_4p4z_comb()
end
if fosi_cesm
%     Locs_CESM_4p4z()
   FOSI_cesm()
   FOSI_cesm_catch()
%     CESM_4p4z()
%     CESM_4p4z_comb()
end

toc
