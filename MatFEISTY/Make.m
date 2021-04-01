%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
spinup_cesm = true;
fosi_cesm = false;

tic

if spinup_cesm
    %Spinup_cesm()
    Spinup_4p4z()
end
if fosi_cesm
%     FOSI_cesm()
    cesm_4p4z()
%     FOSI_cesm_climatol()
%     FOSI_cesm_varTemp()
%     FOSI_cesm_varFood()
end

toc
