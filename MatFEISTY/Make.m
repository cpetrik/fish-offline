%% FEISTY Make file

clear all
close all

%%%%!! EXPERIMENTS
spinup_cesm = false;
fosi_cesm = true;

tic

if spinup_cesm
    Spinup_cesm()
end
if fosi_cesm
    %FOSI_cesm()
    FOSI_cesm_climatol()
    FOSI_cesm_varTemp()
    FOSI_cesm_varFood()
end

toc
