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
    FOSI_cesm()
end

toc
