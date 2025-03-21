%% FEISTY Make file

clear 
close all

%%%%!! EXPERIMENTS
testc       = false;
spinup_cesm = false;
spinup_hr   = false;
hr_cesm     = false;
fosi_cesm   = true;
dple_cesm   = false;
cesm_4p2z   = false;

tic

if testc
    %test_case()
    test_locs3()
end
if spinup_cesm
%     Locs_CESM_4p4z_spinup()
%     Spinup_cesm()
%     Spinup_FOSI_climatol()
%     Spinup_FOSI_varFood()
%     Spinup_FOSI_varTemp()
%     Spinup_4p4z()
%     Spinup_4p4z_comb()
    Spinup_4p2z_1deg()
end
if spinup_hr
    Spinup_cesm_hr_CCE()
end
if hr_cesm
    HR_cesm()
    HR_cesm_catch()
end
if fosi_cesm
%     Locs_CESM_4p4z()
%     FOSI_cesm()
%     FOSI_cesm_catch()
    FOSI_cesm_nu()
    FOSI_cesm_rec_gamIn()
%     FOSI_cesm_search()
end
if cesm_4p2z
%     CESM_4p4z()
%     CESM_4p4z_comb()
    CESM_4p2z_1deg()
    CESM_4p2z_1deg_catch()
end
if dple_cesm
    Interp_run_DPLE_cesm()
end

toc
