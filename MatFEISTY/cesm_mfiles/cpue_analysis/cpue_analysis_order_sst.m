% CPUE analysis order of operations

%% 1. Calc anom timeseries for period of interest

%e.g.
calc_anom_ts_FishMIP_cpue_lme_82_15.m
calc_anom_ts_sst_lme_specyrs.m
calc_var_feisty_cesm_fosi_biom_lme_specyrs.m
calc_var_feisty_cesm_fosi_catch_lme_specyrs.m
calc_var_feisty_cesm_fosi_nu_lme_specyrs.m
calc_var_feisty_cesm_fosi_obsfish_biom_lme_specyrs.m
calc_var_feisty_cesm_fosi_obsfish_catch_lme_specyrs.m
calc_var_feisty_cesm_fosi_obsfish_nu_lme_specyrs.m
calc_var_inputs_cesm_fosi_lme_specyrs.m

%% 2. Calc corrs of each type of driver (satellite, inputs, fish)

%e.g. for sst 1982-2015
corr_lme_sstyr_inputs_cpue15
corr_lme_sstyr_feisty_cpue15
corr_lme_sstyr_obsfish_cpue15

corr_lme_sstyr_inputs_catch
corr_lme_sstyr_feisty_catch
corr_lme_sstyr_obsfish_catch

%% 3. Calc most significant corr across all drivers

%e.g. for sst 1982-2015
corr_lme_posfood_sstyr_catch_maxcorrtab.m
corr_lme_posfood_sstyr_inputs_catch_maxcorrtab.m
corr_lme_posfood_sstyr_inputs_feisty_catch_maxcorrtab.m
corr_lme_posfood_sstyr_inputs_obsfish_catch_maxcorrtab.m

corr_lme_posfood_sstyr_cpue15_maxcorrtab.m
corr_lme_posfood_sstyr_inputs_cpue15_maxcorrtab.m
corr_lme_posfood_sstyr_inputs_feisty_cpue15_maxcorrtab.m
corr_lme_posfood_sstyr_inputs_obsfish_cpue15_maxcorrtab.m

%% 4. Calc most significant lag of each driver 

%e.g.
corr_lme_posfood_sstyr_inputs_feisty_catch_mostsiglagtab.m
corr_lme_posfood_sstyr_inputs_obsfish_catch_mostsiglagtab.m

corr_lme_posfood_sstyr15_inputs_feisty_cpue15_mostsiglagtab.m
corr_lme_posfood_sstyr15_inputs_obsfish_cpue15_mostsiglagtab.m

%% 5. Create tables of # of LMEs with most sig driver 
% Table 1 of ms

table_lme_sstyr_input_posfood_feisty_obsfish_catch_sig.m
table_lme_sstyr_input_posfood_feisty_obsfish_cpue15_sig.m

%% 6. Plot - bar graphs and maps of all fishes combined 
% Figs 1 and 2 of ms

plot_corr_lme_sstyr_cpue15_comp_driver_maxcorrs
plot_corr_lme_sstyr_catch_comp_driver_maxcorrs.m


%% 7. Map R2 of each fn type, driver comp 
% Fig 3 in ms for all fishes
% fn types in Supp

map_corrR2_lme_sstyr_cpue15_comp_driver_maxcorrs
map_corrR2_lme_sstyr_catch_comp_driver_maxcorrs

%% 8. Plot - bar graphs and maps of each fn type 
% Supp figs

plot_corr_lme_sstyr_drivers_feisty_catch_types.m
plot_corr_lme_sstyr_drivers_feisty_cpue_types.m
plot_corr_lme_sstyr_drivers_obsfish_catch_types.m
plot_corr_lme_sstyr_drivers_obsfish_cpue_types.m

%% 9. 8-plot maps of corr coeff of all drivers (8) 
% for each fn type
% Supp figs

map_corrcoef_lme_posfood_sstyr_inputs_feisty_cpue15_types_v2
map_corrcoef_lme_posfood_sstyr_inputs_feisty_catch_types_v2

%% 10. 4-plots of feisty & obsfish results, CPUE and Catch together 
% bar plot and map of most corr driver
% not in ms, in ppts

plot_corr_lme_sstyr_feisty_obsfish_cpue15_catch.m

%% 11. 8-plot maps of max driver grouped as T, R, F
% fn types and const/obs fishing together
% Supp figs

map_lme_sstyr_cpue15_catch_groupdriver_maxcorrs

%% 12.4-plot map of R corr value of sim yeild with catch
% by fn type
% Supp fig

map_R2_lme_sstyr_obsfish2015_yield_catch_types.m
