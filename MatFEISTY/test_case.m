%%%%!! RUN HISTORIC FOR ALL LOCATIONS
% $ matlab -nodisplay -nosplash - nodesktop
% >> run('test_case.m')
function test_case()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Make core parameters/constants
param = make_params_testcase();

%! Idealized bathymetry
load('./input_files/Grid_test_forcing.mat',...
    'GRD');
param.NX = length(GRD.Z);
param.ID = 1:param.NX;
NX = param.NX;
ID = 1:param.NX;

%! Idealized forcing
load('./input_files/Data_cyclic_test_forcing.mat',...
    'ESM');

%! How long to run the model
YEARS = 1;
DAYS = 365;
% DAYS = 2;

%! Create a directory for output
exper = 'v3_move_updateB_';
[fname,simname] = sub_fname_testcase_exper(param,exper);

%! Initialize
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_testcase(ID,DAYS);

%! Storage variables
biomass             = NaN*ones(DAYS,NX,9);
T_habitat           = NaN*ones(DAYS,NX,9);
ingestion_rate      = NaN*ones(DAYS,NX,9);
predation_flux      = NaN*ones(DAYS,NX,9);
predation_rate      = NaN*ones(DAYS,NX,9);
metabolism_rate     = NaN*ones(DAYS,NX,9);
mortality_rate      = NaN*ones(DAYS,NX,9);
energy_avail_rate   = NaN*ones(DAYS,NX,9);
growth_rate         = NaN*ones(DAYS,NX,9);
reproduction_rate   = NaN*ones(DAYS,NX,9);
recruitment_flux    = NaN*ones(DAYS,NX,9);
fish_catch_rate     = NaN*ones(DAYS,NX,9);

%% %%%%%%%%%%%%%%%%%%%% Run the Model
MNT = 0;
%! Run model with no fishing
for YR = 1:YEARS % years
    for DAY = 1:param.DT:DAYS % days

        %%%! Future time step
        DY = int64(ceil(DAY));
        [num2str(YR),' , ', num2str(mod(DY,365))]

        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio_1meso(DY,ESM,GRD,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);

        %! Store
        biomass(DY,:,:) = [Sml_f.bio, Sml_p.bio, Sml_d.bio,...
                          Med_f.bio, Med_p.bio, Med_d.bio,...
                          Lrg_p.bio, Lrg_d.bio, BENT.mass];
        T_habitat(DY,:,1:8) = [Sml_f.thab, Sml_p.thab, Sml_d.thab,...
                               Med_f.thab, Med_p.thab, Med_d.thab,...
                               Lrg_p.thab, Lrg_d.thab];
        ingestion_rate(DY,:,1:8) = [Sml_f.I, Sml_p.I, Sml_d.I,...
                                   Med_f.I, Med_p.I, Med_d.I,...
                                   Lrg_p.I, Lrg_d.I];
        predation_flux(DY,:,1:8) = [Sml_f.die, Sml_p.die, Sml_d.die,...
                                   Med_f.die, Med_p.die, Med_d.die,...
                                   Lrg_p.die, Lrg_d.die];
        predation_rate(DY,:,1:8) = [Sml_f.pred, Sml_p.pred, Sml_d.pred,...
                                   Med_f.pred, Med_p.pred, Med_d.pred,...
                                   Lrg_p.pred, Lrg_d.pred];
        metabolism_rate(DY,:,1:8) = [Sml_f.met, Sml_p.met, Sml_d.met,...
                                   Med_f.met, Med_p.met, Med_d.met,...
                                   Lrg_p.met, Lrg_d.met];
        mortality_rate(DY,:,1:8) = [Sml_f.nmort, Sml_p.nmort, Sml_d.nmort,...
                                   Med_f.nmort, Med_p.nmort, Med_d.nmort,...
                                   Lrg_p.nmort, Lrg_d.nmort];
        energy_avail_rate(DY,:,1:8) = [Sml_f.nu, Sml_p.nu, Sml_d.nu,...
                                      Med_f.nu, Med_p.nu, Med_d.nu,...
                                      Lrg_p.nu, Lrg_d.nu];
        growth_rate(DY,:,1:8) = [Sml_f.gamma, Sml_p.gamma, Sml_d.gamma,...
                                 Med_f.gamma, Med_p.gamma, Med_d.gamma,...
                                 Lrg_p.gamma, Lrg_d.gamma];
        reproduction_rate(DY,:,1:8) = [Sml_f.rep, Sml_p.rep, Sml_d.rep,...
                                      Med_f.rep, Med_p.rep, Med_d.rep,...
                                      Lrg_p.rep, Lrg_d.rep];
        recruitment_flux(DY,:,1:8) = [Sml_f.rec, Sml_p.rec, Sml_d.rec,...
                                      Med_f.rec, Med_p.rec, Med_d.rec,...
                                      Lrg_p.rec, Lrg_d.rec];
        fish_catch_rate(DY,:,1:8) = [Sml_f.fmort, Sml_p.fmort, Sml_d.fmort,...
                                     Med_f.fmort, Med_p.fmort, Med_d.fmort,...
                                     Lrg_p.fmort, Lrg_d.fmort];

    end %Days

end %Years

%%% Save
save([fname '_test_case.mat'],...
    'biomass','T_habitat','ingestion_rate','predation_flux','predation_rate',...
    'metabolism_rate','mortality_rate','energy_avail_rate','growth_rate',...
    'reproduction_rate','recruitment_flux','fish_catch_rate')

end
