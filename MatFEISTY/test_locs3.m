%%%%!! RUN CLIMATOL FOR 3 LOCATIONS: EBS, HOT, PUP
function test_locs3()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Make core parameters/constants
param = make_params_testcase();

%! Idealized bathymetry & forcing
load(['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/fish-offline/MatFEISTY/input_files/',...
    'feisty_input_climatol_daily_locs3.mat'],'GRD','ESM');
param.NX = length(GRD.Z);
param.ID = 1:param.NX;
NX = param.NX;
ID = 1:param.NX;

%! How long to run the model
YEARS = 200;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
exper = 'locs3_';
[fname,simname] = sub_fname_testcase_exper(param,exper);

%! Initialize
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish_spin(ID,DAYS);

%! Storage variables
biom      = NaN*ones(DAYS,NX,9);
T_hab     = NaN*ones(DAYS,NX,9);
ingest    = NaN*ones(DAYS,NX,9);
pred_flux = NaN*ones(DAYS,NX,9);
pred_rate = NaN*ones(DAYS,NX,9);
metab     = NaN*ones(DAYS,NX,9);
mort      = NaN*ones(DAYS,NX,9);
energy    = NaN*ones(DAYS,NX,9);
growth    = NaN*ones(DAYS,NX,9);
repro     = NaN*ones(DAYS,NX,9);
recruit   = NaN*ones(DAYS,NX,9);
frate     = NaN*ones(DAYS,NX,9);

biomass             = NaN*ones(12*YEARS,NX,9);
T_habitat           = NaN*ones(12*YEARS,NX,9);
ingestion_rate      = NaN*ones(12*YEARS,NX,9);
predation_flux      = NaN*ones(12*YEARS,NX,9);
predation_rate      = NaN*ones(12*YEARS,NX,9);
metabolism_rate     = NaN*ones(12*YEARS,NX,9);
mortality_rate      = NaN*ones(12*YEARS,NX,9);
energy_avail_rate   = NaN*ones(12*YEARS,NX,9);
growth_rate         = NaN*ones(12*YEARS,NX,9);
reproduction_rate   = NaN*ones(12*YEARS,NX,9);
recruitment_flux    = NaN*ones(12*YEARS,NX,9);
fish_catch_rate     = NaN*ones(12*YEARS,NX,9);

%% %%%%%%%%%%%%%%%%%%%% Run the Model
MNT = 0;
%! Run model with no fishing
for YR = 1:YEARS % years
    num2str(YR)
        for DAY = 1:param.DT:DAYS % days
        
        %%%! Future time step
        DY = int64(ceil(DAY));
        %[num2str(YR),' , ', num2str(mod(DY,365))]
        
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio_1meso(DY,ESM,GRD,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);
        
        %! Store
        biom(DY,:,:) = [Sml_f.bio, Sml_p.bio, Sml_d.bio,...
                          Med_f.bio, Med_p.bio, Med_d.bio,...
                          Lrg_p.bio, Lrg_d.bio, BENT.mass];
        T_hab(DY,:,1:8) = [Sml_f.thab, Sml_p.thab, Sml_d.thab,...
                               Med_f.thab, Med_p.thab, Med_d.thab,...
                               Lrg_p.thab, Lrg_d.thab];
        ingest(DY,:,1:8) = [Sml_f.I, Sml_p.I, Sml_d.I,...
                                   Med_f.I, Med_p.I, Med_d.I,...
                                   Lrg_p.I, Lrg_d.I];
        pred_flux(DY,:,1:8) = [Sml_f.die, Sml_p.die, Sml_d.die,...
                                   Med_f.die, Med_p.die, Med_d.die,...
                                   Lrg_p.die, Lrg_d.die];
        pred_rate(DY,:,1:8) = [Sml_f.pred, Sml_p.pred, Sml_d.pred,...
                                   Med_f.pred, Med_p.pred, Med_d.pred,...
                                   Lrg_p.pred, Lrg_d.pred];
        metab(DY,:,1:8) = [Sml_f.met, Sml_p.met, Sml_d.met,...
                                   Med_f.met, Med_p.met, Med_d.met,...
                                   Lrg_p.met, Lrg_d.met];
        mort(DY,:,1:8) = [Sml_f.nmort, Sml_p.nmort, Sml_d.nmort,...
                                   Med_f.nmort, Med_p.nmort, Med_d.nmort,...
                                   Lrg_p.nmort, Lrg_d.nmort];
        energy(DY,:,1:8) = [Sml_f.nu, Sml_p.nu, Sml_d.nu,...
                                      Med_f.nu, Med_p.nu, Med_d.nu,...
                                      Lrg_p.nu, Lrg_d.nu];
        growth(DY,:,1:8) = [Sml_f.gamma, Sml_p.gamma, Sml_d.gamma,...
                                 Med_f.gamma, Med_p.gamma, Med_d.gamma,...
                                 Lrg_p.gamma, Lrg_d.gamma];
        repro(DY,:,1:8) = [Sml_f.rep, Sml_p.rep, Sml_d.rep,...
                                      Med_f.rep, Med_p.rep, Med_d.rep,...
                                      Lrg_p.rep, Lrg_d.rep];
        recruit(DY,:,1:8) = [Sml_f.rec, Sml_p.rec, Sml_d.rec,...
                                      Med_f.rec, Med_p.rec, Med_d.rec,...
                                      Lrg_p.rec, Lrg_d.rec];
        frate(DY,:,1:8) = [Sml_f.fmort, Sml_p.fmort, Sml_d.fmort,...
                                     Med_f.fmort, Med_p.fmort, Med_d.fmort,...
                                     Lrg_p.fmort, Lrg_d.fmort];
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH); % end of the month
    for i = 1:12
        MNT = MNT+1; % Update monthly ticker
        biomass(MNT,:,:) = mean(biom(a(i):b(i),:,:),1);
        T_habitat(MNT,:,:) = mean(T_hab(a(i):b(i),:,:),1);
        ingestion_rate(MNT,:,:) = mean(ingest(a(i):b(i),:,:),1);
        predation_flux(MNT,:,:) = mean(pred_flux(a(i):b(i),:,:),1);
        predation_rate(MNT,:,:) = mean(pred_rate(a(i):b(i),:,:),1);
        metabolism_rate(MNT,:,:) = mean(metab(a(i):b(i),:,:),1);
        mortality_rate(MNT,:,:) = mean(mort(a(i):b(i),:,:),1);
        energy_avail_rate(MNT,:,:) = mean(energy(a(i):b(i),:,:),1);
        growth_rate(MNT,:,:) = mean(growth(a(i):b(i),:,:),1);
        reproduction_rate(MNT,:,:) = mean(repro(a(i):b(i),:,:),1);
        recruitment_flux(MNT,:,:) = mean(recruit(a(i):b(i),:,:),1);
        fish_catch_rate(MNT,:,:) = mean(frate(a(i):b(i),:,:),1);
    end
    
end %Years

%%% Save
save([fname '.mat'],...
    'biomass','T_habitat','ingestion_rate','predation_flux','predation_rate',...
    'metabolism_rate','mortality_rate','energy_avail_rate','growth_rate',...
    'reproduction_rate','recruitment_flux','fish_catch_rate')

end