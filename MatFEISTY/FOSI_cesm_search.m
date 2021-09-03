%%%%!! RUN HISTORIC FOR ALL LOCATIONS
function FOSI_cesm_search()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Grid
load('/Volumes/MIP/GCM_DATA/CESM/FOSI/Data_grid_POP_gx1v6.mat','GRD');
param.NX = GRD.N;
param.ID = 1:param.NX;
NX = param.NX;
ID = 1:param.NX;

%! How long to run the model
YEARS = 68;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

% % Orig run for init file
%! Set fishing rate
paramI.frate = 0.3;
paramI.dfrate = paramI.frate/365.0;

%! Make core parameters/constants
paramI = make_parameters_1meso(paramI);

%! Create a directory for output
experI = 'v12_';
[fnameI,simnameI] = sub_fname_cesm_fosi_exper(paramI,experI);
init_sim = [experI simnameI];

%% Loop over search parameter
mz = 0.1:0.1:0.9;
for n=1:length(mz)
    %! Make core parameters/constants with new permutation
    param.frate = 0.3;
    param.dfrate = param.frate/365.0;
    mzp = mz(n);
    param = make_parameters_1meso_mzpref(param,mzp);
    
    %! Create a directory for output
    tmzp = num2str(1000+int64(100*mzp));
    exper = ['sMZ' tmzp(2:end) '_'];
    [fname,simname] = sub_fname_cesm_fosi_exper(param,exper);
    
    %! Initialize
    load(['/Volumes/MIP/NC/CESM_MAPP/',simnameI '/Last_mo_spin_' init_sim '.mat']);
    BENT.mass = BENT.bio;
    [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);
    
    %! Storage variables
    % Dims
    nt = 12*YEARS;
    
    S_Bent_bio = zeros(NX,DAYS);
    S_Mzoo_frac = zeros(NX,DAYS);
    
    S_Sml_f = zeros(NX,DAYS);
    S_Sml_p = zeros(NX,DAYS);
    S_Sml_d = zeros(NX,DAYS);
    S_Med_f = zeros(NX,DAYS);
    S_Med_p = zeros(NX,DAYS);
    S_Med_d = zeros(NX,DAYS);
    S_Lrg_p = zeros(NX,DAYS);
    S_Lrg_d = zeros(NX,DAYS);
    
    S_Med_f_fish = zeros(NX,DAYS);
    S_Med_p_fish = zeros(NX,DAYS);
    S_Med_d_fish = zeros(NX,DAYS);
    S_Lrg_p_fish = zeros(NX,DAYS);
    S_Lrg_d_fish = zeros(NX,DAYS);
    
    Spinup_Sml_f_bio = NaN*ones(NX,nt);
    Spinup_Sml_p_bio = NaN*ones(NX,nt);
    Spinup_Sml_d_bio = NaN*ones(NX,nt);
    Spinup_Med_f_bio = NaN*ones(NX,nt);
    Spinup_Med_p_bio = NaN*ones(NX,nt);
    Spinup_Med_d_bio = NaN*ones(NX,nt);
    Spinup_Lrg_p_bio = NaN*ones(NX,nt);
    Spinup_Lrg_d_bio = NaN*ones(NX,nt);
    Spinup_Bent_bio = NaN*ones(NX,nt);
    Spinup_Mzoo_frac = NaN*ones(NX,nt);
    
    Spinup_Med_f_yield = NaN*ones(NX,nt);
    Spinup_Med_p_yield = NaN*ones(NX,nt);
    Spinup_Med_d_yield = NaN*ones(NX,nt);
    Spinup_Lrg_p_yield = NaN*ones(NX,nt);
    Spinup_Lrg_d_yield = NaN*ones(NX,nt);
    
    %% %%%%%%%%%%%%%%%%%%%% Run the Model
    MNT = 0;
    %! Run model with no fishing
    for YR = 1:YEARS % years
        %! Load a year's ESM data
        ti = num2str(YR);
        [num2str(n) ',' num2str(ti)]
        load(['/Volumes/MIP/GCM_DATA/CESM/FOSI/Data_cesm_fosi_loss_v5_daily_',ti,'.mat'],'ESM');
        
        for DAY = 1:param.DT:DAYS % days
            
            %%%! Future time step
            DY = int64(ceil(DAY));
            %         [num2str(YR),' , ', num2str(mod(DY,365))]
            [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
                sub_futbio_1meso_mzprefM(DY,ESM,GRD,Sml_f,Sml_p,Sml_d,...
                Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);
            
            %! Store
            S_Bent_bio(:,DY) = BENT.mass;
            S_Mzoo_frac(:,DY) = ENVR.fZm;
            
            S_Sml_f(:,DY) = Sml_f.bio;
            S_Sml_p(:,DY) = Sml_p.bio;
            S_Sml_d(:,DY) = Sml_d.bio;
            S_Med_f(:,DY) = Med_f.bio;
            S_Med_p(:,DY) = Med_p.bio;
            S_Med_d(:,DY) = Med_d.bio;
            S_Lrg_p(:,DY) = Lrg_p.bio;
            S_Lrg_d(:,DY) = Lrg_d.bio;
            
            S_Med_f_fish(:,DY) = Med_f.caught;
            S_Med_p_fish(:,DY) = Med_p.caught;
            S_Med_d_fish(:,DY) = Med_d.caught;
            S_Lrg_p_fish(:,DY) = Lrg_p.caught;
            S_Lrg_d_fish(:,DY) = Lrg_d.caught;
            
        end %Days
        
        %! Calculate monthly means and save
        aa = (cumsum(MNTH)+1);
        a = [1,aa(1:end-1)]; % start of the month
        b = cumsum(MNTH); % end of the month
        for i = 1:12
            MNT = MNT+1;     % Update monthly ticker
            
            Spinup_Bent_bio(:,MNT) = mean(S_Bent_bio(:,a(i):b(i)),2);
            Spinup_Mzoo_frac(:,MNT) = mean(S_Mzoo_frac(:,a(i):b(i)),2);
            Spinup_Sml_f_bio(:,MNT) = mean(S_Sml_f(:,a(i):b(i)),2);
            Spinup_Sml_p_bio(:,MNT) = mean(S_Sml_p(:,a(i):b(i)),2);
            Spinup_Sml_d_bio(:,MNT) = mean(S_Sml_d(:,a(i):b(i)),2);
            Spinup_Med_f_bio(:,MNT) = mean(S_Med_f(:,a(i):b(i)),2);
            Spinup_Med_p_bio(:,MNT) = mean(S_Med_p(:,a(i):b(i)),2);
            Spinup_Med_d_bio(:,MNT) = mean(S_Med_d(:,a(i):b(i)),2);
            Spinup_Lrg_p_bio(:,MNT) = mean(S_Lrg_p(:,a(i):b(i)),2);
            Spinup_Lrg_d_bio(:,MNT) = mean(S_Lrg_d(:,a(i):b(i)),2);
            
            Spinup_Med_f_yield(:,MNT) = mean(S_Med_f_fish(:,a(i):b(i)),2);
            Spinup_Med_p_yield(:,MNT) = mean(S_Med_p_fish(:,a(i):b(i)),2);
            Spinup_Med_d_yield(:,MNT) = mean(S_Med_d_fish(:,a(i):b(i)),2);
            Spinup_Lrg_p_yield(:,MNT) = mean(S_Lrg_p_fish(:,a(i):b(i)),2);
            Spinup_Lrg_d_yield(:,MNT) = mean(S_Lrg_d_fish(:,a(i):b(i)),2);
            
        end %Monthly mean
        
    end %Years
    
    %% Take mean of all years and save
    FOSI_Bent_bio = mean(Spinup_Bent_bio,2);
    FOSI_Mzoo_frac = mean(Spinup_Mzoo_frac,2);
    FOSI_Sml_f_bio = mean(Spinup_Sml_f_bio,2);
    FOSI_Sml_p_bio = mean(Spinup_Sml_p_bio,2);
    FOSI_Sml_d_bio = mean(Spinup_Sml_d_bio,2);
    FOSI_Med_f_bio = mean(Spinup_Med_f_bio,2);
    FOSI_Med_p_bio = mean(Spinup_Med_p_bio,2);
    FOSI_Med_d_bio = mean(Spinup_Med_d_bio,2);
    FOSI_Lrg_p_bio = mean(Spinup_Lrg_p_bio,2);
    FOSI_Lrg_d_bio = mean(Spinup_Lrg_d_bio,2);
    
    FOSI_Med_f_yield = mean(Spinup_Med_f_yield,2);
    FOSI_Med_p_yield = mean(Spinup_Med_p_yield,2);
    FOSI_Med_d_yield = mean(Spinup_Med_d_yield,2);
    FOSI_Lrg_p_yield = mean(Spinup_Lrg_p_yield,2);
    FOSI_Lrg_d_yield = mean(Spinup_Lrg_d_yield,2);
    
    %%% Save
    save([fname,'_means.mat'],...
        'FOSI_Sml_f_bio','FOSI_Sml_p_bio','FOSI_Sml_d_bio','FOSI_Med_f_bio',...
        'FOSI_Med_p_bio','FOSI_Med_d_bio','FOSI_Lrg_p_bio','FOSI_Lrg_d_bio',...
        'FOSI_Bent_bio','FOSI_Mzoo_frac','FOSI_Med_f_yield','FOSI_Med_p_yield',...
        'FOSI_Med_d_yield','FOSI_Lrg_p_yield','FOSI_Lrg_d_yield')
    
    
end
