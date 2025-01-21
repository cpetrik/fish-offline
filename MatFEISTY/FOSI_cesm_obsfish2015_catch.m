%%%%!! RUN FOSI FOR ALL LOCATIONS
function FOSI_cesm_obsfish2015_catch()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Make core parameters/constants 
param = make_parameters_1meso_obsfish();

% Assessment method of estimating
spath = '/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/assessment/';
load([spath 'FOSI_POP_gx1v6_noSeas_fmort_ID_annual_1948_2015_tempSc_assessment.mat'],...
    'fmD','fmP','fmF');

% Set fishing rate as 1st year for fname
param.frate = nan;
param.frateF = fmF(:,1);
param.frateP = fmP(:,1);
param.frateD = fmD(:,1);
param.dfrateF = param.frateF/365.0;
param.dfrateP = param.frateP/365.0;
param.dfrateD = param.frateD/365.0;

%! Grid
load('/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/Data_grid_POP_gx1v6_noSeas.mat','GRD');
param.NX = GRD.N;
param.ID = 1:param.NX;
NX = param.NX;
ID = 1:param.NX;

%! How long to run the model
YEARS = 68;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];

%! Create a directory for output
exper = 'v15_';
[fname,simname] = sub_fname_cesm_fosi_exper(param,exper);

%! Storage variables
S_Med_f = zeros(NX,DAYS);
S_Med_p = zeros(NX,DAYS);
S_Med_d = zeros(NX,DAYS);
S_Lrg_p = zeros(NX,DAYS);
S_Lrg_d = zeros(NX,DAYS);

N_Med_f = zeros(NX,DAYS);
N_Lrg_p = zeros(NX,DAYS);
N_Lrg_d = zeros(NX,DAYS);

%! Initialize
init_sim = [exper 'obsfish_' simname];
load(['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/',simname '/FOSI/Last_mo_spin_' init_sim '.mat']);
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_med_f = [fname,'2015_catch_nu_med_f.nc'];
file_med_p = [fname,'2015_catch_nu_med_p.nc'];
file_med_d = [fname,'2015_catch_nu_med_d.nc'];
file_lrg_p = [fname,'2015_catch_nu_lrg_p.nc'];
file_lrg_d = [fname,'2015_catch_nu_lrg_d.nc'];

ncidMF = netcdf.create(file_med_f,'NC_WRITE');
ncidMP = netcdf.create(file_med_p,'NC_WRITE');
ncidMD = netcdf.create(file_med_d,'NC_WRITE');
ncidLP = netcdf.create(file_lrg_p,'NC_WRITE');
ncidLD = netcdf.create(file_lrg_d,'NC_WRITE');

%! Dims of netcdf file
nt = 12*YEARS;
netcdf.setDefaultFormat('NC_FORMAT_64BIT');

%% ! Def vars of netcdf file
['Defining netcdfs, takes ~5 minutes ... ']
xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidnuMF    = netcdf.defVar(ncidMF,'nu','double',[xy_dim,time_dim]);
vidbioMF    = netcdf.defVar(ncidMF,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidbioMP    = netcdf.defVar(ncidMP,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidbioMD    = netcdf.defVar(ncidMD,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidnuLP    = netcdf.defVar(ncidLP,'nu','double',[xy_dim,time_dim]);
vidbioLP    = netcdf.defVar(ncidLP,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidnuLD    = netcdf.defVar(ncidLD,'nu','double',[xy_dim,time_dim]);
vidbioLD    = netcdf.defVar(ncidLD,'yield','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);

%% %%%%%%%%%%%%%%%%%%%% Run the Model
MNT = 0;
%! Run model with no fishing
for YR = 1:YEARS % years
    %! Load a year's ESM data
    ti = num2str(YR)
    load(['/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/Data_cesm_fosi_v7_daily_',ti,'.mat'],'ESM');

    param.frateF = fmF(:,YR);
    param.frateP = fmP(:,YR);
    param.frateD = fmD(:,YR);
    param.dfrateF = param.frateF/365.0;
    param.dfrateP = param.frateP/365.0;
    param.dfrateD = param.frateD/365.0;
    
    for DAY = 1:param.DT:DAYS % days
        
        %%%! Future time step
        DY = int64(ceil(DAY));
%         [num2str(YR),' , ', num2str(mod(DY,365))]
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio_1meso_obsfish(DY,ESM,GRD,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);
        
        %! Store
        S_Med_f(:,DY) = Med_f.caught;
        S_Med_p(:,DY) = Med_p.caught;
        S_Med_d(:,DY) = Med_d.caught;
        S_Lrg_p(:,DY) = Lrg_p.caught;
        S_Lrg_d(:,DY) = Lrg_d.caught;

        N_Med_f(:,DY) = Med_f.nu;
        N_Lrg_p(:,DY) = Lrg_p.nu;
        N_Lrg_d(:,DY) = Lrg_d.nu;
        
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH);    % end of the month
    for i = 1:12
        MNT = MNT+1;     % Update monthly ticker
        
        %! Put vars of netcdf file
        netcdf.putVar(ncidMF,vidbioMF,[0 MNT-1],[NX 1],mean(S_Med_f(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidbioMP,[0 MNT-1],[NX 1],mean(S_Med_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidbioMD,[0 MNT-1],[NX 1],mean(S_Med_d(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidbioLP,[0 MNT-1],[NX 1],mean(S_Lrg_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidbioLD,[0 MNT-1],[NX 1],mean(S_Lrg_d(:,a(i):b(i)),2));

        netcdf.putVar(ncidMF,vidnuMF,[0 MNT-1],[NX 1],mean(N_Med_f(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidnuLP,[0 MNT-1],[NX 1],mean(N_Lrg_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidnuLD,[0 MNT-1],[NX 1],mean(N_Lrg_d(:,a(i):b(i)),2));
       
    end %Monthly mean
    
end %Years

%%
%! Close save
netcdf.close(ncidMF);
netcdf.close(ncidMP);
netcdf.close(ncidMD);
netcdf.close(ncidLP);
netcdf.close(ncidLD);

end
