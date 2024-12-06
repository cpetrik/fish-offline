%%%%!! RUN HISTORIC FOR ALL LOCATIONS
function FOSI_cesm_death_vars()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Make core parameters/constants 
param = make_parameters_1meso(); 

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
S_Sml_f_mort = zeros(NX,DAYS);
S_Sml_p_mort = zeros(NX,DAYS);
S_Sml_d_mort = zeros(NX,DAYS);
S_Med_f_mort = zeros(NX,DAYS);
S_Med_p_mort = zeros(NX,DAYS);
S_Med_d_mort = zeros(NX,DAYS);
S_Lrg_p_mort = zeros(NX,DAYS);
S_Lrg_d_mort = zeros(NX,DAYS);

S_Sml_f_die = zeros(NX,DAYS);
S_Sml_p_die = zeros(NX,DAYS);
S_Sml_d_die = zeros(NX,DAYS);
S_Med_f_die = zeros(NX,DAYS);
S_Med_p_die = zeros(NX,DAYS);
S_Med_d_die = zeros(NX,DAYS);
S_Lrg_p_die = zeros(NX,DAYS);
S_Lrg_d_die = zeros(NX,DAYS);

%! Initialize
init_sim = [exper simname];
load(['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/',simname '/FOSI/Last_mo_spin_' init_sim '.mat']);
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_sml_f = [fname,'_die_nmort_sml_f.nc'];
file_sml_p = [fname,'_die_nmort_sml_p.nc'];
file_sml_d = [fname,'_die_nmort_sml_d.nc'];
file_med_f = [fname,'_die_nmort_med_f.nc'];
file_med_p = [fname,'_die_nmort_med_p.nc'];
file_med_d = [fname,'_die_nmort_med_d.nc'];
file_lrg_p = [fname,'_die_nmort_lrg_p.nc'];
file_lrg_d = [fname,'_die_nmort_lrg_d.nc'];

ncidSF = netcdf.create(file_sml_f,'NC_WRITE');
ncidSP = netcdf.create(file_sml_p,'NC_WRITE');
ncidSD = netcdf.create(file_sml_d,'NC_WRITE');
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
xy_dim      = netcdf.defDim(ncidSF,'nid',NX);
time_dim    = netcdf.defDim(ncidSF,'ntime',nt);
vidmortSF     = netcdf.defVar(ncidSF,'mort','double',[xy_dim,time_dim]);
viddieSF    = netcdf.defVar(ncidSF,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidmortSP     = netcdf.defVar(ncidSP,'mort','double',[xy_dim,time_dim]);
viddieSP    = netcdf.defVar(ncidSP,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidmortSD     = netcdf.defVar(ncidSD,'mort','double',[xy_dim,time_dim]);
viddieSD    = netcdf.defVar(ncidSD,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidmortMF     = netcdf.defVar(ncidMF,'mort','double',[xy_dim,time_dim]);
viddieMF    = netcdf.defVar(ncidMF,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidmortMP     = netcdf.defVar(ncidMP,'mort','double',[xy_dim,time_dim]);
viddieMP    = netcdf.defVar(ncidMP,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidmortMD     = netcdf.defVar(ncidMD,'mort','double',[xy_dim,time_dim]);
viddieMD    = netcdf.defVar(ncidMD,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidmortLP     = netcdf.defVar(ncidLP,'mort','double',[xy_dim,time_dim]);
viddieLP    = netcdf.defVar(ncidLP,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidmortLD     = netcdf.defVar(ncidLD,'mort','double',[xy_dim,time_dim]);
viddieLD    = netcdf.defVar(ncidLD,'die','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLD);

%% %%%%%%%%%%%%%%%%%%%% Run the Model
MNT = 0;
%! Run model with no fishing
for YR = 1:YEARS % years
    %! Load a year's ESM data
    ti = num2str(YR)
    load(['/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/Data_cesm_fosi_v7_daily_',ti,'.mat'],'ESM');
    
    for DAY = 1:param.DT:DAYS % days
        
        %%%! Future time step
        DY = int64(ceil(DAY));
%         [num2str(YR),' , ', num2str(mod(DY,365))]
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
            sub_futbio_1meso(DY,ESM,GRD,Sml_f,Sml_p,Sml_d,...
            Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);
        
        %! Store
        S_Sml_f_die(:,DY) = Sml_f.die;
        S_Sml_p_die(:,DY) = Sml_p.die;
        S_Sml_d_die(:,DY) = Sml_d.die;
        S_Med_f_die(:,DY) = Med_f.die;
        S_Med_p_die(:,DY) = Med_p.die;
        S_Med_d_die(:,DY) = Med_d.die;
        S_Lrg_p_die(:,DY) = Lrg_p.die;
        S_Lrg_d_die(:,DY) = Lrg_d.die;

        S_Sml_f_mort(:,DY) = Sml_f.nmort;
        S_Sml_p_mort(:,DY) = Sml_p.nmort;
        S_Sml_d_mort(:,DY) = Sml_d.nmort;
        S_Med_f_mort(:,DY) = Med_f.nmort;
        S_Med_p_mort(:,DY) = Med_p.nmort;
        S_Med_d_mort(:,DY) = Med_d.nmort;
        S_Lrg_p_mort(:,DY) = Lrg_p.nmort;
        S_Lrg_d_mort(:,DY) = Lrg_d.nmort;

    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH);    % end of the month
    for i = 1:12
        MNT = MNT+1;     % Update monthly ticker
        
        %! Put vars of netcdf file
        netcdf.putVar(ncidSF,vidmortSF,[0 MNT-1],[NX 1],mean(S_Sml_f_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,vidmortSP,[0 MNT-1],[NX 1],mean(S_Sml_p_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,vidmortSD,[0 MNT-1],[NX 1],mean(S_Sml_d_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidmortMF,[0 MNT-1],[NX 1],mean(S_Med_f_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidmortMP,[0 MNT-1],[NX 1],mean(S_Med_p_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidmortMD,[0 MNT-1],[NX 1],mean(S_Med_d_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidmortLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_mort(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidmortLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_mort(:,a(i):b(i)),2));
        
        netcdf.putVar(ncidSF,viddieSF,[0 MNT-1],[NX 1],mean(S_Sml_f_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,viddieSP,[0 MNT-1],[NX 1],mean(S_Sml_p_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,viddieSD,[0 MNT-1],[NX 1],mean(S_Sml_d_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,viddieMF,[0 MNT-1],[NX 1],mean(S_Med_f_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,viddieMP,[0 MNT-1],[NX 1],mean(S_Med_p_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,viddieMD,[0 MNT-1],[NX 1],mean(S_Med_d_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,viddieLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_die(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,viddieLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_die(:,a(i):b(i)),2));
       
    end %Monthly mean
    
end %Years

%%
%! Close save
netcdf.close(ncidSF);
netcdf.close(ncidSP);
netcdf.close(ncidSD);
netcdf.close(ncidMF);
netcdf.close(ncidMP);
netcdf.close(ncidMD);
netcdf.close(ncidLP);
netcdf.close(ncidLD);

end
