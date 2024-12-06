%%%%!! RUN HISTORIC FOR ALL LOCATIONS
function FOSI_cesm_rec_gamOut()

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
S_Sml_f = zeros(NX,DAYS);
S_Sml_p = zeros(NX,DAYS);
S_Sml_d = zeros(NX,DAYS);
S_Med_f = zeros(NX,DAYS);
S_Med_p = zeros(NX,DAYS);
S_Med_d = zeros(NX,DAYS);
S_Lrg_p = zeros(NX,DAYS);
S_Lrg_d = zeros(NX,DAYS);

G_Sml_f = zeros(NX,DAYS);
G_Sml_p = zeros(NX,DAYS);
G_Sml_d = zeros(NX,DAYS);
G_Med_f = zeros(NX,DAYS);
G_Med_p = zeros(NX,DAYS);
G_Med_d = zeros(NX,DAYS);
G_Lrg_p = zeros(NX,DAYS);
G_Lrg_d = zeros(NX,DAYS);

%! Initialize
init_sim = [exper simname];
load(['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/',simname '/FOSI/Last_mo_spin_' init_sim '.mat']);
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_sml_f = [fname,'_rec_gamma_sml_f.nc'];
file_sml_p = [fname,'_rec_gamma_sml_p.nc'];
file_sml_d = [fname,'_rec_gamma_sml_d.nc'];
file_med_p = [fname,'_rec_gamma_med_p.nc'];
file_med_d = [fname,'_rec_gamma_med_d.nc'];
file_med_f = [fname,'_rec_gamma_med_f.nc'];
file_lrg_p = [fname,'_rec_gamma_lrg_p.nc'];
file_lrg_d = [fname,'_rec_gamma_lrg_d.nc'];

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
vidrecSF    = netcdf.defVar(ncidSF,'rec','double',[xy_dim,time_dim]);
vidgamSF    = netcdf.defVar(ncidSF,'gamma','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidrecSP    = netcdf.defVar(ncidSP,'rec','double',[xy_dim,time_dim]);
vidgamSP    = netcdf.defVar(ncidSP,'gamma','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidrecSD    = netcdf.defVar(ncidSD,'rec','double',[xy_dim,time_dim]);
vidgamSD    = netcdf.defVar(ncidSD,'gamma','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidrecMP    = netcdf.defVar(ncidMP,'rec','double',[xy_dim,time_dim]);
vidgamMP    = netcdf.defVar(ncidMP,'gamma','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidrecMD    = netcdf.defVar(ncidMD,'rec','double',[xy_dim,time_dim]);
vidgamMD    = netcdf.defVar(ncidMD,'gamma','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidrecMF    = netcdf.defVar(ncidMF,'rec','double',[xy_dim,time_dim]);
vidgamMF    = netcdf.defVar(ncidMF,'gamma','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidrecLP    = netcdf.defVar(ncidLP,'rec','double',[xy_dim,time_dim]);
vidgamLP    = netcdf.defVar(ncidLP,'gamma','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidrecLD    = netcdf.defVar(ncidLD,'rec','double',[xy_dim,time_dim]);
vidgamLD    = netcdf.defVar(ncidLD,'gamma','double',[xy_dim,time_dim]);
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
        S_Sml_f(:,DY) = Sml_f.rec;
        S_Sml_p(:,DY) = Sml_p.rec;
        S_Sml_d(:,DY) = Sml_d.rec;
        S_Med_f(:,DY) = Med_f.rec;
        S_Med_p(:,DY) = Med_p.rec;
        S_Med_d(:,DY) = Med_d.rec;
        S_Lrg_p(:,DY) = Lrg_p.rec;
        S_Lrg_d(:,DY) = Lrg_d.rec;
        
        G_Sml_f(:,DY) = Sml_f.gamma;
        G_Sml_p(:,DY) = Sml_p.gamma;
        G_Sml_d(:,DY) = Sml_d.gamma;
        G_Med_f(:,DY) = Med_f.gamma;
        G_Med_p(:,DY) = Med_p.gamma;
        G_Med_d(:,DY) = Med_d.gamma;
        G_Lrg_p(:,DY) = Lrg_p.gamma;
        G_Lrg_d(:,DY) = Lrg_d.gamma;
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH);    % end of the month
    for i = 1:12
        MNT = MNT+1;     % Update monthly ticker
        
        %! Put vars of netcdf file
        netcdf.putVar(ncidMF,vidrecMF,[0 MNT-1],[NX 1],mean(S_Med_f(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidrecLP,[0 MNT-1],[NX 1],mean(S_Lrg_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidrecLD,[0 MNT-1],[NX 1],mean(S_Lrg_d(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidgamMF,[0 MNT-1],[NX 1],mean(G_Med_f(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidgamLP,[0 MNT-1],[NX 1],mean(G_Lrg_p(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidgamLD,[0 MNT-1],[NX 1],mean(G_Lrg_d(:,a(i):b(i)),2));
       
    end %Monthly mean
    
end %Years

%%
%! Close save
netcdf.close(ncidMF);
netcdf.close(ncidLP);
netcdf.close(ncidLD);

end
