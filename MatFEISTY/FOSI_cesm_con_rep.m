%%%%!! RUN HISTORIC FOR ALL LOCATIONS
function FOSI_cesm_con_rep()

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
S_Sml_f_con = zeros(NX,DAYS);
S_Sml_p_con = zeros(NX,DAYS);
S_Sml_d_con = zeros(NX,DAYS);
S_Med_f_con = zeros(NX,DAYS);
S_Med_p_con = zeros(NX,DAYS);
S_Med_d_con = zeros(NX,DAYS);
S_Lrg_p_con = zeros(NX,DAYS);
S_Lrg_d_con = zeros(NX,DAYS);

S_Med_f_rep = zeros(NX,DAYS);
S_Lrg_p_rep = zeros(NX,DAYS);
S_Lrg_d_rep = zeros(NX,DAYS);

%! Initialize
init_sim = [exper simname];
load(['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/',simname '/FOSI/Last_mo_spin_' init_sim '.mat']);
BENT.mass = BENT.bio;
[Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

%%%%%%%%%%%%%%% Setup NetCDF save
%! Setup netcdf path to store to
file_sml_f = [fname,'_con_rep_sml_f.nc'];
file_sml_p = [fname,'_con_rep_sml_p.nc'];
file_sml_d = [fname,'_con_rep_sml_d.nc'];
file_med_f = [fname,'_con_rep_med_f.nc'];
file_med_p = [fname,'_con_rep_med_p.nc'];
file_med_d = [fname,'_con_rep_med_d.nc'];
file_lrg_p = [fname,'_con_rep_lrg_p.nc'];
file_lrg_d = [fname,'_con_rep_lrg_d.nc'];

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
vidconSF    = netcdf.defVar(ncidSF,'con','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSF);

xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
vidconSP    = netcdf.defVar(ncidSP,'con','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSP);

xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
vidconSD    = netcdf.defVar(ncidSD,'con','double',[xy_dim,time_dim]);
netcdf.endDef(ncidSD);

xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
vidconMF    = netcdf.defVar(ncidMF,'con','double',[xy_dim,time_dim]);
vidrepMF    = netcdf.defVar(ncidMF,'rep','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMF);

xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
vidconMP    = netcdf.defVar(ncidMP,'con','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMP);

xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
vidconMD    = netcdf.defVar(ncidMD,'con','double',[xy_dim,time_dim]);
netcdf.endDef(ncidMD);

xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
vidconLP    = netcdf.defVar(ncidLP,'con','double',[xy_dim,time_dim]);
vidrepLP    = netcdf.defVar(ncidLP,'rep','double',[xy_dim,time_dim]);
netcdf.endDef(ncidLP);

xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
vidconLD    = netcdf.defVar(ncidLD,'con','double',[xy_dim,time_dim]);
vidrepLD    = netcdf.defVar(ncidLD,'rep','double',[xy_dim,time_dim]);
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
        S_Sml_f_con(:,DY) = Sml_f.I;
        S_Sml_p_con(:,DY) = Sml_p.I;
        S_Sml_d_con(:,DY) = Sml_d.I;
        S_Med_f_con(:,DY) = Med_f.I;
        S_Med_p_con(:,DY) = Med_p.I;
        S_Med_d_con(:,DY) = Med_d.I;
        S_Lrg_p_con(:,DY) = Lrg_p.I;
        S_Lrg_d_con(:,DY) = Lrg_d.I;
        
        S_Med_f_rep(:,DY) = Med_f.rep;
        S_Lrg_p_rep(:,DY) = Lrg_p.rep;
        S_Lrg_d_rep(:,DY) = Lrg_d.rep;
        
    end %Days
    
    %! Calculate monthly means and save
    aa = (cumsum(MNTH)+1);
    a = [1,aa(1:end-1)]; % start of the month
    b = cumsum(MNTH);    % end of the month
    for i = 1:12
        MNT = MNT+1;     % Update monthly ticker
        
        %! Put vars of netcdf file
        netcdf.putVar(ncidMF,vidrepMF,[0 MNT-1],[NX 1],mean(S_Med_f_rep(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidrepLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_rep(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidrepLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_rep(:,a(i):b(i)),2));        

        netcdf.putVar(ncidSF,vidconSF,[0 MNT-1],[NX 1],mean(S_Sml_f_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidSP,vidconSP,[0 MNT-1],[NX 1],mean(S_Sml_p_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidSD,vidconSD,[0 MNT-1],[NX 1],mean(S_Sml_d_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidMF,vidconMF,[0 MNT-1],[NX 1],mean(S_Med_f_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidMP,vidconMP,[0 MNT-1],[NX 1],mean(S_Med_p_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidMD,vidconMD,[0 MNT-1],[NX 1],mean(S_Med_d_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidLP,vidconLP,[0 MNT-1],[NX 1],mean(S_Lrg_p_con(:,a(i):b(i)),2));
        netcdf.putVar(ncidLD,vidconLD,[0 MNT-1],[NX 1],mean(S_Lrg_d_con(:,a(i):b(i)),2));
       
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
