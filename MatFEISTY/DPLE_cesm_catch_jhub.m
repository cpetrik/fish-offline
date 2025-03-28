%%%%!! RUN HISTORIC FOR ALL LOCATIONS
function DPLE_cesm_catch_jhub()

%%%%%%%%%%%%%%% Initialize Model Variables
%! Make core parameters/constants
param = make_parameters_1meso();

%! Grid
spath = '/glade/scratch/cpetrik/fish-offline/dailies/';
Cdir = '/glade/u/home/cpetrik/fish-offline/MatFEISTY/input_files/';
load([Cdir 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');
param.NX = GRD.N;
param.ID = 1:param.NX;
NX = param.NX;
ID = 1:param.NX;

%! How long to run the model
YEARS = 10;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];
StartYr = 1954:2017; %will loop over
submem = 1:40;

for sim = 1:length(StartYr)
    for mem = 1:length(submem) %will loop over
        Member = submem(mem);
        
        %! Create a directory for output
        exper = ['v14_Y' num2str(StartYr(sim)) '_M' num2str(Member) '_' ];
        [fname,simname] = sub_fname_dple_exper_jhub(param,exper);
        
        %! Storage variables
        S_Med_f = zeros(NX,DAYS);
        S_Med_p = zeros(NX,DAYS);
        S_Med_d = zeros(NX,DAYS);
        S_Lrg_p = zeros(NX,DAYS);
        S_Lrg_d = zeros(NX,DAYS);
        
        %! Initialize
        % this will have to be the biomass from the last mo before the start year
        % from the FOSI
        init_sim = ['v14_All_fish03_' simname];
        load([Cdir,simname,'/Last_mo_FOSI_' num2str(StartYr) '_' init_sim '.mat']);
        BENT.mass = BENT.bio;
        [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);
        
        %%%%%%%%%%%%%%% Setup NetCDF save
        %! Setup netcdf path to store to
        file_med_f = [fname,'_catch_med_f.nc'];
        file_med_p = [fname,'_catch_med_p.nc'];
        file_med_d = [fname,'_catch_med_d.nc'];
        file_lrg_p = [fname,'_catch_lrg_p.nc'];
        file_lrg_d = [fname,'_catch_lrg_d.nc'];
        
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
        vidcatchMF    = netcdf.defVar(ncidMF,'yield','double',[xy_dim,time_dim]);
        netcdf.endDef(ncidMF);
        
        xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
        time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
        vidcatchMP    = netcdf.defVar(ncidMP,'yield','double',[xy_dim,time_dim]);
        netcdf.endDef(ncidMP);
        
        xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
        time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
        vidcatchMD    = netcdf.defVar(ncidMD,'yield','double',[xy_dim,time_dim]);
        netcdf.endDef(ncidMD);
        
        xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
        time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
        vidcatchLP    = netcdf.defVar(ncidLP,'yield','double',[xy_dim,time_dim]);
        netcdf.endDef(ncidLP);
        
        xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
        time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
        vidcatchLD    = netcdf.defVar(ncidLD,'yield','double',[xy_dim,time_dim]);
        netcdf.endDef(ncidLD);
        
        %% %%%%%%%%%%%%%%%%%%%% Run the Model
        MNT = 0;
        %! Run model with no fishing
        for YR = 1:YEARS % years
            %! Load a year's ESM data
            ti = num2str(YR);
            ['Y',num2str(YR),' M',num2str(Member)]
            load([spath 'Data_cesm_dple_daily_Y',num2str(StartYr(sim)),'_M',...
                num2str(Member),'_LY',ti,'.mat'],'ESM')
            
            for DAY = 1:param.DT:DAYS % days
                
                %%%! Future time step
                DY = int64(ceil(DAY));
                %         [num2str(YR),' , ', num2str(mod(DY,365))]
                [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
                    sub_futbio_1meso(DY,ESM,GRD,Sml_f,Sml_p,Sml_d,...
                    Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,param);
                
                %! Store
                S_Med_f(:,DY) = Med_f.caught;
                S_Med_p(:,DY) = Med_p.caught;
                S_Med_d(:,DY) = Med_d.caught;
                S_Lrg_p(:,DY) = Lrg_p.caught;
                S_Lrg_d(:,DY) = Lrg_d.caught;
                
            end %Days
            
            %! Calculate monthly means and save
            aa = (cumsum(MNTH)+1);
            a = [1,aa(1:end-1)]; % start of the month
            b = cumsum(MNTH);    % end of the month
            for i = 1:12
                MNT = MNT+1;     % Update monthly ticker
                
                %! Put vars of netcdf file
                netcdf.putVar(ncidMF,vidcatchMF,[0 MNT-1],[NX 1],mean(S_Med_f(:,a(i):b(i)),2));
                netcdf.putVar(ncidMP,vidcatchMP,[0 MNT-1],[NX 1],mean(S_Med_p(:,a(i):b(i)),2));
                netcdf.putVar(ncidMD,vidcatchMD,[0 MNT-1],[NX 1],mean(S_Med_d(:,a(i):b(i)),2));
                netcdf.putVar(ncidLP,vidcatchLP,[0 MNT-1],[NX 1],mean(S_Lrg_p(:,a(i):b(i)),2));
                netcdf.putVar(ncidLD,vidcatchLD,[0 MNT-1],[NX 1],mean(S_Lrg_d(:,a(i):b(i)),2));
                
            end %Monthly mean
            
        end %Years
        
        %%
        %! Close save
        netcdf.close(ncidMF);
        netcdf.close(ncidMP);
        netcdf.close(ncidMD);
        netcdf.close(ncidLP);
        netcdf.close(ncidLD);
        
    end %member loop
end %calendar year
end
