% Code for running DPLE simulations
% Order of operations:
% 1. Calc FOSI climatology
% 2. GET DPLE time info
% Year, Lead year, and Member loop starts here
% 3. Convert DPLE from anom to raw
% 4. Interpolate from monthly to daily
% 5. Run FEISTY DPLE
% Stop loop

clear 
close all

%%
Cdir = '/glade/u/home/cpetrik/fish-offline/MatFEISTY/input_files/';
fpath='/glade/scratch/kristenk/fish-offline/';
%fpath = '/glade/campaign/cesm/development/bgcwg/projects/DPLE-FEISTY/';
%spath = '/glade/scratch/cpetrik/fish-offline/dailies/';

% LOAD GRIDDATA
load([Cdir 'gridspec_POP_gx1v6_noSeas.mat'],'mask');
load([Cdir 'Data_grid_POP_gx1v6_noSeas.mat'],'GRD');

%% 1. Calc FOSI climatology ----------------------------
[clim_Tp,clim_Tb,clim_POC,clim_zooC,clim_loss] = fosi_clim_for_dple(Cdir);

%% 2. GET DPLE time info ----------------------------
ncdisp([fpath 'DPLE-FIESTY-forcing_zooC_150m.nc'])

% Pelagic temperature
ncid = netcdf.open([fpath 'DPLE-FIESTY-forcing_TEMP_150m.nc'],'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
for i = 1:(nvars-1)
    varname = netcdf.inqVar(ncid, i-1);
    eval([ varname ' = netcdf.getVar(ncid,i-1);']);
    eval([ varname '(' varname ' == 9.969209968386869e+36) = NaN;']);
end
[ni,nj] = size(TLONG);
whos M L Y

%% Times
% FOSI is 1948-2015
% DPLE is 1954-2017
% leadtimes = 1:24;
firstyear = 1954; %of DPLE initializations
lastyear  = 2015; %of DPLE initializations
startmonth = 11;  %of DPLE initializations

% TIME PERIOD FOR INTERPOLATION
mos = length(L);
mstart = 3:12:mos; %starts at Nov
mend = 14:12:mos;  %ends at Dec
mstart = [1 mstart];
mend = [2 mend];
nyrs = ceil(mos/12);
simy = 1:nyrs;

Tdays=1:365;

%% Set up FEISTY  --------------------------------------------
% Initialize Model Variables
%! Make core parameters/constants
param = make_parameters_1meso();

%! Grid
param.NX = GRD.N;
param.ID = 1:param.NX;
NX = param.NX;
ID = 1:param.NX;

%! How long to run the model
StartYr = 1954:2017;
DAYS = 365;
MNTH = [31,28,31,30,31,30,31,31,30,31,30,31];


%% START LOOP ===============================================================
% Loop over initialization year
for yr=1:length(Y) %Y(62)

    % Loop over member
    for mem=1:length(M)

        im = M(mem);
        iy = yr;

        %% More FEISTY set up
        %! Create a directory for output
        exper = ['v15_Y' num2str(StartYr(yr)) '_M' num2str(im) '_' ];
        [fname,simname] = sub_fname_dple_exper_jhub(param,exper);

        %! Storage variables
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

        P_Sml_f = zeros(NX,DAYS);
        P_Sml_p = zeros(NX,DAYS);
        P_Sml_d = zeros(NX,DAYS);
        P_Med_f = zeros(NX,DAYS);
        P_Med_p = zeros(NX,DAYS);
        P_Med_d = zeros(NX,DAYS);
        P_Lrg_p = zeros(NX,DAYS);
        P_Lrg_d = zeros(NX,DAYS);

        C_Med_f = zeros(NX,DAYS);
        C_Med_p = zeros(NX,DAYS);
        C_Med_d = zeros(NX,DAYS);
        C_Lrg_p = zeros(NX,DAYS);
        C_Lrg_d = zeros(NX,DAYS);

        %%%%%%%%%%%%%%% Setup NetCDF save
        %! Setup netcdf path to store to
        file_sml_f = [fname,'_sml_f.nc'];
        file_sml_p = [fname,'_sml_p.nc'];
        file_sml_d = [fname,'_sml_d.nc'];
        file_med_f = [fname,'_med_f.nc'];
        file_med_p = [fname,'_med_p.nc'];
        file_med_d = [fname,'_med_d.nc'];
        file_lrg_p = [fname,'_lrg_p.nc'];
        file_lrg_d = [fname,'_lrg_d.nc'];
        file_bent  = [fname,'_bent.nc'];
        file_mzoo  = [fname,'_mzoo.nc'];

        ncidSF = netcdf.create(file_sml_f,'NC_WRITE');
        ncidSP = netcdf.create(file_sml_p,'NC_WRITE');
        ncidSD = netcdf.create(file_sml_d,'NC_WRITE');
        ncidMF = netcdf.create(file_med_f,'NC_WRITE');
        ncidMP = netcdf.create(file_med_p,'NC_WRITE');
        ncidMD = netcdf.create(file_med_d,'NC_WRITE');
        ncidLP = netcdf.create(file_lrg_p,'NC_WRITE');
        ncidLD = netcdf.create(file_lrg_d,'NC_WRITE');
        ncidB  = netcdf.create(file_bent,'NC_WRITE');
        ncidMZ = netcdf.create(file_mzoo,'NC_WRITE');

        %! Dims of netcdf file
        nt = 12*YEARS;
        netcdf.setDefaultFormat('NC_FORMAT_64BIT');

        %% ! Def vars of netcdf file
        ['Defining netcdfs, takes ~5 minutes ... ']
        xy_dim      = netcdf.defDim(ncidSF,'nid',NX);
        time_dim    = netcdf.defDim(ncidSF,'ntime',nt);
        vidbioSF    = netcdf.defVar(ncidSF,'biomass','double',[xy_dim,time_dim]);
        vidprodSF   = netcdf.defVar(ncidSF,'prod','double',[xy_dim,time_dim]);
        netcdf.endDef(ncidSF);

        xy_dim      = netcdf.defDim(ncidSP,'nid',NX);
        time_dim    = netcdf.defDim(ncidSP,'ntime',nt);
        vidbioSP    = netcdf.defVar(ncidSP,'biomass','double',[xy_dim,time_dim]);
        vidprodSP   = netcdf.defVar(ncidSP,'prod','double',[xy_dim,time_dim]);
        netcdf.endDef(ncidSP);

        xy_dim      = netcdf.defDim(ncidSD,'nid',NX);
        time_dim    = netcdf.defDim(ncidSD,'ntime',nt);
        vidbioSD    = netcdf.defVar(ncidSD,'biomass','double',[xy_dim,time_dim]);
        vidprodSD   = netcdf.defVar(ncidSD,'prod','double',[xy_dim,time_dim]);
        netcdf.endDef(ncidSD);

        xy_dim      = netcdf.defDim(ncidMF,'nid',NX);
        time_dim    = netcdf.defDim(ncidMF,'ntime',nt);
        vidbioMF    = netcdf.defVar(ncidMF,'biomass','double',[xy_dim,time_dim]);
        vidprodMF   = netcdf.defVar(ncidMF,'prod','double',[xy_dim,time_dim]);
        vidcatchMF    = netcdf.defVar(ncidMF,'yield','double',[xy_dim,time_dim]);
        netcdf.endDef(ncidMF);

        xy_dim      = netcdf.defDim(ncidMP,'nid',NX);
        time_dim    = netcdf.defDim(ncidMP,'ntime',nt);
        vidbioMP    = netcdf.defVar(ncidMP,'biomass','double',[xy_dim,time_dim]);
        vidprodMP   = netcdf.defVar(ncidMP,'prod','double',[xy_dim,time_dim]);
        vidcatchMP    = netcdf.defVar(ncidMP,'yield','double',[xy_dim,time_dim]);
        netcdf.endDef(ncidMP);

        xy_dim      = netcdf.defDim(ncidMD,'nid',NX);
        time_dim    = netcdf.defDim(ncidMD,'ntime',nt);
        vidbioMD    = netcdf.defVar(ncidMD,'biomass','double',[xy_dim,time_dim]);
        vidprodMD   = netcdf.defVar(ncidMD,'prod','double',[xy_dim,time_dim]);
        vidcatchMD    = netcdf.defVar(ncidMD,'yield','double',[xy_dim,time_dim]);
        netcdf.endDef(ncidMD);

        xy_dim      = netcdf.defDim(ncidLP,'nid',NX);
        time_dim    = netcdf.defDim(ncidLP,'ntime',nt);
        vidbioLP    = netcdf.defVar(ncidLP,'biomass','double',[xy_dim,time_dim]);
        vidprodLP   = netcdf.defVar(ncidLP,'prod','double',[xy_dim,time_dim]);
        vidcatchLP    = netcdf.defVar(ncidLP,'yield','double',[xy_dim,time_dim]);
        netcdf.endDef(ncidLP);

        xy_dim      = netcdf.defDim(ncidLD,'nid',NX);
        time_dim    = netcdf.defDim(ncidLD,'ntime',nt);
        vidbioLD    = netcdf.defVar(ncidLD,'biomass','double',[xy_dim,time_dim]);
        vidprodLD   = netcdf.defVar(ncidLD,'prod','double',[xy_dim,time_dim]);
        vidcatchLD    = netcdf.defVar(ncidLD,'yield','double',[xy_dim,time_dim]);
        netcdf.endDef(ncidLD);

        xy_dim     = netcdf.defDim(ncidB,'nid',NX);
        time_dim   = netcdf.defDim(ncidB,'ntime',nt);
        vidbioB    = netcdf.defVar(ncidB,'biomass','double',[xy_dim,time_dim]);
        vidTB      = netcdf.defVar(ncidB,'time','double',time_dim);
        netcdf.endDef(ncidB);

        xy_dim      = netcdf.defDim(ncidMZ,'nid',NX);
        time_dim    = netcdf.defDim(ncidMZ,'ntime',nt);
        vidfracMZ   = netcdf.defVar(ncidMZ,'fraction','double',[xy_dim,time_dim]);
        vidTMZ      = netcdf.defVar(ncidMZ,'time','double',time_dim);
        netcdf.endDef(ncidMZ);


        %% 3. Convert DPLE from anom to raw ----------------------------
        [TEMP_150m,TEMP_bottom,POC_FLUX_IN_bottom,LzooC_150m,Lzoo_loss_150m] = dple_anom_to_raw(fpath,im,yr,ni,nj,L);


        %% Loop over 10 yrs + 2 months (run each year individually)
        MNT = 0;
        for y = 1:nyrs
            lyr = simy(y);
            ['Y',num2str(Y(iy)),'_M',num2str(im),'_LY',num2str(lyr)]


            %% 4. Interpolate from monthly to daily ----------------------------
            % 1ST MONTH IS NOV, ADJUSTED

            if y==1
                range = mstart(y):(mend(y)+1);
                Time=315:30:395;
            elseif y==nyrs
                range = (mstart(y)-1):mend(y);
                Time=-15:30:365;
            else
                range = (mstart(y)-1):(mend(y)+1);
                Time=-15:30:395;
            end

            ESM = daily_interp_dple_member(range,Time,Tdays,...
                TEMP_150m,TEMP_bottom,POC_FLUX_IN_bottom,LzooC_150m,Lzoo_loss_150m);


            %% 5. Run FEISTY DPLE -----------------------------------
            if y==1
                %! Initialize FEISTY
                % FEISTY INIT IS DEC NOT OCT, FIX!!!
                % this will have to be the biomass from the last mo before the start year
                % from the FOSI
                init_sim = ['v14_All_fish03_' simname];
                load([Cdir,simname,'/Last_mo_FOSI_' num2str(StartYr(yr)) '_' init_sim '.mat']);
                BENT.mass = BENT.bio;
                [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT] = sub_init_fish(ID,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

                frange = 305:param.DT:DAYS;
            else
                frange = 1:param.DT:DAYS;
            end

            for DAY = frange % days

                %%%! Future time step
                DY = int64(ceil(DAY));
                [Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT,ENVR] = ...
                    sub_futbio_1meso(DY,ESM,GRD,Sml_f,Sml_p,Sml_d,...
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

                P_Sml_f(:,DY) = Sml_f.prod;
                P_Sml_p(:,DY) = Sml_p.prod;
                P_Sml_d(:,DY) = Sml_d.prod;
                P_Med_f(:,DY) = Med_f.prod;
                P_Med_p(:,DY) = Med_p.prod;
                P_Med_d(:,DY) = Med_d.prod;
                P_Lrg_p(:,DY) = Lrg_p.prod;
                P_Lrg_d(:,DY) = Lrg_d.prod;

                C_Med_f(:,DY) = Med_f.caught;
                C_Med_p(:,DY) = Med_p.caught;
                C_Med_d(:,DY) = Med_d.caught;
                C_Lrg_p(:,DY) = Lrg_p.caught;
                C_Lrg_d(:,DY) = Lrg_d.caught;

            end %Days

            %! Calculate monthly means and save
            aa = (cumsum(MNTH)+1);
            a = [1,aa(1:end-1)]; % start of the month
            b = cumsum(MNTH);    % end of the month
            for i = 1:12
                MNT = MNT+1;     % Update monthly ticker

                %! Put vars of netcdf file
                netcdf.putVar(ncidB,vidbioB,[0 MNT-1],[NX 1],mean(S_Bent_bio(:,a(i):b(i)),2));
                netcdf.putVar(ncidB,vidTB,MNT-1,1,MNT);

                netcdf.putVar(ncidMZ,vidfracMZ,[0 MNT-1],[NX 1],mean(S_Mzoo_frac(:,a(i):b(i)),2));
                netcdf.putVar(ncidMZ,vidTMZ,MNT-1,1,MNT);

                netcdf.putVar(ncidSF,vidbioSF,[0 MNT-1],[NX 1],mean(S_Sml_f(:,a(i):b(i)),2));
                netcdf.putVar(ncidSP,vidbioSP,[0 MNT-1],[NX 1],mean(S_Sml_p(:,a(i):b(i)),2));
                netcdf.putVar(ncidSD,vidbioSD,[0 MNT-1],[NX 1],mean(S_Sml_d(:,a(i):b(i)),2));
                netcdf.putVar(ncidMF,vidbioMF,[0 MNT-1],[NX 1],mean(S_Med_f(:,a(i):b(i)),2));
                netcdf.putVar(ncidMP,vidbioMP,[0 MNT-1],[NX 1],mean(S_Med_p(:,a(i):b(i)),2));
                netcdf.putVar(ncidMD,vidbioMD,[0 MNT-1],[NX 1],mean(S_Med_d(:,a(i):b(i)),2));
                netcdf.putVar(ncidLP,vidbioLP,[0 MNT-1],[NX 1],mean(S_Lrg_p(:,a(i):b(i)),2));
                netcdf.putVar(ncidLD,vidbioLD,[0 MNT-1],[NX 1],mean(S_Lrg_d(:,a(i):b(i)),2));

                netcdf.putVar(ncidSF,vidprodSF,[0 MNT-1],[NX 1],mean(P_Sml_f(:,a(i):b(i)),2));
                netcdf.putVar(ncidSP,vidprodSP,[0 MNT-1],[NX 1],mean(P_Sml_p(:,a(i):b(i)),2));
                netcdf.putVar(ncidSD,vidprodSD,[0 MNT-1],[NX 1],mean(P_Sml_d(:,a(i):b(i)),2));
                netcdf.putVar(ncidMF,vidprodMF,[0 MNT-1],[NX 1],mean(P_Med_f(:,a(i):b(i)),2));
                netcdf.putVar(ncidMP,vidprodMP,[0 MNT-1],[NX 1],mean(P_Med_p(:,a(i):b(i)),2));
                netcdf.putVar(ncidMD,vidprodMD,[0 MNT-1],[NX 1],mean(P_Med_d(:,a(i):b(i)),2));
                netcdf.putVar(ncidLP,vidprodLP,[0 MNT-1],[NX 1],mean(P_Lrg_p(:,a(i):b(i)),2));
                netcdf.putVar(ncidLD,vidprodLD,[0 MNT-1],[NX 1],mean(P_Lrg_d(:,a(i):b(i)),2));

                netcdf.putVar(ncidMF,vidcatchMF,[0 MNT-1],[NX 1],mean(C_Med_f(:,a(i):b(i)),2));
                netcdf.putVar(ncidMP,vidcatchMP,[0 MNT-1],[NX 1],mean(C_Med_p(:,a(i):b(i)),2));
                netcdf.putVar(ncidMD,vidcatchMD,[0 MNT-1],[NX 1],mean(C_Med_d(:,a(i):b(i)),2));
                netcdf.putVar(ncidLP,vidcatchLP,[0 MNT-1],[NX 1],mean(C_Lrg_p(:,a(i):b(i)),2));
                netcdf.putVar(ncidLD,vidcatchLD,[0 MNT-1],[NX 1],mean(C_Lrg_d(:,a(i):b(i)),2));

            end %Monthly mean

        end %Lead year loop

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
        netcdf.close(ncidB);
        netcdf.close(ncidMZ);

    end %loop over M members

end %loop over Y years
