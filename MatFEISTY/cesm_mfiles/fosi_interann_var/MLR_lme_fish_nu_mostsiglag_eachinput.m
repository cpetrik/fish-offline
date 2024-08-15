% Forcing-fish mult linear regressions
% divided both ts by 2SD
% Find most signif time lag for driver
% Keep that coef to reduce complexity of heatmap

clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
%spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
spath ='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';

% Coef from MLR with interactions, every driver 0-2 lags
load([spath 'LME_alltypes_nu_mlr_coeffs_ALLdiv2SD_alllags.mat'])

%% forcing ---------------------------------------------------------
cnam = {'Det','TB','TP','Zmeso','ZmLoss','DetTB','ZmesoTP','ZmLossTP'};
lid = LME;

did = 1:3;
bid = 4:6;
pid = 7:9;
zid = 10:12;
mid = 13:15;
dtb = 16:18;
ztp = [19,21,23];
mtp = [20,22,24];

LFtab = nan*ones(length(lid),length(cnam));
LPtab = nan*ones(length(lid),length(cnam));
LDtab = nan*ones(length(lid),length(cnam));
LAtab = nan*ones(length(lid),length(cnam));

%%
for L = 1:length(lid)

    %LME
    i = L;
    ilme = LME(i);

    %% Forage
    %Det
    maxC = max(abs(Fcoef(i,did)));
    if(~isnan(maxC))
        fid = find(abs(Fcoef(i,:))==maxC);
        LFtab(L,1) = Fcoef(i,fid);
    end
    clear maxC
    %TB
    maxC = max(abs(Fcoef(i,bid)));
    if(~isnan(maxC))
        fid = find(abs(Fcoef(i,:))==maxC);
        LFtab(L,2) = Fcoef(i,fid);
    end
    clear maxC
    %TP
    maxC = max(abs(Fcoef(i,pid)));
    if(~isnan(maxC))
        fid = find(abs(Fcoef(i,:))==maxC);
        LFtab(L,3) = Fcoef(i,fid);
    end
    clear maxC
    %Zmeso
    maxC = max(abs(Fcoef(i,zid)));
    if(~isnan(maxC))
        fid = find(abs(Fcoef(i,:))==maxC);
        LFtab(L,4) = Fcoef(i,fid);
    end
    clear maxC
    %ZmLoss
    maxC = max(abs(Fcoef(i,mid)));
    if(~isnan(maxC))
        fid = find(abs(Fcoef(i,:))==maxC);
        LFtab(L,5) = Fcoef(i,fid);
    end
    clear maxC
    %Det:TB
    maxC = max(abs(Fcoef(i,dtb)));
    if(~isnan(maxC))
        fid = find(abs(Fcoef(i,:))==maxC);
        LFtab(L,6) = Fcoef(i,fid);
    end
    clear maxC
    %Zmeso:TP
    maxC = max(abs(Fcoef(i,ztp)));
    if(~isnan(maxC))
        fid = find(abs(Fcoef(i,:))==maxC);
        LFtab(L,7) = Fcoef(i,fid);
    end
    clear maxC
    %ZmLoss:TP
    maxC = max(abs(Fcoef(i,mtp)));
    if(~isnan(maxC))
        fid = find(abs(Fcoef(i,:))==maxC);
        LFtab(L,8) = Fcoef(i,fid);
    end
    clear maxC

    %% Pelagic
    %Det
    maxC = max(abs(Pcoef(i,did)));
    if(~isnan(maxC))
        fid = find(abs(Pcoef(i,:))==maxC);
        LPtab(L,1) = Pcoef(i,fid);
    end
    clear maxC
    %TB
    maxC = max(abs(Pcoef(i,bid)));
    if(~isnan(maxC))
        fid = find(abs(Pcoef(i,:))==maxC);
        LPtab(L,2) = Pcoef(i,fid);
    end
    clear maxC
    %TP
    maxC = max(abs(Pcoef(i,pid)));
    if(~isnan(maxC))
        fid = find(abs(Pcoef(i,:))==maxC);
        LPtab(L,3) = Pcoef(i,fid);
    end
    clear maxC
    %Zmeso
    maxC = max(abs(Pcoef(i,zid)));
    if(~isnan(maxC))
        fid = find(abs(Pcoef(i,:))==maxC);
        LPtab(L,4) = Pcoef(i,fid);
    end
    clear maxC
    %ZmLoss
    maxC = max(abs(Pcoef(i,mid)));
    if(~isnan(maxC))
        fid = find(abs(Pcoef(i,:))==maxC);
        LPtab(L,5) = Pcoef(i,fid);
    end
    clear maxC
    %Det:TB
    maxC = max(abs(Pcoef(i,dtb)));
    if(~isnan(maxC))
        fid = find(abs(Pcoef(i,:))==maxC);
        LPtab(L,6) = Pcoef(i,fid);
    end
    clear maxC
    %Zmeso:TP
    maxC = max(abs(Pcoef(i,ztp)));
    if(~isnan(maxC))
        fid = find(abs(Pcoef(i,:))==maxC);
        LPtab(L,7) = Pcoef(i,fid);
    end
    clear maxC
    %ZmLoss:TP
    maxC = max(abs(Pcoef(i,mtp)));
    if(~isnan(maxC))
        fid = find(abs(Pcoef(i,:))==maxC);
        LPtab(L,8) = Pcoef(i,fid);
    end
    clear maxC


    %% Demersal
    %Det
    maxC = max(abs(Dcoef(i,did)));
    if(~isnan(maxC))
        fid = find(abs(Dcoef(i,:))==maxC);
        LDtab(L,1) = Dcoef(i,fid);
    end
    clear maxC
    %TB
    maxC = max(abs(Dcoef(i,bid)));
    if(~isnan(maxC))
        fid = find(abs(Dcoef(i,:))==maxC);
        LDtab(L,2) = Dcoef(i,fid);
    end
    clear maxC
    %TP
    maxC = max(abs(Dcoef(i,pid)));
    if(~isnan(maxC))
        fid = find(abs(Dcoef(i,:))==maxC);
        LDtab(L,3) = Dcoef(i,fid);
    end
    clear maxC
    %Zmeso
    maxC = max(abs(Dcoef(i,zid)));
    if(~isnan(maxC))
        fid = find(abs(Dcoef(i,:))==maxC);
        LDtab(L,4) = Dcoef(i,fid);
    end
    clear maxC
    %ZmLoss
    maxC = max(abs(Dcoef(i,mid)));
    if(~isnan(maxC))
        fid = find(abs(Dcoef(i,:))==maxC);
        LDtab(L,5) = Dcoef(i,fid);
    end
    clear maxC
    %Det:TB
    maxC = max(abs(Dcoef(i,dtb)));
    if(~isnan(maxC))
        fid = find(abs(Dcoef(i,:))==maxC);
        LDtab(L,6) = Dcoef(i,fid);
    end
    clear maxC
    %Zmeso:TP
    maxC = max(abs(Dcoef(i,ztp)));
    if(~isnan(maxC))
        fid = find(abs(Dcoef(i,:))==maxC);
        LDtab(L,7) = Dcoef(i,fid);
    end
    clear maxC
    %ZmLoss:TP
    maxC = max(abs(Dcoef(i,mtp)));
    if(~isnan(maxC))
        fid = find(abs(Dcoef(i,:))==maxC);
        LDtab(L,8) = Dcoef(i,fid);
    end
    clear maxC

    %% All
    %Det
    maxC = max(abs(Acoef(i,did)));
    if(~isnan(maxC))
        fid = find(abs(Acoef(i,:))==maxC);
        LAtab(L,1) = Acoef(i,fid);
    end
    clear maxC
    %TB
    maxC = max(abs(Acoef(i,bid)));
    if(~isnan(maxC))
        fid = find(abs(Acoef(i,:))==maxC);
        LAtab(L,2) = Acoef(i,fid);
    end
    clear maxC
    %TP
    maxC = max(abs(Acoef(i,pid)));
    if(~isnan(maxC))
        fid = find(abs(Acoef(i,:))==maxC);
        LAtab(L,3) = Acoef(i,fid);
    end
    clear maxC
    %Zmeso
    maxC = max(abs(Acoef(i,zid)));
    if(~isnan(maxC))
        fid = find(abs(Acoef(i,:))==maxC);
        LAtab(L,4) = Acoef(i,fid);
    end
    clear maxC
    %ZmLoss
    maxC = max(abs(Acoef(i,mid)));
    if(~isnan(maxC))
        fid = find(abs(Acoef(i,:))==maxC);
        LAtab(L,5) = Acoef(i,fid);
    end
    clear maxC
    %Det:TB
    maxC = max(abs(Acoef(i,dtb)));
    if(~isnan(maxC))
        fid = find(abs(Acoef(i,:))==maxC);
        LAtab(L,6) = Acoef(i,fid);
    end
    clear maxC
    %Zmeso:TP
    maxC = max(abs(Acoef(i,ztp)));
    if(~isnan(maxC))
        fid = find(abs(Acoef(i,:))==maxC);
        LAtab(L,7) = Acoef(i,fid);
    end
    clear maxC
    %ZmLoss:TP
    maxC = max(abs(Acoef(i,mtp)));
    if(~isnan(maxC))
        fid = find(abs(Acoef(i,:))==maxC);
        LAtab(L,8) = Acoef(i,fid);
    end
    clear maxC


end %LME

%%
lname = cellstr(num2str(lid));
% cnam, lname, drivers2
Ftab1 = array2table(LFtab,"RowNames",lname);
Ftab1.Properties.VariableNames = cnam;

Ptab1 = array2table(LPtab,"RowNames",lname);
Ptab1.Properties.VariableNames = cnam;

Dtab1 = array2table(LDtab,"RowNames",lname);
Dtab1.Properties.VariableNames = cnam;

Atab1 = array2table(LAtab,"RowNames",lname);
Atab1.Properties.VariableNames = cnam;

%%
writetable(Ftab1,[spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_alllags_reduc_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_alllags_reduc_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_alllags_reduc_D.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Atab1,[spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_alllags_reduc_A.csv'],...
    'Delimiter',',','WriteRowNames',true);


save([spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_alllags_reduc.mat'],...
    'LFtab','LPtab','LDtab','LAtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','lid');

