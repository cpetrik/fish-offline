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
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
%spath ='/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';

% Coef from MLR with interactions, every driver 0-2 lags
load([spath 'LME_FishMIP_cme_alltypes_mlr_coeffs_ALLdiv2SD_alllags_v2_noint.mat'])

%% forcing ---------------------------------------------------------
cnam = {'Det','TB','TP','Zmeso'};
lid = Acoef(:,21);

did = 2:4;
bid = 5:7;
pid = 8:10;
zid = 11:13;

LFtab = nan*ones(length(lid),length(cnam));
LPtab = nan*ones(length(lid),length(cnam));
LDtab = nan*ones(length(lid),length(cnam));
LAtab = nan*ones(length(lid),length(cnam));

%%
for L = 1:length(lid)

    %LME
    i = L;
    ilme = num2str(i);

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
    clear maxC fid


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
    

end %LME

%% Add R2 info back
LFtab(:,5) = Fcoef(:,14);
LPtab(:,5) = Pcoef(:,14);
LDtab(:,5) = Dcoef(:,14);
LAtab(:,5) = Acoef(:,14);

LFtab(:,6) = Fcoef(:,21);
LPtab(:,6) = Pcoef(:,21);
LDtab(:,6) = Dcoef(:,21);
LAtab(:,6) = Acoef(:,21);

%%
cnam = {'Det','TB','TP','ZmLoss','R2','LME'};

lname = cellstr(num2str(lid));

% cnam, lname, drivers2
Ftab1 = array2table(LFtab);
Ftab1.Properties.VariableNames = cnam;

Ptab1 = array2table(LPtab);
Ptab1.Properties.VariableNames = cnam;

Dtab1 = array2table(LDtab);
Dtab1.Properties.VariableNames = cnam;

Atab1 = array2table(LAtab);
Atab1.Properties.VariableNames = cnam;

%%
writetable(Ftab1,[spath,'LMEs_mlr_cme_drivers_ALLdiv2SD_alllags3_noint_reduc_F.csv'],...
    'Delimiter',',','WriteRowNames',false);
writetable(Ptab1,[spath,'LMEs_mlr_cme_drivers_ALLdiv2SD_alllags3_noint_reduc_P.csv'],...
    'Delimiter',',','WriteRowNames',false);
writetable(Dtab1,[spath,'LMEs_mlr_cme_drivers_ALLdiv2SD_alllags3_noint_reduc_D.csv'],...
    'Delimiter',',','WriteRowNames',false);
writetable(Atab1,[spath,'LMEs_mlr_cme_drivers_ALLdiv2SD_alllags3_noint_reduc_A.csv'],...
    'Delimiter',',','WriteRowNames',false);


save([spath,'LMEs_mlr_cme_drivers_ALLdiv2SD_alllags3_noint_reduc.mat'],...
    'LFtab','LPtab','LDtab','LAtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','lid');

