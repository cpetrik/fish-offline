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
load([spath 'LME_biom_mlr_coeffs_ALLdiv2SD_alllags3_R2_05_cluster.mat'])

%% forcing ---------------------------------------------------------
lid = Abiom(:,1);

did = 3:5;
bid = 6:8;
pid = 9:11;
zid = 12:14;

cnam = {'Det','TB','TP','ZmLoss'};

LFtab = nan*ones(length(lid),length(cnam));
LPtab = nan*ones(length(lid),length(cnam));
LDtab = nan*ones(length(lid),length(cnam));
LAtab = nan*ones(length(lid),length(cnam));
LBtab = nan*ones(length(lid),length(cnam));

%%
for L = 1:length(lid)

    %LME
    i = L;
    ilme = num2str(i);

    %% Forage
    %Det
    maxC = max(abs(Fbiom(i,did)));
    if(~isnan(maxC))
        fid = find(abs(Fbiom(i,:))==maxC);
        LFtab(L,1) = Fbiom(i,fid);
    end
    clear maxC
    %TB
    maxC = max(abs(Fbiom(i,bid)));
    if(~isnan(maxC))
        fid = find(abs(Fbiom(i,:))==maxC);
        LFtab(L,2) = Fbiom(i,fid);
    end
    clear maxC
    %TP
    maxC = max(abs(Fbiom(i,pid)));
    if(~isnan(maxC))
        fid = find(abs(Fbiom(i,:))==maxC);
        LFtab(L,3) = Fbiom(i,fid);
    end
    clear maxC
    %Zmeso
    maxC = max(abs(Fbiom(i,zid)));
    if(~isnan(maxC))
        fid = find(abs(Fbiom(i,:))==maxC);
        LFtab(L,4) = Fbiom(i,fid);
    end
    clear maxC
    

    %% Pelagic
    %Det
    maxC = max(abs(Pbiom(i,did)));
    if(~isnan(maxC))
        fid = find(abs(Pbiom(i,:))==maxC);
        LPtab(L,1) = Pbiom(i,fid);
    end
    clear maxC
    %TB
    maxC = max(abs(Pbiom(i,bid)));
    if(~isnan(maxC))
        fid = find(abs(Pbiom(i,:))==maxC);
        LPtab(L,2) = Pbiom(i,fid);
    end
    clear maxC
    %TP
    maxC = max(abs(Pbiom(i,pid)));
    if(~isnan(maxC))
        fid = find(abs(Pbiom(i,:))==maxC);
        LPtab(L,3) = Pbiom(i,fid);
    end
    clear maxC
    %Zmeso
    maxC = max(abs(Pbiom(i,zid)));
    if(~isnan(maxC))
        fid = find(abs(Pbiom(i,:))==maxC);
        LPtab(L,4) = Pbiom(i,fid);
    end
    clear maxC
    


    %% Demersal
    %Det
    maxC = max(abs(Dbiom(i,did)));
    if(~isnan(maxC))
        fid = find(abs(Dbiom(i,:))==maxC);
        LDtab(L,1) = Dbiom(i,fid);
    end
    clear maxC
    %TB
    maxC = max(abs(Dbiom(i,bid)));
    if(~isnan(maxC))
        fid = find(abs(Dbiom(i,:))==maxC);
        LDtab(L,2) = Dbiom(i,fid);
    end
    clear maxC
    %TP
    maxC = max(abs(Dbiom(i,pid)));
    if(~isnan(maxC))
        fid = find(abs(Dbiom(i,:))==maxC);
        LDtab(L,3) = Dbiom(i,fid);
    end
    clear maxC
    %Zmeso
    maxC = max(abs(Dbiom(i,zid)));
    if(~isnan(maxC))
        fid = find(abs(Dbiom(i,:))==maxC);
        LDtab(L,4) = Dbiom(i,fid);
    end
    clear maxC
    

    %% All
    %Det
    maxC = max(abs(Abiom(i,did)));
    if(~isnan(maxC))
        fid = find(abs(Abiom(i,:))==maxC);
        LAtab(L,1) = Abiom(i,fid);
    end
    clear maxC
    %TB
    maxC = max(abs(Abiom(i,bid)));
    if(~isnan(maxC))
        fid = find(abs(Abiom(i,:))==maxC);
        LAtab(L,2) = Abiom(i,fid);
    end
    clear maxC
    %TP
    maxC = max(abs(Abiom(i,pid)));
    if(~isnan(maxC))
        fid = find(abs(Abiom(i,:))==maxC);
        LAtab(L,3) = Abiom(i,fid);
    end
    clear maxC
    %Zmeso
    maxC = max(abs(Abiom(i,zid)));
    if(~isnan(maxC))
        fid = find(abs(Abiom(i,:))==maxC);
        LAtab(L,4) = Abiom(i,fid);
    end
    clear maxC
    

    %% Benthos
    %Det
    maxC = max(abs(Bbiom(i,did)));
    if(~isnan(maxC))
        fid = find(abs(Bbiom(i,:))==maxC);
        LBtab(L,1) = Bbiom(i,fid);
    end
    clear maxC
    %TB
    maxC = max(abs(Bbiom(i,bid)));
    if(~isnan(maxC))
        fid = find(abs(Bbiom(i,:))==maxC);
        LBtab(L,2) = Bbiom(i,fid);
    end
    clear maxC
    %TP
    maxC = max(abs(Bbiom(i,pid)));
    if(~isnan(maxC))
        fid = find(abs(Bbiom(i,:))==maxC);
        LBtab(L,3) = Bbiom(i,fid);
    end
    clear maxC
    %Zmeso
    maxC = max(abs(Bbiom(i,zid)));
    if(~isnan(maxC))
        fid = find(abs(Bbiom(i,:))==maxC);
        LBtab(L,4) = Bbiom(i,fid);
    end
    clear maxC
    

end %LME

%% Add R2 and LME back
cnam = {'Det','TB','TP','ZmLoss','R2','LME'};

LFtab(:,5) = Fbiom(:,15);
LPtab(:,5) = Pbiom(:,15);
LDtab(:,5) = Dbiom(:,15);
LAtab(:,5) = Abiom(:,15);
LBtab(:,5) = Bbiom(:,15);

LFtab(:,6) = Fbiom(:,1);
LPtab(:,6) = Pbiom(:,1);
LDtab(:,6) = Dbiom(:,1);
LAtab(:,6) = Abiom(:,1);
LBtab(:,6) = Bbiom(:,1);

%%
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

Btab1 = array2table(LBtab);
Btab1.Properties.VariableNames = cnam;

%%
writetable(Ftab1,[spath,'LMEs_mlr_biom_drivers_ALLdiv2SD_alllags3_noint_reduc_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[spath,'LMEs_mlr_biom_drivers_ALLdiv2SD_alllags3_noint_reduc_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[spath,'LMEs_mlr_biom_drivers_ALLdiv2SD_alllags3_noint_reduc_D.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Atab1,[spath,'LMEs_mlr_biom_drivers_ALLdiv2SD_alllags3_noint_reduc_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Btab1,[spath,'LMEs_mlr_biom_drivers_ALLdiv2SD_alllags3_noint_reduc_B.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'LMEs_mlr_biom_drivers_ALLdiv2SD_alllags3_noint_reduc.mat'],...
    'LFtab','LPtab','LDtab','LAtab','LBtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','Btab1','lid');

