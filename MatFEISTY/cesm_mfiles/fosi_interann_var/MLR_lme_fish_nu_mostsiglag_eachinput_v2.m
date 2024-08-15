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
load([spath 'LME_nu_mlr_coeffs_ALLdiv2SD_alllags3_R2_07_cluster.mat'])

%% forcing ---------------------------------------------------------
lid = Anu(:,1);

did = 3:5;
bid = 6:8;
pid = 9:11;
zid = 12:14;

cnam = {'Det','TB','TP','ZmLoss'};

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
    maxC = max(abs(Fnu(i,did)));
    if(~isnan(maxC))
        fid = find(abs(Fnu(i,:))==maxC);
        LFtab(L,1) = Fnu(i,fid);
    end
    clear maxC
    %TB
    maxC = max(abs(Fnu(i,bid)));
    if(~isnan(maxC))
        fid = find(abs(Fnu(i,:))==maxC);
        LFtab(L,2) = Fnu(i,fid);
    end
    clear maxC
    %TP
    maxC = max(abs(Fnu(i,pid)));
    if(~isnan(maxC))
        fid = find(abs(Fnu(i,:))==maxC);
        LFtab(L,3) = Fnu(i,fid);
    end
    clear maxC
    %Zmeso
    maxC = max(abs(Fnu(i,zid)));
    if(~isnan(maxC))
        fid = find(abs(Fnu(i,:))==maxC);
        LFtab(L,4) = Fnu(i,fid);
    end
    clear maxC

    %% Pelagic
    %Det
    maxC = max(abs(Pnu(i,did)));
    if(~isnan(maxC))
        fid = find(abs(Pnu(i,:))==maxC);
        LPtab(L,1) = Pnu(i,fid);
    end
    clear maxC
    %TB
    maxC = max(abs(Pnu(i,bid)));
    if(~isnan(maxC))
        fid = find(abs(Pnu(i,:))==maxC);
        LPtab(L,2) = Pnu(i,fid);
    end
    clear maxC
    %TP
    maxC = max(abs(Pnu(i,pid)));
    if(~isnan(maxC))
        fid = find(abs(Pnu(i,:))==maxC);
        LPtab(L,3) = Pnu(i,fid);
    end
    clear maxC
    %Zmeso
    maxC = max(abs(Pnu(i,zid)));
    if(~isnan(maxC))
        fid = find(abs(Pnu(i,:))==maxC);
        LPtab(L,4) = Pnu(i,fid);
    end
    clear maxC
    

    %% Demersal
    %Det
    maxC = max(abs(Dnu(i,did)));
    if(~isnan(maxC))
        fid = find(abs(Dnu(i,:))==maxC);
        LDtab(L,1) = Dnu(i,fid);
    end
    clear maxC
    %TB
    maxC = max(abs(Dnu(i,bid)));
    if(~isnan(maxC))
        fid = find(abs(Dnu(i,:))==maxC);
        LDtab(L,2) = Dnu(i,fid);
    end
    clear maxC
    %TP
    maxC = max(abs(Dnu(i,pid)));
    if(~isnan(maxC))
        fid = find(abs(Dnu(i,:))==maxC);
        LDtab(L,3) = Dnu(i,fid);
    end
    clear maxC
    %Zmeso
    maxC = max(abs(Dnu(i,zid)));
    if(~isnan(maxC))
        fid = find(abs(Dnu(i,:))==maxC);
        LDtab(L,4) = Dnu(i,fid);
    end
    clear maxC
    

    %% All
    %Det
    maxC = max(abs(Anu(i,did)));
    if(~isnan(maxC))
        fid = find(abs(Anu(i,:))==maxC);
        LAtab(L,1) = Anu(i,fid);
    end
    clear maxC
    %TB
    maxC = max(abs(Anu(i,bid)));
    if(~isnan(maxC))
        fid = find(abs(Anu(i,:))==maxC);
        LAtab(L,2) = Anu(i,fid);
    end
    clear maxC
    %TP
    maxC = max(abs(Anu(i,pid)));
    if(~isnan(maxC))
        fid = find(abs(Anu(i,:))==maxC);
        LAtab(L,3) = Anu(i,fid);
    end
    clear maxC
    %Zmeso
    maxC = max(abs(Anu(i,zid)));
    if(~isnan(maxC))
        fid = find(abs(Anu(i,:))==maxC);
        LAtab(L,4) = Anu(i,fid);
    end
    clear maxC
    

end %LME

%% Add R2 info back
LFtab(:,5) = Fnu(:,15);
LPtab(:,5) = Pnu(:,15);
LDtab(:,5) = Dnu(:,15);
LAtab(:,5) = Anu(:,15);

LFtab(:,6) = Fnu(:,1);
LPtab(:,6) = Pnu(:,1);
LDtab(:,6) = Dnu(:,1);
LAtab(:,6) = Anu(:,1);

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
writetable(Ftab1,[spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_alllags3_noint_reduc_F.csv'],...
    'Delimiter',',','WriteRowNames',false);
writetable(Ptab1,[spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_alllags3_noint_reduc_P.csv'],...
    'Delimiter',',','WriteRowNames',false);
writetable(Dtab1,[spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_alllags3_noint_reduc_D.csv'],...
    'Delimiter',',','WriteRowNames',false);
writetable(Atab1,[spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_alllags3_noint_reduc_A.csv'],...
    'Delimiter',',','WriteRowNames',false);


save([spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_alllags3_noint_reduc.mat'],...
    'LFtab','LPtab','LDtab','LAtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','lid');

