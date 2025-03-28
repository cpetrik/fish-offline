% Find most signif driver in forcing-fish mult linear regressions
% divide both ts by 2SD


clear
close all

ppath = "/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/corrs/";

%% FOSI input forcing
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% Fish data
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];

sims = {'v15_All_fish03';'v15_climatol';'v15_varFood';'v15_varTemp'};
mod = sims{1};

% Anoms with linear trend removed
load([spath 'LME_biom_mlr_coeff_pvals_ALLdiv2SD_lag3.mat'])

%%
Fcoef = LMEFmlrcoeffsALLdiv2SDlag3;
Fpval = LMEFmlrpvalsALLdiv2SDlag3;
Pcoef = LMEPmlrcoeffsALLdiv2SDlag3;
Ppval = LMEPmlrpvalsALLdiv2SDlag3;
Dcoef = LMEDmlrcoeffsALLdiv2SDlag3;
Dpval = LMEDmlrpvalsALLdiv2SDlag3;
Acoef = LMEAmlrcoeffsALLdiv2SDlag3;
Apval = LMEAmlrpvalsALLdiv2SDlag3;
Bcoef = LMEBmlrcoeffsALLdiv2SDlag3;
Bpval = LMEBmlrpvalsALLdiv2SDlag3;

%% forcing ---------------------------------------------------------
cnam = {'coef','p','idriver','driver'};
lid = LME;

LFtab = nan*ones(length(lid),3);
LPtab = nan*ones(length(lid),3);
LDtab = nan*ones(length(lid),3);
LAtab = nan*ones(length(lid),3);
LBtab = nan*ones(length(lid),3);

LFt = cell(length(lid),1);
LPt = cell(length(lid),1);
LDt = cell(length(lid),1);
LAt = cell(length(lid),1);
LBt = cell(length(lid),1);

%%
for L = 1:length(lid)

    %LME
    i = L;
    ilme = LME(i);

    %%
    maxC = max(abs(Fcoef(i,:)));
    if(~isnan(maxC))
        fid = find(abs(Fcoef(i,:))==maxC);
        LFtab(L,1) = Fcoef(i,fid);
        LFtab(L,2) = Fpval(i,fid);
        LFtab(L,3) = (fid);
        LFt(L) = cDriver(fid);
    end
    clear maxC

    maxC = max(abs(Pcoef(i,:)));
    if(~isnan(maxC))
        pid = find(abs(Pcoef(i,:))==maxC);
        LPtab(L,1) = Pcoef(i,pid);
        LPtab(L,2) = Ppval(i,pid);
        LPtab(L,3) = (pid);
        LPt(L) = cDriver(pid);
    end
    clear maxC

    maxC = max(abs(Dcoef(i,:)));
    if(~isnan(maxC))
        did = find(abs(Dcoef(i,:))==maxC);
        LDtab(L,1) = Dcoef(i,did);
        LDtab(L,2) = Dpval(i,did);
        LDtab(L,3) = (did);
        LDt(L) = cDriver(did);
    end
    clear maxC

    maxC = max(abs(Acoef(i,:)));
    if(~isnan(maxC))
        aid = find(abs(Acoef(i,:))==maxC);
        LAtab(L,1) = Acoef(i,aid);
        LAtab(L,2) = Apval(i,aid);
        LAtab(L,3) = (aid);
        LAt(L) = cDriver(aid);
    end
    clear maxC

    maxC = max(abs(Bcoef(i,:)));
    if(~isnan(maxC))
        cid = find(abs(Bcoef(i,:))==maxC);
        LBtab(L,1) = Bcoef(i,cid);
        LBtab(L,2) = Bpval(i,cid);
        LBtab(L,3) = (cid);
        LBt(L) = cDriver(cid);
    end
    clear maxC

end %LME

%%
lname = cellstr(num2str(lid));
% cnam, lname, drivers2
Ftab1 = array2table(LFtab,"RowNames",lname);
Ftab1(:,4) = LFt;
Ftab1.Properties.VariableNames = cnam;

Ptab1 = array2table(LPtab,"RowNames",lname);
Ptab1(:,4) = LPt;
Ptab1.Properties.VariableNames = cnam;

Dtab1 = array2table(LDtab,"RowNames",lname);
Dtab1(:,4) = LDt;
Dtab1.Properties.VariableNames = cnam;

Atab1 = array2table(LAtab,"RowNames",lname);
Atab1(:,4) = LAt;
Atab1.Properties.VariableNames = cnam;

Btab1 = array2table(LBtab,"RowNames",lname);
Btab1(:,4) = LBt;
Btab1.Properties.VariableNames = cnam;

%%
writetable(Ftab1,[spath,'LMEs_mlr_drivers_ALLdiv2SD_maxcoef_F.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[spath,'LMEs_mlr_drivers_ALLdiv2SD_maxcoef_P.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[spath,'LMEs_mlr_drivers_ALLdiv2SD_maxcoef_D.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Atab1,[spath,'LMEs_mlr_drivers_ALLdiv2SD_maxcoef_A.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Btab1,[spath,'LMEs_mlr_drivers_ALLdiv2SD_maxcoef_B.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'LMEs_mlr_drivers_ALLdiv2SD_maxcoefs.mat'],...
    'LFtab','LPtab','LDtab','LAtab','LBtab',...
    'Ftab1','Ptab1','Dtab1','Atab1','Btab1','lid');

