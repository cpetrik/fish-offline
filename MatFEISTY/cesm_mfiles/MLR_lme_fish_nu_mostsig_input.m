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
load([spath 'LME_nu_mlr_coeff_pvals_ALLdiv2SD.mat'])

%%
Fcoef = LMEFnumlrcoeffsALLdiv2SD;
Fpval = LMEFnumlrpvalsALLdiv2SD;
Pcoef = LMEPnumlrcoeffsALLdiv2SD;
Ppval = LMEPnumlrpvalsALLdiv2SD;
Dcoef = LMEDnumlrcoeffsALLdiv2SD;
Dpval = LMEDnumlrpvalsALLdiv2SD;
Acoef = LMEAnumlrcoeffsALLdiv2SD;
Apval = LMEAnumlrpvalsALLdiv2SD;

%% forcing ---------------------------------------------------------
cnam = {'coef','p','idriver','driver'};

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
lid = LME;

LFtab = nan*ones(length(lid),3);
LPtab = nan*ones(length(lid),3);
LDtab = nan*ones(length(lid),3);
LAtab = nan*ones(length(lid),3);

LFt = cell(length(lid),1);
LPt = cell(length(lid),1);
LDt = cell(length(lid),1);
LAt = cell(length(lid),1);

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
        LFt(L) = cDrivers(fid);
    end
    clear maxC

    maxC = max(abs(Pcoef(i,:)));
    if(~isnan(maxC))
        pid = find(abs(Pcoef(i,:))==maxC);
        LPtab(L,1) = Pcoef(i,pid);
        LPtab(L,2) = Ppval(i,pid);
        LPtab(L,3) = (pid);
        LPt(L) = cDrivers(pid);
    end
    clear maxC

    maxC = max(abs(Dcoef(i,:)));
    if(~isnan(maxC))
        did = find(abs(Dcoef(i,:))==maxC);
        LDtab(L,1) = Dcoef(i,did);
        LDtab(L,2) = Dpval(i,did);
        LDtab(L,3) = (did);
        LDt(L) = cDrivers(did);
    end
    clear maxC

    maxC = max(abs(Acoef(i,:)));
    if(~isnan(maxC))
        aid = find(abs(Acoef(i,:))==maxC);
        LAtab(L,1) = Acoef(i,aid);
        LAtab(L,2) = Apval(i,aid);
        LAtab(L,3) = (aid);
        LAt(L) = cDrivers(aid);
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

%%
writetable(Ftab1,[spath,'LMEs_mlr_drivers_ALLdiv2SD_maxcoef_Fnu.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Ptab1,[spath,'LMEs_mlr_drivers_ALLdiv2SD_maxcoef_Pnu.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Dtab1,[spath,'LMEs_mlr_drivers_ALLdiv2SD_maxcoef_Dnu.csv'],...
    'Delimiter',',','WriteRowNames',true);
writetable(Atab1,[spath,'LMEs_mlr_drivers_ALLdiv2SD_maxcoef_Anu.csv'],...
    'Delimiter',',','WriteRowNames',true);

save([spath,'LMEs_mlr_nu_drivers_ALLdiv2SD_maxcoefs.mat'],...
    'LFtab','LPtab','LDtab','LAtab',...
    'Ftab1','Ptab1','Dtab1','Atab1',...
    'lid','LME','cDrivers','sDrivers','cnam');

