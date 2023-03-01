% Corr LME biomass of FEISTY experiments
% with full (and each other?)
% CESM FOSI

clear
close all

%% Map data
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';
load([cpath 'gridspec_POP_gx1v6_noSeas.mat']);
load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
load([cpath 'LME-mask-POP_gx1v6.mat']);

[ni,nj]=size(TLONG);
ID = GRD.ID;

plotminlat=-90; %Set these bounds for your data
plotmaxlat=90;
plotminlon=-280;
plotmaxlon=80;
latlim=[plotminlat plotmaxlat];
lonlim=[plotminlon plotmaxlon];

load coastlines;

%% Fish data
%cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_nmort1_BE08_noCC_RE00100';
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];

pp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
ppath = [pp cfile '/corrs/'];
if (~isfolder(ppath))
    mkdir(ppath)
end

sims = {'v15_All_fish03_';'v15_climatol_';'v15_varTemp_';'v15_varFood_'};

%% Full var
mod = sims{1};

load([fpath 'FEISTY_FOSI_',mod,'lme_ann_mean_anoms.mat'],'as','am','al','af','ap','ad','aa','ab'); % Anoms with linear trend removed
vas = as;
vam = am;
val = al;
vaf = af;
vap = ap;
vad = ad;
vaa = aa;
vab = ab;

clear as am al af ap ad aa ab

%% Clim
mod = sims{2};

load([fpath 'FEISTY_FOSI_',mod,'lme_ann_mean_anoms.mat'],'as','am','al','af','ap','ad','aa','ab'); % Anoms with linear trend removed
cas = as;
cam = am;
cal = al;
caf = af;
cap = ap;
cad = ad;
caa = aa;
cab = ab;

clear as am al af ap ad aa ab

%% varTemp
mod = sims{3};

load([fpath 'FEISTY_FOSI_',mod,'lme_ann_mean_anoms.mat'],'as','am','al','af','ap','ad','aa','ab'); % Anoms with linear trend removed
tas = as;
tam = am;
tal = al;
taf = af;
tap = ap;
tad = ad;
taa = aa;
tab = ab;

clear as am al af ap ad aa ab

%% varFood
mod = sims{4};

load([fpath 'FEISTY_FOSI_',mod,'lme_ann_mean_anoms.mat'],'as','am','al','af','ap','ad','aa','ab'); % Anoms with linear trend removed
fas = as;
fam = am;
fal = al;
faf = af;
fap = ap;
fad = ad;
faa = aa;
fab = ab;

clear as am al af ap ad aa ab

%% Loop over all LMEs and all Climate - Use anomalies of fish mean biomass
% US LMEs only
lid = [54,1:2,10,3,5:7,65];
lname = {'CHK','EBS','GAK','HI','CCE','GMX','SE','NE','AI'};

colororder({'k','b'})
close all

%% Corrs

rS = nan*ones(66,3);
rM = rS;
rL = rS;
rF = rS;
rP = rS;
rD = rS;
rA = rS;
rB = rS;
pS = rS;
pM = rS;
pL = rS;
pF = rS;
pP = rS;
pD = rS;
pA = rS;
pB = rS;

for i=1:66 %length(lid)
%     lme = lid(i);
%     ltex = lname{i};
    lme = i;

    %% Corr
    % Clim
    [rs,ps] = corrcoef(cas(lme,:),vas(lme,:));
    rS(i,1) = rs(1,2); pS(i,1) = ps(1,2);

    [rm,pm] = corrcoef(cam(lme,:),vam(lme,:));
    rM(i,1) = rm(1,2); pM(i,1) = pm(1,2);

    [rl,pl] = corrcoef(cal(lme,:),val(lme,:));
    rL(i,1) = rl(1,2); pL(i,1) = pl(1,2);

    [rf,pf] = corrcoef(caf(lme,:),vaf(lme,:));
    rF(i,1) = rf(1,2); pF(i,1) = pf(1,2);

    [rp,pp] = corrcoef(cap(lme,:),vap(lme,:));
    rP(i,1) = rp(1,2); pP(i,1) = pp(1,2);

    [rd,pd] = corrcoef(cad(lme,:),vad(lme,:));
    rD(i,1) = rd(1,2); pD(i,1) = pd(1,2);

    [ra,pa] = corrcoef(caa(lme,:),vaa(lme,:));
    rA(i,1) = ra(1,2); pA(i,1) = pa(1,2);

    [rb,pb] = corrcoef(cab(lme,:),vab(lme,:));
    rB(i,1) = rb(1,2); pB(i,1) = pb(1,2);

    % varTemp
    [rs,ps] = corrcoef(tas(lme,:),vas(lme,:));
    rS(i,2) = rs(1,2); pS(i,2) = ps(1,2);

    [rm,pm] = corrcoef(tam(lme,:),vam(lme,:));
    rM(i,2) = rm(1,2); pM(i,2) = pm(1,2);

    [rl,pl] = corrcoef(tal(lme,:),val(lme,:));
    rL(i,2) = rl(1,2); pL(i,2) = pl(1,2);

    [rf,pf] = corrcoef(taf(lme,:),vaf(lme,:));
    rF(i,2) = rf(1,2); pF(i,2) = pf(1,2);

    [rp,pp] = corrcoef(tap(lme,:),vap(lme,:));
    rP(i,2) = rp(1,2); pP(i,2) = pp(1,2);

    [rd,pd] = corrcoef(tad(lme,:),vad(lme,:));
    rD(i,2) = rd(1,2); pD(i,2) = pd(1,2);

    [ra,pa] = corrcoef(taa(lme,:),vaa(lme,:));
    rA(i,2) = ra(1,2); pA(i,2) = pa(1,2);

    [rb,pb] = corrcoef(tab(lme,:),vab(lme,:));
    rB(i,2) = rb(1,2); pB(i,2) = pb(1,2);

    % varFood
    [rs,ps] = corrcoef(fas(lme,:),vas(lme,:));
    rS(i,3) = rs(1,2); pS(i,3) = ps(1,2);

    [rm,pm] = corrcoef(fam(lme,:),vam(lme,:));
    rM(i,3) = rm(1,2); pM(i,3) = pm(1,2);

    [rl,pl] = corrcoef(fal(lme,:),val(lme,:));
    rL(i,3) = rl(1,2); pL(i,3) = pl(1,2);

    [rf,pf] = corrcoef(faf(lme,:),vaf(lme,:));
    rF(i,3) = rf(1,2); pF(i,3) = pf(1,2);

    [rp,pp] = corrcoef(fap(lme,:),vap(lme,:));
    rP(i,3) = rp(1,2); pP(i,3) = pp(1,2);

    [rd,pd] = corrcoef(fad(lme,:),vad(lme,:));
    rD(i,3) = rd(1,2); pD(i,3) = pd(1,2);

    [ra,pa] = corrcoef(faa(lme,:),vaa(lme,:));
    rA(i,3) = ra(1,2); pA(i,3) = pa(1,2);

    [rb,pb] = corrcoef(fab(lme,:),vab(lme,:));
    rB(i,3) = rb(1,2); pB(i,3) = pb(1,2);

end


%%
save([fpath 'LME_corrs_FOSI_fished_expers.mat'],...
    'rS','rM','rL','rF','rP','rD','rA','rB',...
    'pS','pM','pL','pF','pP','pD','pA','pB','sims');

%%
% writecell(sigS,[fpath 'LME_fosi_fished_',mod,'sigS_climate.csv']);
% writecell(sigM,[fpath 'LME_fosi_fished_',mod,'sigM_climate.csv']);
% writecell(sigL,[fpath 'LME_fosi_fished_',mod,'sigL_climate.csv']);
% writecell(sigF,[fpath 'LME_fosi_fished_',mod,'sigF_climate.csv']);
% writecell(sigP,[fpath 'LME_fosi_fished_',mod,'sigP_climate.csv']);
% writecell(sigD,[fpath 'LME_fosi_fished_',mod,'sigD_climate.csv']);
% writecell(sigA,[fpath 'LME_fosi_fished_',mod,'sigA_climate.csv']);
% writecell(sigB,[fpath 'LME_fosi_fished_',mod,'sigB_climate.csv']);
