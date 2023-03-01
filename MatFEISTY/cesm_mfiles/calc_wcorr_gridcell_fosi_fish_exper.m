% Map correlations of FEISTY FOSI experiments
% with full var

clear
close all

%% FOSI input forcing
%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% Fish data
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

load([fpath 'FEISTY_FOSI_',mod,'ann_mean_anoms.mat'],'as','am','al','af','ap','ad','aa','ab'); % Anoms with linear trend removed
vas = as;
vam = am;
val = al;
vaf = af;
vap = ap;
vad = ad;
vaa = aa;
vab = ab;

clear as am al af ap ad aa ab

%% correlations by grid cell
rS = nan*ones(length(ID),3);
rM = rS;
rL = rS;
rF = rS;
rP = rS;
rD = rS;
rA = rS;
rB = rS;
pS = rS;
pM = pS;
pL = pS;
pF = pS;
pP = pS;
pD = pS;
pA = pS;
pB = pS;

%% Corr

% Loop over sims
for j=2:4
    mod = sims{j};

    load([fpath 'FEISTY_FOSI_',mod,'ann_mean_anoms.mat'],'as','am','al','af','ap','ad','aa','ab'); % Anoms with linear trend removed


    %Fish
    [rho,pval] = (corr((as(1:28604,:))',(vas(1:28604,:))'));
    rS1 = diag(rho);
    p1  = diag(pval);
    [rho,pval] = (corr((as(28605:57208,:))',(vas(28605:57208,:))'));
    rS2 = diag(rho);
    p2  = diag(pval);
    [rho,pval] = (corr((as(57209:end,:))',(vas(57209:end,:))'));
    rS3 = diag(rho);
    p3  = diag(pval);
    rS(:,j) = [rS1;rS2;rS3];
    pS(:,j) = [p1;p2;p3];
    clear rS1 rS2 rS3 p1 p2 p3

    [rho,pval] = (corr((am(1:28604,:))',(vam(1:28604,:))'));
    rS1 = diag(rho);
    p1  = diag(pval);
    [rho,pval] = (corr((am(28605:57208,:))',(vam(28605:57208,:))'));
    rS2 = diag(rho);
    p2  = diag(pval);
    [rho,pval] = (corr((am(57209:end,:))',(vam(57209:end,:))'));
    rS3 = diag(rho);
    p3  = diag(pval);
    rM(:,j) = [rS1;rS2;rS3];
    pM(:,j) = [p1;p2;p3];
    clear rS1 rS2 rS3 p1 p2 p3

    [rho,pval] = (corr((al(1:28604,:))',(val(1:28604,:))'));
    rS1 = diag(rho);
    p1  = diag(pval);
    [rho,pval] = (corr((al(28605:57208,:))',(val(28605:57208,:))'));
    rS2 = diag(rho);
    p2  = diag(pval);
    [rho,pval] = (corr((al(57209:end,:))',(val(57209:end,:))'));
    rS3 = diag(rho);
    p3  = diag(pval);
    rL(:,j) = [rS1;rS2;rS3];
    pL(:,j) = [p1;p2;p3];
    clear rS1 rS2 rS3 p1 p2 p3

    [rho,pval] = (corr((af(1:28604,:))',(vaf(1:28604,:))'));
    rS1 = diag(rho);
    p1  = diag(pval);
    [rho,pval] = (corr((af(28605:57208,:))',(vaf(28605:57208,:))'));
    rS2 = diag(rho);
    p2  = diag(pval);
    [rho,pval] = (corr((af(57209:end,:))',(vaf(57209:end,:))'));
    rS3 = diag(rho);
    p3  = diag(pval);
    rF(:,j) = [rS1;rS2;rS3];
    pF(:,j) = [p1;p2;p3];
    clear rS1 rS2 rS3 p1 p2 p3

    [rho,pval] = (corr((ap(1:28604,:))',(vap(1:28604,:))'));
    rS1 = diag(rho);
    p1  = diag(pval);
    [rho,pval] = (corr((ap(28605:57208,:))',(vap(28605:57208,:))'));
    rS2 = diag(rho);
    p2  = diag(pval);
    [rho,pval] = (corr((ap(57209:end,:))',(vap(57209:end,:))'));
    rS3 = diag(rho);
    p3  = diag(pval);
    rP(:,j) = [rS1;rS2;rS3];
    pP(:,j) = [p1;p2;p3];
    clear rS1 rS2 rS3 p1 p2 p3

    [rho,pval] = (corr((ad(1:28604,:))',(vad(1:28604,:))'));
    rS1 = diag(rho);
    p1  = diag(pval);
    [rho,pval] = (corr((ad(28605:57208,:))',(vad(28605:57208,:))'));
    rS2 = diag(rho);
    p2  = diag(pval);
    [rho,pval] = (corr((ad(57209:end,:))',(vad(57209:end,:))'));
    rS3 = diag(rho);
    p3  = diag(pval);
    rD(:,j) = [rS1;rS2;rS3];
    pD(:,j) = [p1;p2;p3];
    clear rS1 rS2 rS3 p1 p2 p3

    [rho,pval] = (corr((aa(1:28604,:))',(vaa(1:28604,:))'));
    rS1 = diag(rho);
    p1  = diag(pval);
    [rho,pval] = (corr((aa(28605:57208,:))',(vaa(28605:57208,:))'));
    rS2 = diag(rho);
    p2  = diag(pval);
    [rho,pval] = (corr((aa(57209:end,:))',(vaa(57209:end,:))'));
    rS3 = diag(rho);
    p3  = diag(pval);
    rA(:,j) = [rS1;rS2;rS3];
    pA(:,j) = [p1;p2;p3];
    clear rS1 rS2 rS3 p1 p2 p3

    [rho,pval] = (corr((ab(1:28604,:))',(vab(1:28604,:))'));
    rS1 = diag(rho);
    p1  = diag(pval);
    [rho,pval] = (corr((ab(28605:57208,:))',(vab(28605:57208,:))'));
    rS2 = diag(rho);
    p2  = diag(pval);
    [rho,pval] = (corr((ab(57209:end,:))',(vab(57209:end,:))'));
    rS3 = diag(rho);
    p3  = diag(pval);
    rB(:,j) = [rS1;rS2;rS3];
    pB(:,j) = [p1;p2;p3];
    clear rS1 rS2 rS3 p1 p2 p3

end

%%
save([fpath 'grid_corrs_fish_FOSI_fished_expers.mat'],...
    'rS','rM','rL','rF','rP','rD','rA','rB',...
    'pS','pM','pL','pF','pP','pD','pA','pB','sims');

