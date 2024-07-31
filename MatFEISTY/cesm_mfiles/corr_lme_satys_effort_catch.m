% Calc corr of catch with effort
% calc only once, then put together with others
% min yrs as sat chl

clear
close all

%% FOSI input forcing

%cpath = '/Volumes/MIP/GCM_DATA/CESM/FOSI/';
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

load([cpath 'Data_grid_POP_gx1v6_noSeas.mat']);
ID = GRD.ID;

%% FEISTY outputs
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

%fpath=['/Volumes/MIP/NC/CESM_MAPP/' cfile '/'];
fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',cfile,'/corrs'];

mod = 'v15_obsfish';

%% Fish data - Catch
ypath='/Volumes/petrik-lab/Feisty/Fish-MIP/Phase3/fishing/';

% Anoms with linear trend removed
load([ypath 'FishMIP_Phase3a_LME_Catch_1948-2015_ann_mean_anoms.mat'],'aall',...
    'af','ap','ad')

aca = aall;
acf = af;
acp = ap;
acd = ad;

clear aall af ap ad

%% Effort
load([ypath 'FishMIP_Phase3a_LME_Effort_1948-2010_ann_mean_anoms.mat'],'aall',...
    'af','ap','ad')

aea = aall;
aef = af;
aep = ap;
aed = ad;

clear aall af ap ad

%% subset effort years
fyr = 1948:2010;
eyr = 1948:2015;
[yr,fid] = intersect(fyr,eyr);

aca    = aca(:,fid);
acf    = acf(:,fid);
acp    = acp(:,fid);
acd    = acd(:,fid);


%% Drivers from satellite obs
load([fpath 'lme_satellite_sst_chl_ann_mean_anoms.mat'])

%% match years
[~,cid] = intersect(yr,cyr);
[~,tid] = intersect(yr,tyr);

%% restrict analysis to only years with satellite chl data
aca    = aca(:,cid);
acf    = acf(:,cid);
acp    = acp(:,cid);
acd    = acd(:,cid);

aea    = aea(:,cid);
aef    = aef(:,cid);
aep    = aep(:,cid);
aed    = aed(:,cid);

%% %Corr of forcing ---------------------------------------------------------
cnam = {'corr','p','lag','idriver','driver'};

% All LMEs except inland seas (23=Baltic, 33=Red Sea, 62=Black Sea)
AA = ones(66,1);
iis = [23;33;62];
AA(iis) = nan;
lid = find(~isnan(AA));

%Lags
yr = 0;  %catch ~ effort that year only

% Drivers
FtabC = nan*ones(length(lid),length(yr));
FtabP = nan*ones(length(lid),length(yr));
PtabC = FtabC;
PtabP = FtabC;
DtabC = FtabC;
DtabP = FtabC;
AtabC = FtabC;
AtabP = FtabC;

%%
yst = 1;
yen = length(cid);

for L = 1:length(lid)

    %LME
    i = lid(L);
    ilme = num2str(i);

    t = 0;

    %                       LME   time           LME   time   
    %Fish
    [rp,pp] = corrcoef((aea(i,yst:yen-t))' , (aca(i,yst+t:yen))');
    AtabC(L) = rp(1,2);
    AtabP(L) = pp(1,2);
    clear rp pp

    [rp,pp] = corrcoef((aef(i,yst:yen-t))' , (acf(i,yst+t:yen))');
    FtabC(L) = rp(1,2);
    FtabP(L) = pp(1,2);
    clear rp pp

    [rp,pp] = corrcoef((aep(i,yst:yen-t))' , (acp(i,yst+t:yen))');
    PtabC(L) = rp(1,2);
    PtabP(L) = pp(1,2);
    clear rp pp

    [rp,pp] = corrcoef((aed(i,yst:yen-t))' , (acd(i,yst+t:yen))');
    DtabC(L) = rp(1,2);
    DtabP(L) = pp(1,2);
    clear rp pp

end %LME

%%
save([spath,'LMEs_corr_catch_satyrs_effort_nolags.mat'],...
    'FtabC','PtabC','DtabC','AtabC',...
    'FtabP','PtabP','DtabP','AtabP','lid');



