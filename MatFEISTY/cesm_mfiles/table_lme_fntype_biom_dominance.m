% Use fish rel biomass to define 3-4 ecosys/foodweb structures
% For all 63 LMEs

clear
close all

%% % ------------------------------------------------------------
cpath = '/Volumes/petrik-lab/Feisty/GCM_Data/CESM/FOSI/';

cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath = ['/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/' cfile '/corrs/'];

mod = 'v15_All_fish03_';

%% inputs
load([cpath 'lme_means_g.e11_LENS.GECOIAF.T62_g16.009.mat'],...
    'lme_tp_fosi','lme_tb_fosi','lme_det_fosi','lme_mz_fosi','lme_mzloss_fosi')

% fish
load([fpath 'LME_fosi_fished_',mod,cfile '.mat'],'lme_area','lme_mtype',...
    'lme_mnu');

lme_nu = lme_mnu;


%% Relative amounts
lme_type = lme_mtype(:,1:3);
lme_btot = sum(lme_type,2);

lbio = lme_type ./ repmat(lme_btot,1,3);

%% Define foodweb types
% F dom: F>=0.4     %1
% D&F dom: P<0.2    %2
% D dom: D>=0.5     %3
% P&F dom: D<0.2    %4
% even              %5

% Forage fish dominated regions and large pelagic fish dominated regions 
% were determined as LMEs where demersals represented <20% of total fish 
% biomass. Demersal dominated regions were LMEs where demersals accounted 
% for >50% of total fish biomass. 

etype = 5*ones(66,1);
etype(lbio(:,1)>=0.4) = 1;
etype(lbio(:,2)<0.2) = 2;
etype(lbio(:,3)>=0.5) = 3;
etype(lbio(:,3)<0.2) = 4;

%NaNs are 23, 33, 62 (inland seas)
etype([23, 33, 62],1) = nan;

alltex = {'Forage',...    % = 1
    'F & D',...          % = 2
    'Demersal',...           % = 3
    'F & P',...           % = 4
    'Even'};             % = 5

fid = (etype==1);
fdid = (etype==2);
did = (etype==3);
fpid = (etype==4);
eid = (etype==5);

etex = cell(66,1);
etex(fid) = {'Forage'};
etex(fdid) = {'F & D'};
etex(did) = {'Demersal'};
etex(fpid) = {'F & P'};
etex(eid) = {'Even'};
etex{23} = '';
etex{33} = '';
etex{62} = '';

et2 = string(etex);

lid = [1:66]';

%% Save table
Dtab = table(lid,et2,'VariableNames',{'LME','Dominance'});
writetable(Dtab,[fpath 'LME_dominance.csv'],'Delimiter',',',...
    'WriteRowNames',false)