% Visualize output of FEISTY Climatology at single locations
% 150 years, monthly means saved

clear
close all

%%

datap = '/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/';

figp = '/Users/cpetrik/Dropbox/Princeton/FEISTY/CODE/Figs/PNG/CESM_MAPP/FOSI/';
fpath = [figp 'param_ensemble/'];
if (~isfolder(figp))
    mkdir(figp)
end

load('/Volumes/petrik-lab/Feisty/GCM_DATA/CESM/FOSI/cesm_grid_id_locs_area_dep.mat','ids','abbrev');
spots = abbrev;
ID = ids;
cols = {'bio','enc_f','enc_p','enc_d','enc_zm','enc_zl','enc_be','con_f',...
    'con_p','con_d','con_zm','con_zl','con_be','I','nu','gamma','die','rep',...
    'rec','clev','prod','pred','nmort','met','caught'};
cols=cols';
spots=spots';

encs = 50:10:70;
cmaxs = 5:5:20;

sname = 'FOSI_v15_locs_Spinup_All_fish03';

%%
SP_mbio=nan*ones(length(cmaxs),length(encs),length(ID));
SF_mbio=SP_mbio;
SD_mbio=SP_mbio;
MP_mbio=SP_mbio;
MF_mbio=SP_mbio;
MD_mbio=SP_mbio;
LP_mbio=SP_mbio;
LD_mbio=SP_mbio;
B_mbio =SP_mbio;

SP_mgr=SP_mbio;
SF_mgr=SP_mbio;
SD_mgr=SP_mbio;
MP_mgr=SP_mbio;
MF_mgr=SP_mbio;
MD_mgr=SP_mbio;
LP_mgr=SP_mbio;
LD_mgr=SP_mbio;

SP_lev=SP_mbio;
SF_lev=SP_mbio;
SD_lev=SP_mbio;
MP_lev=SP_mbio;
MF_lev=SP_mbio;
MD_lev=SP_mbio;
LP_lev=SP_mbio;
LD_lev=SP_mbio;

SP_prod=SP_mbio;
SF_prod=SP_mbio;
SD_prod=SP_mbio;
MP_prod=SP_mbio;
MF_prod=SP_mbio;
MD_prod=SP_mbio;
LP_prod=SP_mbio;
LD_prod=SP_mbio;

SP_gge=SP_mbio;
SF_gge=SP_mbio;
SD_gge=SP_mbio;
MP_gge=SP_mbio;
MF_gge=SP_mbio;
MD_gge=SP_mbio;
LP_gge=SP_mbio;
LD_gge=SP_mbio;

%% ------------------------ GROUP TOGETHER ----------------------------
for c=1%:length(cmaxs)
    for e=1%:length(encs)

        gam = encs(e);
        h = cmaxs(c);
        tcfn = num2str(h);
        tefn = num2str(round(gam));

        cfile = ['Dc_Lam700_enc',tefn,'-b200_m400-b175-k086_c',tcfn,...
            '-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100'];
        fpath=[datap cfile '/FOSI/'];
        sfile = [fpath sname];
        load(sfile);

        %%
        SP = S_Sml_p;
        SF = S_Sml_f;
        SD = S_Sml_d;
        MP = S_Med_p;
        MF = S_Med_f;
        MD = S_Med_d;
        LP = S_Lrg_p;
        LD = S_Lrg_d;
        CO = S_LTL;

        t=1:size(SP,1);
        % Last year
        lyr=t((end-12+1):end);
        % First year
        %lyr=1:12;

        %% Final mean biomass in each size

        SP_mbio(c,e,:)=squeeze(mean(SP(lyr,1,:)))';
        SF_mbio(c,e,:)=squeeze(mean(SF(lyr,1,:)))';
        SD_mbio(c,e,:)=squeeze(mean(SD(lyr,1,:)))';
        MP_mbio(c,e,:)=squeeze(mean(MP(lyr,1,:)))';
        MF_mbio(c,e,:)=squeeze(mean(MF(lyr,1,:)))';
        MD_mbio(c,e,:)=squeeze(mean(MD(lyr,1,:)))';
        LP_mbio(c,e,:)=squeeze(mean(LP(lyr,1,:)))';
        LD_mbio(c,e,:)=squeeze(mean(LD(lyr,1,:)))';
        B_mbio (c,e,:)=squeeze(mean(CO(lyr,1,:)))';

        %% Growth rate (nu - energy for biomass production)
        SP_mgr(c,e,:)=squeeze(nanmean(SP(lyr,15,:)))';
        SF_mgr(c,e,:)=squeeze(nanmean(SF(lyr,15,:)))';
        SD_mgr(c,e,:)=squeeze(nanmean(SD(lyr,15,:)))';
        MP_mgr(c,e,:)=squeeze(nanmean(MP(lyr,15,:)))';
        MF_mgr(c,e,:)=squeeze(nanmean(MF(lyr,15,:)))';
        MD_mgr(c,e,:)=squeeze(nanmean(MD(lyr,15,:)))';
        LP_mgr(c,e,:)=squeeze(nanmean(LP(lyr,15,:)))';
        LD_mgr(c,e,:)=squeeze(nanmean(LD(lyr,15,:)))';

        %% Feeding level = con/cmax
        SP_lev(c,e,:)=squeeze(nanmean(SP(lyr,20,:)))';
        SF_lev(c,e,:)=squeeze(nanmean(SF(lyr,20,:)))';
        SD_lev(c,e,:)=squeeze(nanmean(SD(lyr,20,:)))';
        MP_lev(c,e,:)=squeeze(nanmean(MP(lyr,20,:)))';
        MF_lev(c,e,:)=squeeze(nanmean(MF(lyr,20,:)))';
        MD_lev(c,e,:)=squeeze(nanmean(MD(lyr,20,:)))';
        LP_lev(c,e,:)=squeeze(nanmean(LP(lyr,20,:)))';
        LD_lev(c,e,:)=squeeze(nanmean(LD(lyr,20,:)))';

        %% Production (= nu * biom)
        SP_prod(c,e,:)=squeeze(nanmean(SP(lyr,21,:)))';
        SF_prod(c,e,:)=squeeze(nanmean(SF(lyr,21,:)))';
        SD_prod(c,e,:)=squeeze(nanmean(SD(lyr,21,:)))';
        MP_prod(c,e,:)=squeeze(nanmean(MP(lyr,21,:)))';
        MF_prod(c,e,:)=squeeze(nanmean(MF(lyr,21,:)))';
        MD_prod(c,e,:)=squeeze(nanmean(MD(lyr,21,:)))';
        LP_prod(c,e,:)=squeeze(nanmean(LP(lyr,21,:)))';
        LD_prod(c,e,:)=squeeze(nanmean(LD(lyr,21,:)))';

        %% Gross growth efficiency (= nu/consump)
        SP_gge(c,e,:)=squeeze(nanmean(SP(lyr,15,:)./SP(lyr,14,:)))';
        SF_gge(c,e,:)=squeeze(nanmean(SF(lyr,15,:)./SF(lyr,14,:)))';
        SD_gge(c,e,:)=squeeze(nanmean(SD(lyr,15,:)./SD(lyr,14,:)))';
        MP_gge(c,e,:)=squeeze(nanmean(MP(lyr,15,:)./MP(lyr,14,:)))';
        MF_gge(c,e,:)=squeeze(nanmean(MF(lyr,15,:)./MF(lyr,14,:)))';
        MD_gge(c,e,:)=squeeze(nanmean(MD(lyr,15,:)./MD(lyr,14,:)))';
        LP_gge(c,e,:)=squeeze(nanmean(LP(lyr,15,:)./LP(lyr,14,:)))';
        LD_gge(c,e,:)=squeeze(nanmean(LD(lyr,15,:)./LD(lyr,14,:)))';

    end
end

%% save

save([datap sname '_aenc_acmax_search.mat'],'spots','cols','ID',...
    'encs','cmaxs',...
    'SP_mbio','SD_mbio','MP_mbio','MF_mbio','MD_mbio','LP_mbio','LD_mbio','B_mbio ',...
    'SP_mgr','SF_mgr','SD_mgr','MP_mgr','MF_mgr','MD_mgr','LP_mgr','LD_mgr',...
    'SP_lev','SF_lev','SD_lev','MP_lev','MF_lev','MD_lev','LP_lev','LD_lev',...
    'SP_prod','SF_prod','SD_prod','MP_prod','MF_prod','MD_prod','LP_prod','LD_prod',...
    'SP_gge','SF_gge','SD_gge','MP_gge','MF_gge','MD_gge','LP_gge','LD_gge')

%% load if already saved
load([datap sname '_aenc_acmax_search.mat']);

%% ------------------------- PLOTS ------------------------------
% Biomass of each type
allF = SF_mbio+MF_mbio;
allP = SP_mbio+MP_mbio+LP_mbio;
allD = SD_mbio+MD_mbio+LD_mbio;
allB = B_mbio;

% Gut fullness
Fflev = (SF_lev+MF_lev)/2;
Pflev = (SP_lev+MP_lev+LP_lev)/3;
Dflev = (SD_lev+MD_lev+LD_lev)/3;
Sflev = (SF_lev+SP_lev+SD_lev)/3;
Mflev = (MF_lev+MP_lev+MD_lev)/3;
Lflev = (LP_lev+LD_lev)/2;
Tflev = (SF_lev+MF_lev+SP_lev+MP_lev+LP_lev+SD_lev+MD_lev+LD_lev)/7;

% GGE
Fgge = (SF_gge+MF_gge)/2;
Pgge = (SP_gge+MP_gge+LP_gge)/3;
Dgge = (SD_gge+MD_gge+LD_gge)/3;
Sgge = (SF_gge+SP_gge+SD_gge)/3;
Mgge = (MF_gge+MP_gge+MD_gge)/3;
Lgge = (LP_gge+LD_gge)/2;
Tgge = (SF_gge+MF_gge+SP_gge+MP_gge+LP_gge+SD_gge+MD_gge+LD_gge)/7;

%%
nk = length(encs);
nj = length(cmaxs);

FPrat = squeeze(allF./(allF+allP));
DPrat = squeeze(allD./(allD+allP));

ays = [encs (encs(end)+10)];
jays = [cmaxs (cmaxs(end)+5)];
[agrid,jgrid]=meshgrid(ays,jays);

allF2 = NaN*ones(nj+1,nk+1,14);
allP2 = NaN*ones(nj+1,nk+1,14);
allD2 = NaN*ones(nj+1,nk+1,14);
FP2 = NaN*ones(nj+1,nk+1,14);
DP2 = NaN*ones(nj+1,nk+1,14);
Tgge2 = NaN*ones(nj+1,nk+1,14);
Tflev2 = NaN*ones(nj+1,nk+1,14);
Sgge2 = NaN*ones(nj+1,nk+1,14);
Sflev2 = NaN*ones(nj+1,nk+1,14);
Mgge2 = NaN*ones(nj+1,nk+1,14);
Mflev2 = NaN*ones(nj+1,nk+1,14);
Lgge2 = NaN*ones(nj+1,nk+1,14);
Lflev2 = NaN*ones(nj+1,nk+1,14);

allF2(1:nj,1:nk,:)=allF;
allP2(1:nj,1:nk,:)=allP;
allD2(1:nj,1:nk,:)=allD;
FP2(1:nj,1:nk,:)=FPrat;
DP2(1:nj,1:nk,:)=DPrat;
Tgge2(1:nj,1:nk,:)=Tgge;
Tflev2(1:nj,1:nk,:)=Tflev;
Sgge2(1:nj,1:nk,:)=Sgge;
Sflev2(1:nj,1:nk,:)=Sflev;
Mgge2(1:nj,1:nk,:)=Mgge;
Mflev2(1:nj,1:nk,:)=Mflev;
Lgge2(1:nj,1:nk,:)=Lgge;
Lflev2(1:nj,1:nk,:)=Lflev;

%%
% colors
cmBP=cbrewer('seq','BuPu',50,'PCHIP');
cmYOR=cbrewer('seq','YlOrRd',50,'PCHIP');

%% Only use 3 domain examples: EBS(8), PUp(14), HOT(11)
sid = [8;14;11];
for s=1:length(sid)
    domain = sid(s);
    loc = spots{domain};
    lname = [loc '_'];
    
    %% Sum mean biom over stages
    f1=figure(1);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(log10(allF2(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmBP)
    caxis([-2 2])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'log10 Mean F Biom (g m^-^2)'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_C')
    end
    
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(log10(allP2(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmBP)
    caxis([-2 2])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('log10 Mean P Biom (g m^-^2)')
    end
    if (s==1)
        ylabel('a_C')
    end
    
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(log10(allD2(:,:,domain))))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmBP)
    caxis([-2 2])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('log10 Mean D Biom (g m^-^2)')
    end
    xlabel('a_E')
    if (s==1)
        ylabel('a_C')
    end
    %stamp(sname)
    
    
    f2=figure(2);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(FP2(:,:,domain)))
    cmocean('balance')
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'Fraction F/(F+P)'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_C')
    end
    
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(DP2(:,:,domain)))
    colorbar('Position',[0.92 0.5 0.025 0.3],'orientation','vertical')
    cmocean('balance')
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('Fraction D/(D+P)')
    end
    xlabel('a_E')
    if (s==1)
        ylabel('a_C')
    end
    %stamp(sname)
    
    %% Feeding level
    f3=figure(3);
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(Tflev2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0.5 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'Mean feeding level'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_C')
    end
    
    % GGE
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(Tgge2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('Mean gross growth efficiency')
    end
    xlabel('a_E')
    if (s==1)
        ylabel('a_C')
    end
    %stamp(sname)
    
    %% Feeding level only
    f4=figure(4);
    % Small
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(Sflev2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'S mean feeding level'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_C')
    end
    
    % Medium
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(Mflev2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('M mean feeding level')
    end
    if (s==1)
        ylabel('a_C')
    end
    
    % Large
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(Lflev2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('L mean feeding level')
    end
    xlabel('a_E')
    if (s==1)
        ylabel('a_C')
    end
    %stamp(sname)
    
    %% GGE only
    f5=figure(5);
    % Small
    subplot(3,3,s)
    pcolor(agrid,jgrid,squeeze(Sgge2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title({loc; 'S mean gross growth efficiency'})
    else
        title({loc; ''})
    end
    if (s==1)
        ylabel('a_C')
    end
    
    % Medium
    subplot(3,3,s+3)
    pcolor(agrid,jgrid,squeeze(Mgge2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('M mean gross growth efficiency')
    end
    if (s==1)
        ylabel('a_C')
    end
    
    % Large
    subplot(3,3,s+6)
    pcolor(agrid,jgrid,squeeze(Lgge2(:,:,domain)))
    colorbar('Position',[0.92 0.35 0.025 0.3],'orientation','vertical')
    colormap(cmYOR)
    caxis([0 1])
    set(gca,'XTick',(10:20:100)+5,'XTickLabel',10:20:100,...
        'YTick',(10:20:100)+5,'YTickLabel',10:20:100)
    if (s==2)
        title('L mean gross growth efficiency')
    end
    xlabel('a_E')
    if (s==1)
        ylabel('a_C')
    end
    %stamp(sname)
    
end %spots
%%
print(f1,'-dpng',[figp sname '_mean_biomass_type_all_3locs.png'])
print(f2,'-dpng',[figp sname '_frac_all_3locs.png'])
print(f3,'-dpng',[figp sname '_flev_gge_all_3locs.png'])
print(f4,'-dpng',[figp sname '_flev_size_all_3locs.png'])
print(f5,'-dpng',[figp sname '_gge_size_all_3locs.png'])



