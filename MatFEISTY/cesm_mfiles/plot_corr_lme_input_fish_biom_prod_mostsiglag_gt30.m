% Plot max corr of driver-fish corrs
% Multipanel barplot of biomass and prod
% For all 63 LMEs
% Only plot fn type if significant and biomass > 30%

clear
close all

%% % ------------------------------------------------------------
cfile = 'Dc_Lam700_enc70-b200_m400-b175-k086_c20-b250_D075_A050_sMZ090_mMZ045_nmort1_BE08_CC80_RE00100';

fpath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/FOSI/'];
spath=['/Volumes/petrik-lab/Feisty/NC/CESM_MAPP/' cfile '/regressions/'];
ppath=['/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Figs/CESM_MAPP/FOSI/',cfile,'/corrs/'];

mod = 'v15_All_fish03_';

%spath = '/Users/cpetrik/Petrik Lab Group Dropbox/Colleen Petrik/Princeton/FEISTY/CODE/Data/FOSI/';

%% Biomass
load([spath,'LMEs_corr_driver_maxcorrs.mat'],...
    'LAtab','LFtab','LPtab','LDtab','lid')

BAtab = LAtab;
BFtab = LFtab;
BPtab = LPtab;
BDtab = LDtab;

clear LAtab LFtab LPtab LDtab

%% Prod
load([spath,'LMEs_nu_corr_driver_maxcorrs.mat'],...
    'LAtab','LFtab','LPtab','LDtab')

PAtab = LAtab;
PFtab = LFtab;
PPtab = LPtab;
PDtab = LDtab;

clear LAtab LFtab LPtab LDtab

%% LMEs with >30% biomass
load([fpath 'LME_fosi_fished_',mod,cfile '.mat'],'lme_area','lme_mtype')

lme_type = lme_mtype(:,1:3);
lme_btot = sum(lme_type,2);

lbio = lme_type ./ repmat(lme_btot,1,3);

%remove NaNs are 23, 33, 62 (inland seas)
iis = setdiff(1:66,[23, 33, 62]);
% Relative biomass
rbio= lbio(iis,:);

%%  ---------------------------------------------------------
cnam = {'coef','p','lag','idriver','driver'};
ctex = {'TP','TB','Det','ZmLoss'};

load('paul_tol_cmaps.mat')

mcol(1,:) = drainbow(12,:)/255; % orange
mcol(2,:) = drainbow(4,:)/255; %dk blue
mcol(3,:) = drainbow(15,:)/255; %grey
mcol(4,:) = drainbow(6,:)/255; %lt blue

%% All
f1 = figure('Units','inches','Position',[1 3 7.5 10]);

%Bottom Left
subplot('Position',[0.075 0.075 0.425 0.2])
% Get fake colors first for legend
for i=1:length(ctex)
    L=lid(i);
    b=bar(L,BDtab(i,1),'EdgeColor','none','FaceColor',mcol(i,:));
    hold on
end
%legend of colors and shapes
lgd = legend(ctex,'Position',[0.25 0.01 0.5 0.025],'Orientation','horizontal');
lgd.AutoUpdate = 'off';
for i=1:length(lid)
    L=lid(i);
    if (BDtab(i,2) <= 0.05)
        b=bar(L,BDtab(i,1),'EdgeColor',mcol(BDtab(i,4),:),'FaceColor',mcol(BDtab(i,4),:));
        hold on
    else
        b=bar(L,BDtab(i,1),'EdgeColor',mcol(BDtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(BDtab(i,3));
    if (BDtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66)
ylabel('Demersal')

%Top Left
subplot('Position',[0.075 0.775 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if (BAtab(i,2) <= 0.05)
        b=bar(L,BAtab(i,1),'EdgeColor',mcol(BAtab(i,4),:),'FaceColor',mcol(BAtab(i,4),:));
        hold on
    else
        b=bar(L,BAtab(i,1),'EdgeColor',mcol(BAtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(BAtab(i,3));
    if (BAtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end

end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('All fishes')
title('Biomass')

%2nd L
subplot('Position',[0.075 0.545 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if (BFtab(i,2) <= 0.05)
        b=bar(L,BFtab(i,1),'EdgeColor',mcol(BFtab(i,4),:),'FaceColor',mcol(BFtab(i,4),:));
        hold on
    else
        b=bar(L,BFtab(i,1),'EdgeColor',mcol(BFtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(BFtab(i,3));
    if (BFtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('Forage')

%3rd L
subplot('Position',[0.075 0.315 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if (BPtab(i,2) <= 0.05)
        b=bar(L,BPtab(i,1),'EdgeColor',mcol(BPtab(i,4),:),'FaceColor',mcol(BPtab(i,4),:));
        hold on
    else
        b=bar(L,BPtab(i,1),'EdgeColor',mcol(BPtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(BPtab(i,3));
    if (BPtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('Lg Pel')


% Top Right
subplot('Position',[0.525 0.775 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if (PAtab(i,2) <= 0.05)
        b=bar(L,PAtab(i,1),'EdgeColor',mcol(PAtab(i,4),:),'FaceColor',mcol(PAtab(i,4),:));
        hold on
    else
        b=bar(L,PAtab(i,1),'EdgeColor',mcol(PAtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(PAtab(i,3));
    if (PAtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end

end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','','YTickLabel','')
%ylabel('All fishes')
title('Production')

%2nd R
subplot('Position',[0.525 0.545 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if (PFtab(i,2) <= 0.05)
        b=bar(L,PFtab(i,1),'EdgeColor',mcol(PFtab(i,4),:),'FaceColor',mcol(PFtab(i,4),:));
        hold on
    else
        b=bar(L,PFtab(i,1),'EdgeColor',mcol(PFtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(PFtab(i,3));
    if (PFtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','','YTickLabel','')
%ylabel('Forage')

%3rd R
subplot('Position',[0.525 0.315 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if (PPtab(i,2) <= 0.05)
        b=bar(L,PPtab(i,1),'EdgeColor',mcol(PPtab(i,4),:),'FaceColor',mcol(PPtab(i,4),:));
        hold on
    else
        b=bar(L,PPtab(i,1),'EdgeColor',mcol(PPtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(PPtab(i,3));
    if (PPtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','','YTickLabel','')
%ylabel('Lg Pel')

%Bottom R
subplot('Position',[0.525 0.075 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if (PDtab(i,2) <= 0.05)
        b=bar(L,PDtab(i,1),'EdgeColor',mcol(PDtab(i,4),:),'FaceColor',mcol(PDtab(i,4),:));
        hold on
    else
        b=bar(L,PDtab(i,1),'EdgeColor',mcol(PDtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(PDtab(i,3));
    if (PDtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66,'YTickLabel','')
%ylabel('Demersal')

%stamp('')
print('-dpng',[ppath 'Bar_LMEs_biom_nu_maxcorr_fntypes.png'])

%% Only biom > 30%
f2 = figure('Units','inches','Position',[1 3 7.5 10]);

%Bottom Left
subplot('Position',[0.075 0.075 0.425 0.2])
% Get fake colors first for legend
dc=[46;9;1;6];
for i=1:length(dc) %1:length(ctex)
    di = dc(i);
    L=lid(di);
    if(rbio(di,3)>0.3)
        b=bar(L,BDtab(di,1),'EdgeColor','none','FaceColor',mcol(BDtab(di,4),:));
        hold on
    end
end
%legend of colors and shapes
lgd = legend(ctex,'Position',[0.25 0.01 0.5 0.025],'Orientation','horizontal');
lgd.AutoUpdate = 'off';
for i=1:length(lid)
    L=lid(i);
    if(rbio(i,3)>0.3)
        if (BDtab(i,2) <= 0.05)
            b=bar(L,BDtab(i,1),'EdgeColor',mcol(BDtab(i,4),:),'FaceColor',mcol(BDtab(i,4),:));
            hold on
        else
            b=bar(L,BDtab(i,1),'EdgeColor',mcol(BDtab(i,4),:),'FaceColor',[1 1 1]);
            hold on
        end
        xtips1 = b(1).XEndPoints;
        ytips1 = b(1).YEndPoints;
        labels = num2str(BDtab(i,3));
        if (BDtab(i,1) < 0)
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','top')
        else
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom')
        end
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66)
ylabel('Demersal')

%Top Left
subplot('Position',[0.075 0.775 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if (BAtab(i,2) <= 0.05)
        b=bar(L,BAtab(i,1),'EdgeColor',mcol(BAtab(i,4),:),'FaceColor',mcol(BAtab(i,4),:));
        hold on
    else
        b=bar(L,BAtab(i,1),'EdgeColor',mcol(BAtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(BAtab(i,3));
    if (BAtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end

end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('All fishes')
title('Biomass')

%2nd L
subplot('Position',[0.075 0.545 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if(rbio(i,1)>0.3)
        if (BFtab(i,2) <= 0.05)
            b=bar(L,BFtab(i,1),'EdgeColor',mcol(BFtab(i,4),:),'FaceColor',mcol(BFtab(i,4),:));
            hold on
        else
            b=bar(L,BFtab(i,1),'EdgeColor',mcol(BFtab(i,4),:),'FaceColor',[1 1 1]);
            hold on
        end
        xtips1 = b(1).XEndPoints;
        ytips1 = b(1).YEndPoints;
        labels = num2str(BFtab(i,3));
        if (BFtab(i,1) < 0)
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','top')
        else
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom')
        end
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('Forage')

%3rd L
subplot('Position',[0.075 0.315 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if(rbio(i,2)>0.3)
        if (BPtab(i,2) <= 0.05)
            b=bar(L,BPtab(i,1),'EdgeColor',mcol(BPtab(i,4),:),'FaceColor',mcol(BPtab(i,4),:));
            hold on
        else
            b=bar(L,BPtab(i,1),'EdgeColor',mcol(BPtab(i,4),:),'FaceColor',[1 1 1]);
            hold on
        end
        xtips1 = b(1).XEndPoints;
        ytips1 = b(1).YEndPoints;
        labels = num2str(BPtab(i,3));
        if (BPtab(i,1) < 0)
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','top')
        else
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom')
        end
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','')
ylabel('Lg Pel')


% Top Right
subplot('Position',[0.525 0.775 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if (PAtab(i,2) <= 0.05)
        b=bar(L,PAtab(i,1),'EdgeColor',mcol(PAtab(i,4),:),'FaceColor',mcol(PAtab(i,4),:));
        hold on
    else
        b=bar(L,PAtab(i,1),'EdgeColor',mcol(PAtab(i,4),:),'FaceColor',[1 1 1]);
        hold on
    end
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels = num2str(PAtab(i,3));
    if (PAtab(i,1) < 0)
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','top')
    else
        text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end

end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','','YTickLabel','')
%ylabel('All fishes')
title('Production')

%2nd R
subplot('Position',[0.525 0.545 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if(rbio(i,1)>0.3)
        if (PFtab(i,2) <= 0.05)
            b=bar(L,PFtab(i,1),'EdgeColor',mcol(PFtab(i,4),:),'FaceColor',mcol(PFtab(i,4),:));
            hold on
        else
            b=bar(L,PFtab(i,1),'EdgeColor',mcol(PFtab(i,4),:),'FaceColor',[1 1 1]);
            hold on
        end
        xtips1 = b(1).XEndPoints;
        ytips1 = b(1).YEndPoints;
        labels = num2str(PFtab(i,3));
        if (PFtab(i,1) < 0)
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','top')
        else
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom')
        end
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','','YTickLabel','')
%ylabel('Forage')

%3rd R
subplot('Position',[0.525 0.315 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if(rbio(i,2)>0.3)
        if (PPtab(i,2) <= 0.05)
            b=bar(L,PPtab(i,1),'EdgeColor',mcol(PPtab(i,4),:),'FaceColor',mcol(PPtab(i,4),:));
            hold on
        else
            b=bar(L,PPtab(i,1),'EdgeColor',mcol(PPtab(i,4),:),'FaceColor',[1 1 1]);
            hold on
        end
        xtips1 = b(1).XEndPoints;
        ytips1 = b(1).YEndPoints;
        labels = num2str(PPtab(i,3));
        if (PPtab(i,1) < 0)
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','top')
        else
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom')
        end
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel','','YTickLabel','')
%ylabel('Lg Pel')

%Bottom R
subplot('Position',[0.525 0.075 0.425 0.2])
for i=1:length(lid)
    L=lid(i);
    if(rbio(i,3)>0.3)
        if (PDtab(i,2) <= 0.05)
            b=bar(L,PDtab(i,1),'EdgeColor',mcol(PDtab(i,4),:),'FaceColor',mcol(PDtab(i,4),:));
            hold on
        else
            b=bar(L,PDtab(i,1),'EdgeColor',mcol(PDtab(i,4),:),'FaceColor',[1 1 1]);
            hold on
        end
        xtips1 = b(1).XEndPoints;
        ytips1 = b(1).YEndPoints;
        labels = num2str(PDtab(i,3));
        if (PDtab(i,1) < 0)
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','top')
        else
            text(xtips1,ytips1,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom')
        end
    end
end
xlim([0 67])
ylim([-1.1 1.1])
set(gca,'XTick',1:3:66,'XTickLabel',1:3:66,'YTickLabel','')
%ylabel('Demersal')

%stamp('')

print('-dpng',[ppath 'Bar_LMEs_biom_nu_maxcorr_fntypes_gt30.png'])


