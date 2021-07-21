YR=42;
ti = num2str(YR);
load(['/Volumes/MIP/GCM_DATA/CESM/FOSI/Data_cesm_fosi_daily_',ti,'.mat'],'ESM');

DY=365;
ENVR.Tp(:,1)  = ESM.Tp(param.ID,DY);
ENVR.Tb(:,1)  = ESM.Tb(param.ID,DY);
ENVR.Zm(:,1)  = ESM.Zm(param.ID,DY);
ENVR.det(:,1) = ESM.det(param.ID,DY);
ENVR.fZm(:,1) = zeros(param.NX,1);
ENVR.fB(:,1)  = zeros(param.NX,1);
ENVR.H(:,1)   = GRD.Z(param.ID);
ENVR.A(:,1)   = GRD.area(param.ID);

hist(ENVR.Tp)
hist(ENVR.Tb)
hist(ENVR.Zm)
hist(log10(ENVR.Zm))
hist(log10(ENVR.det))

zoo_loss = ESM.dZm(param.ID,DY);
hist(log10(zoo_loss))

Tfn = exp(-4000 .* ( (1./(ENVR.Tp+273.15)) - (1./303.15) ));
hist(Tfn)

Zprime = max((ENVR.Zm - 0.01),0);
hist(log10(Zprime))

Lzoo_quad = zoo_loss - (Tfn .* 0.1 .* Zprime);
%hist(log10(Lzoo_quad)) %Error using histc
hist((Lzoo_quad))
%hist(log10(Lzoo_quad)+eps) %Error using histc
 
iz = find(Lzoo_quad==0);
in = find(Lzoo_quad<0);
hist((Lzoo_quad(in)))

ENVR.dZm(:,1) = Lzoo_quad;
ENVR.dZm  = sub_neg(ENVR.dZm);

hist((zoo_loss(in)))

%% use full year from ESM
zoo_loss = ESM.dZm(param.ID,:);
Tfn = exp(-4000 .* ( (1./(ESM.Tp+273.15)) - (1./303.15) ));

plot(1:365,mean(Tfn)); hold on
plot(1:365,median(Tfn));
%%
Zprime = max((ESM.Zm - 0.01),0);
hist(Zprime(:))
%%
hist(log10(Zprime(:)))
%%
plot(1:365,mean(Zprime));
%%
plot(1:365,log10(mean(Zprime)));
%%
Lzoo_quad = zoo_loss - (Tfn .* 0.1 .* Zprime);
plot(1:365,zeros(365,1)); hold on;
plot(1:365,Lzoo_quad(in(1),:)); hold on;
plot(1:365,Lzoo_quad(in(2),:)); hold on;
plot(1:365,Lzoo_quad(in(3),:)); hold on;
plot(1:365,Lzoo_quad(in(150),:))

figure
plot(1:365,zeros(365,1)); hold on;
plot(1:365,Lzoo_quad(in(1),:)); hold on;
plot(1:365,ESM.dZm(in(2),:)); hold on;
plot(1:365,ESM.dZm(in(3),:)); hold on;
plot(1:365,ESM.dZm(in(150),:)); hold on;

%%
hist(ESM.dZm(:))

%%
clear ESM
YR=43;
ti = num2str(YR);
load(['/Volumes/MIP/GCM_DATA/CESM/FOSI/Data_cesm_fosi_daily_',ti,'.mat'],'ESM');




