% this file plots the figures
clear all;
currdir = pwd;
restoredefaultpath;
% replicate Figure 1 in KM (2012)

fixed_rf=0;
label_R = 'cmy';
prior_type = 'minnesota';
ndraws  =1e6;
prior_only=0;
spec= ['rgibbs_',label_R,'prior_only_',num2str(prior_only),'prior_',prior_type,'fixed_rf_',num2str(fixed_rf),'_ndraws',num2str(ndraws)];


%L =load('matfiles/rgibbs_cmyprior_only_0prior_minnesotafixed_rf_0_ndraws100.mat')

L =load(['matfiles/',filesep,spec,'.mat']);




nhorizons = 21;%size(L.L,1);
nvar      = size(L.L,2);
nshocks   = size(L.L,3);



L_q50=nan(nhorizons,nvar,nshocks); % store conjF quantile 50th
L_q16=nan(nhorizons,nvar,nshocks); % store conjF quantile 16th
L_q84=nan(nhorizons,nvar,nshocks); % store conjF quantile 84th


for ii=1:nhorizons
    for jj=1:nvar
        for kk=1:nshocks
            


        if jj==15

            L_q50(ii,15,kk) = quantile(L.L(ii,15,kk,:)-L.L(ii,11,kk,:),0.5); % real wage 15 - 11
            L_q16(ii,15,kk) = quantile(L.L(ii,15,kk,:)-L.L(ii,11,kk,:),0.16);
            L_q84(ii,15,kk) = quantile(L.L(ii,15,kk,:)-L.L(ii,11,kk,:),0.84);

        else


            L_q50(ii,jj,kk) = quantile(L.L(ii,jj,kk,:),0.5);
            L_q16(ii,jj,kk) = quantile(L.L(ii,jj,kk,:),0.16);
            L_q84(ii,jj,kk) = quantile(L.L(ii,jj,kk,:),0.84);

        end

        end
    end
end




% figure settings
ftsizeaxis   = 11;
ftsizexlabel = 11;
ftsizetitle  = 11;
ftsizelegend = 11;
ftlinewidth  = 1.0;
medianwidth  = 1.0;


plot_uv_scale = 1;
horizon       = nhorizons-1;

if ~exist('pngfiles', 'dir')
    mkdir('pngfiles');
    disp('Created directory: pngfiles');
else
    disp('Directory pngfiles already exists');
end


close all
%% figure 1: demand

hFig = figure('name','Demand');
set(hFig, 'Position', [20 20 900 450])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'};
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(2,3,i)
variable   = demand_variables_index(i);
shock      = 1;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;

[~,~]=jbfill([0:horizon],auv16,buv84,rgb('crimson'),rgb('crimson'),0,0.5);

hold on
%hq16=plot([0:horizon],auv16,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
%hold on
%hq84=plot([0:horizon],buv84,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Quarters','FontSize',ftsizexlabel,'Interpreter','latex')
if i==6
yticks([-0.5 -0.25 0 0.25 0.5])
end
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i),'Interpreter','latex')
box off
legend off
grid on

end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_demand_slides.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_demand_slides.png'],'-dpng');





%% figure 2: investment

hFig = figure('name','Invesment');
set(hFig, 'Position', [20 20 600 300])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(2,3,i)
variable   = demand_variables_index(i);
shock      = 2;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;

[~,~]=jbfill([0:horizon],auv16,buv84,rgb('crimson'),rgb('crimson'),0,0.5);

hold on
%hq16=plot([0:horizon],auv16,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
%hold on
%hq84=plot([0:horizon],buv84,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Quarters','FontSize',ftsizexlabel,'Interpreter','latex')
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i),'Interpreter','latex')
box off
legend off
grid on

end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_investment.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_investment.png'],'-dpng');






%% figure 3: financial

hFig = figure('name','Financial');
set(hFig, 'Position', [20 20 700/2 900/2])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 3;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;

[~,~]=jbfill([0:horizon],auv16,buv84,rgb('crimson'),rgb('crimson'),0,0.5);

hold on
%hq16=plot([0:horizon],auv16,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
%hold on
%hq84=plot([0:horizon],buv84,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Quarters','FontSize',ftsizexlabel,'Interpreter','latex')
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i),'Interpreter','latex')
box off
legend off
grid on

end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_financial.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_financial.png'],'-dpng');








%% figure 4: monetary

hFig = figure('name','Monetary');
set(hFig, 'Position', [20 20 700 900])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 4;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;

[~,~]=jbfill([0:horizon],auv16,buv84,rgb('crimson'),rgb('crimson'),0,0.5);

hold on
%hq16=plot([0:horizon],auv16,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
%hold on
%hq84=plot([0:horizon],buv84,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Quarters','FontSize',ftsizexlabel,'Interpreter','latex')
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i),'Interpreter','latex')
box off
legend off
grid on

end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_monetary.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_monetary.png'],'-dpng');






%% figure 5: government Spending

hFig = figure('name','Government Spending');
set(hFig, 'Position', [20 20 700 900])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 5;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;

[~,~]=jbfill([0:horizon],auv16,buv84,rgb('crimson'),rgb('crimson'),0,0.5);

hold on
%hq16=plot([0:horizon],auv16,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
%hold on
%hq84=plot([0:horizon],buv84,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Quarters','FontSize',ftsizexlabel,'Interpreter','latex')
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i),'Interpreter','latex')
box off
legend off
grid on

end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_Gspending.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_Gspending.png'],'-dpng');






%% figure 6: technology

hFig = figure('name','technology');
set(hFig, 'Position', [20 20 700 900])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 6;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;

[~,~]=jbfill([0:horizon],auv16,buv84,rgb('crimson'),rgb('crimson'),0,0.5);

hold on
%hq16=plot([0:horizon],auv16,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
%hold on
%hq84=plot([0:horizon],buv84,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Quarters','FontSize',ftsizexlabel,'Interpreter','latex')
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i),'Interpreter','latex')
box off
legend off
grid on

end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_technology.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_technology.png'],'-dpng');







%% figure 6: technology

hFig = figure('name','technology');
set(hFig, 'Position', [20 20 700 900])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 6;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;

[~,~]=jbfill([0:horizon],auv16,buv84,rgb('crimson'),rgb('crimson'),0,0.5);

hold on
%hq16=plot([0:horizon],auv16,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
%hold on
%hq84=plot([0:horizon],buv84,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Quarters','FontSize',ftsizexlabel,'Interpreter','latex')
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i),'Interpreter','latex')
box off
legend off
grid on

end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_technology.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_technology.png'],'-dpng');






%% figure 7: laborsupply

hFig = figure('name','laborsupply');
set(hFig, 'Position', [20 20 700 900])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 7;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;

[~,~]=jbfill([0:horizon],auv16,buv84,rgb('crimson'),rgb('crimson'),0,0.5);

hold on
%hq16=plot([0:horizon],auv16,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
%hold on
%hq84=plot([0:horizon],buv84,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Quarters','FontSize',ftsizexlabel,'Interpreter','latex')
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i),'Interpreter','latex')
box off
legend off
grid on

end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_laborsupply.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_laborsupply.png'],'-dpng');








%% figure 7: wagebargaining

hFig = figure('name','wagebargaining');
set(hFig, 'Position', [20 20 700 900])

demand_variables_labels = {'Real GDP','PCE Price Index','Fed Funds Rate','Nonresidential Investment','Unemployment','Real Wage'}
demand_variables_index = [1,11,25,4,19,15];

for i=1:length(demand_variables_index)
subplot(3,2,i)
variable   = demand_variables_index(i);
shock      = 8;
auv16=L_q16(:,variable,shock)'*plot_uv_scale;
auv50=L_q50(:,variable,shock)'*plot_uv_scale;
buv84=L_q84(:,variable,shock)'*plot_uv_scale;

[~,~]=jbfill([0:horizon],auv16,buv84,rgb('crimson'),rgb('crimson'),0,0.5);

hold on
%hq16=plot([0:horizon],auv16,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
%hold on
%hq84=plot([0:horizon],buv84,'linestyle','-','linewidth',2,'color',rgb('royalblue'));
hold on
hmed=plot([0:horizon],auv50,'linestyle','-','linewidth',2,'color',rgb('crimson'));
hold on
hzero=plot([0:horizon],auv50*0,'k-.','linewidth',0.5);
hold on
xlabel('Quarters','FontSize',ftsizexlabel,'Interpreter','latex')
set(gca,'FontSize',ftsizeaxis)
set(gca,'LineWidth',ftlinewidth)
title(demand_variables_labels(i),'Interpreter','latex')
box off
legend off
grid on

end

set(gcf, 'PaperPositionMode', 'auto');
print(['epsfiles/',num2str(spec),'_shock_wagebargaining.eps'],'-depsc');
print(['pngfiles/',num2str(spec),'_shock_wagebargaining.png'],'-dpng');


