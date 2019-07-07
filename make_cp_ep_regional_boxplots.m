% Make plots of Ahlstrom regional GPP/NEP response to EP/CP ENSO

% Parameters
wdth = 0.25;
offs = 0.15;
scale = 0.001;
clr = [103,0,31
    178,24,43
    214,96,77
    244,165,130
    253,219,199
    224,224,224
    186,186,186
    135,135,135
    77,77,77
    26,26,26]/255;

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 5];


%% GPP
subplot(2,1,1)
plot([0 9], [0 0], 'k-')
hold on;

% MsTMIP
load ./data/cp_ep_gpp_mstmip_regional.mat;
clear *_tropical_* *_extratropical* *_grass_* *_semiarid_* *_tundra_*;

% Ensemble
b1 = dberg_box(1-offs,CP_GPP_amazon_annual_beta*scale,clr(4,:),'w',wdth,10);
dberg_box(2-offs,CP_GPP_africa_annual_beta*scale,clr(4,:),'w',wdth,10);
dberg_box(4-offs,CP_GPP_austr_annual_beta*scale,clr(4,:),'w',wdth,10);
dberg_box(3-offs,CP_GPP_casia_annual_beta*scale,clr(4,:),'w',wdth,10);
dberg_box(5-offs,CP_GPP_eastus_annual_beta*scale,clr(4,:),'w',wdth,10);
dberg_box(6-offs,CP_GPP_europe_annual_beta*scale,clr(4,:),'w',wdth,10);
dberg_box(7-offs,CP_GPP_sahel_annual_beta*scale,clr(4,:),'w',wdth,10);
dberg_box(8-offs,CP_GPP_westna_annual_beta*scale,clr(4,:),'w',wdth,10);

b2 = dberg_box(1+offs,EP_GPP_amazon_annual_beta*scale,clr(8,:),'w',wdth,10);
dberg_box(2+offs,EP_GPP_africa_annual_beta*scale,clr(8,:),'w',wdth,10);
dberg_box(4+offs,EP_GPP_austr_annual_beta*scale,clr(8,:),'w',wdth,10);
dberg_box(3+offs,EP_GPP_casia_annual_beta*scale,clr(8,:),'w',wdth,10);
dberg_box(5+offs,EP_GPP_eastus_annual_beta*scale,clr(8,:),'w',wdth,10);
dberg_box(6+offs,EP_GPP_europe_annual_beta*scale,clr(8,:),'w',wdth,10);
dberg_box(7+offs,EP_GPP_sahel_annual_beta*scale,clr(8,:),'w',wdth,10);
dberg_box(8+offs,EP_GPP_westna_annual_beta*scale,clr(8,:),'w',wdth,10);

% Means and CIs
plot([1-offs+offs/3 1-offs+offs/3],...
    [CP_GPP_amazon_annual_mean_beta(1)-CP_GPP_amazon_annual_mean_beta_CI(1) CP_GPP_amazon_annual_mean_beta(1)+CP_GPP_amazon_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([2-offs+offs/3 2-offs+offs/3],...
    [CP_GPP_africa_annual_mean_beta(1)-CP_GPP_africa_annual_mean_beta_CI(1) CP_GPP_africa_annual_mean_beta(1)+CP_GPP_africa_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([4-offs+offs/3 4-offs+offs/3],...
    [CP_GPP_austr_annual_mean_beta(1)-CP_GPP_austr_annual_mean_beta_CI(1) CP_GPP_austr_annual_mean_beta(1)+CP_GPP_austr_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([3-offs+offs/3 3-offs+offs/3],...
    [CP_GPP_casia_annual_mean_beta(1)-CP_GPP_casia_annual_mean_beta_CI(1) CP_GPP_casia_annual_mean_beta(1)+CP_GPP_casia_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([5-offs+offs/3 5-offs+offs/3],...
    [CP_GPP_eastus_annual_mean_beta(1)-CP_GPP_eastus_annual_mean_beta_CI(1) CP_GPP_eastus_annual_mean_beta(1)+CP_GPP_eastus_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([6-offs+offs/3 6-offs+offs/3],...
    [CP_GPP_europe_annual_mean_beta(1)-CP_GPP_europe_annual_mean_beta_CI(1) CP_GPP_europe_annual_mean_beta(1)+CP_GPP_europe_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([7-offs+offs/3 7-offs+offs/3],...
    [CP_GPP_sahel_annual_mean_beta(1)-CP_GPP_sahel_annual_mean_beta_CI(1) CP_GPP_sahel_annual_mean_beta(1)+CP_GPP_sahel_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([8-offs+offs/3 8-offs+offs/3],...
    [CP_GPP_westna_annual_mean_beta(1)-CP_GPP_westna_annual_mean_beta_CI(1) CP_GPP_westna_annual_mean_beta(1)+CP_GPP_westna_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));

plot([1+offs+offs/3 1+offs+offs/3],...
    [EP_GPP_amazon_annual_mean_beta(1)-EP_GPP_amazon_annual_mean_beta_CI(1) EP_GPP_amazon_annual_mean_beta(1)+EP_GPP_amazon_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([2+offs+offs/3 2+offs+offs/3],...
    [EP_GPP_africa_annual_mean_beta(1)-EP_GPP_africa_annual_mean_beta_CI(1) EP_GPP_africa_annual_mean_beta(1)+EP_GPP_africa_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([4+offs+offs/3 4+offs+offs/3],...
    [EP_GPP_austr_annual_mean_beta(1)-EP_GPP_austr_annual_mean_beta_CI(1) EP_GPP_austr_annual_mean_beta(1)+EP_GPP_austr_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([3+offs+offs/3 3+offs+offs/3],...
    [EP_GPP_casia_annual_mean_beta(1)-EP_GPP_casia_annual_mean_beta_CI(1) EP_GPP_casia_annual_mean_beta(1)+EP_GPP_casia_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([5+offs+offs/3 5+offs+offs/3],...
    [EP_GPP_eastus_annual_mean_beta(1)-EP_GPP_eastus_annual_mean_beta_CI(1) EP_GPP_eastus_annual_mean_beta(1)+EP_GPP_eastus_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([6+offs+offs/3 6+offs+offs/3],...
    [EP_GPP_europe_annual_mean_beta(1)-EP_GPP_europe_annual_mean_beta_CI(1) EP_GPP_europe_annual_mean_beta(1)+EP_GPP_europe_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([7+offs+offs/3 7+offs+offs/3],...
    [EP_GPP_sahel_annual_mean_beta(1)-EP_GPP_sahel_annual_mean_beta_CI(1) EP_GPP_sahel_annual_mean_beta(1)+EP_GPP_sahel_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([8+offs+offs/3 8+offs+offs/3],...
    [EP_GPP_westna_annual_mean_beta(1)-EP_GPP_westna_annual_mean_beta_CI(1) EP_GPP_westna_annual_mean_beta(1)+EP_GPP_westna_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));

scatter(1-offs+offs/3,CP_GPP_amazon_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(2-offs+offs/3,CP_GPP_africa_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(4-offs+offs/3,CP_GPP_austr_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(3-offs+offs/3,CP_GPP_casia_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(5-offs+offs/3,CP_GPP_eastus_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(6-offs+offs/3,CP_GPP_europe_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(7-offs+offs/3,CP_GPP_sahel_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(8-offs+offs/3,CP_GPP_westna_annual_mean_beta*scale,15,clr(2,:),'filled','^');

scatter(1+offs+offs/3,EP_GPP_amazon_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(2+offs+offs/3,EP_GPP_africa_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(4+offs+offs/3,EP_GPP_austr_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(3+offs+offs/3,EP_GPP_casia_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(5+offs+offs/3,EP_GPP_eastus_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(6+offs+offs/3,EP_GPP_europe_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(7+offs+offs/3,EP_GPP_sahel_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(8+offs+offs/3,EP_GPP_westna_annual_mean_beta*scale,15,clr(10,:),'filled','^');

% LUE
clear CP* EP* ep_* cp_* models;
load ./data/cp_ep_gpp_lue_regional.mat;
clear *_tropical_* *_extratropical* *_grass_* *_semiarid_* *_tundra_*;

plot([1-offs-offs/3 1-offs-offs/3],...
    [CP_GPP_amazon_annual_mean_beta-CP_GPP_amazon_annual_mean_beta_CI CP_GPP_amazon_annual_mean_beta+CP_GPP_amazon_annual_mean_beta_CI]/1000,...
    '-','Color',clr(2,:));
plot([2-offs-offs/3 2-offs-offs/3],...
    [CP_GPP_africa_annual_mean_beta-CP_GPP_africa_annual_mean_beta_CI CP_GPP_africa_annual_mean_beta+CP_GPP_africa_annual_mean_beta_CI]/1000,...
    '-','Color',clr(2,:));
plot([4-offs-offs/3 4-offs-offs/3],...
    [CP_GPP_austr_annual_mean_beta-CP_GPP_austr_annual_mean_beta_CI CP_GPP_austr_annual_mean_beta+CP_GPP_austr_annual_mean_beta_CI]/1000,...
    '-','Color',clr(2,:));
plot([3-offs-offs/3 3-offs-offs/3],...
    [CP_GPP_casia_annual_mean_beta-CP_GPP_casia_annual_mean_beta_CI CP_GPP_casia_annual_mean_beta+CP_GPP_casia_annual_mean_beta_CI]/1000,...
    '-','Color',clr(2,:));
plot([5-offs-offs/3 5-offs-offs/3],...
    [CP_GPP_eastus_annual_mean_beta-CP_GPP_eastus_annual_mean_beta_CI CP_GPP_eastus_annual_mean_beta+CP_GPP_eastus_annual_mean_beta_CI]/1000,...
    '-','Color',clr(2,:));
plot([6-offs-offs/3 6-offs-offs/3],...
    [CP_GPP_europe_annual_mean_beta-CP_GPP_europe_annual_mean_beta_CI CP_GPP_europe_annual_mean_beta+CP_GPP_europe_annual_mean_beta_CI]/1000,...
    '-','Color',clr(2,:));
plot([7-offs-offs/3 7-offs-offs/3],...
    [CP_GPP_sahel_annual_mean_beta-CP_GPP_sahel_annual_mean_beta_CI CP_GPP_sahel_annual_mean_beta+CP_GPP_sahel_annual_mean_beta_CI]/1000,...
    '-','Color',clr(2,:));
plot([8-offs-offs/3 8-offs-offs/3],...
    [CP_GPP_westna_annual_mean_beta-CP_GPP_westna_annual_mean_beta_CI CP_GPP_westna_annual_mean_beta+CP_GPP_westna_annual_mean_beta_CI]/1000,...
    '-','Color',clr(2,:));

plot([1+offs-offs/3 1+offs-offs/3],...
    [EP_GPP_amazon_annual_mean_beta-EP_GPP_amazon_annual_mean_beta_CI EP_GPP_amazon_annual_mean_beta+EP_GPP_amazon_annual_mean_beta_CI]/1000,...
    '-','Color',clr(10,:));
plot([2+offs-offs/3 2+offs-offs/3],...
    [EP_GPP_africa_annual_mean_beta-EP_GPP_africa_annual_mean_beta_CI EP_GPP_africa_annual_mean_beta+EP_GPP_africa_annual_mean_beta_CI]/1000,...
    '-','Color',clr(10,:));
plot([4+offs-offs/3 4+offs-offs/3],...
    [EP_GPP_austr_annual_mean_beta-EP_GPP_austr_annual_mean_beta_CI EP_GPP_austr_annual_mean_beta+EP_GPP_austr_annual_mean_beta_CI]/1000,...
    '-','Color',clr(10,:));
plot([3+offs-offs/3 3+offs-offs/3],...
    [EP_GPP_casia_annual_mean_beta-EP_GPP_casia_annual_mean_beta_CI EP_GPP_casia_annual_mean_beta+EP_GPP_casia_annual_mean_beta_CI]/1000,...
    '-','Color',clr(10,:));
plot([5+offs-offs/3 5+offs-offs/3],...
    [EP_GPP_eastus_annual_mean_beta-EP_GPP_eastus_annual_mean_beta_CI EP_GPP_eastus_annual_mean_beta+EP_GPP_eastus_annual_mean_beta_CI]/1000,...
    '-','Color',clr(10,:));
plot([6+offs-offs/3 6+offs-offs/3],...
    [EP_GPP_europe_annual_mean_beta-EP_GPP_europe_annual_mean_beta_CI EP_GPP_europe_annual_mean_beta+EP_GPP_europe_annual_mean_beta_CI]/1000,...
    '-','Color',clr(10,:));
plot([7+offs-offs/3 7+offs-offs/3],...
    [EP_GPP_sahel_annual_mean_beta-EP_GPP_sahel_annual_mean_beta_CI EP_GPP_sahel_annual_mean_beta+EP_GPP_sahel_annual_mean_beta_CI]/1000,...
    '-','Color',clr(10,:));
plot([8+offs-offs/3 8+offs-offs/3],...
    [EP_GPP_westna_annual_mean_beta-EP_GPP_westna_annual_mean_beta_CI EP_GPP_westna_annual_mean_beta+EP_GPP_westna_annual_mean_beta_CI]/1000,...
    '-','Color',clr(10,:));

scatter(1-offs-offs/3,CP_GPP_amazon_annual_mean_beta*scale,15,clr(2,:),'x');
scatter(2-offs-offs/3,CP_GPP_africa_annual_mean_beta*scale,15,clr(2,:),'x');
scatter(4-offs-offs/3,CP_GPP_austr_annual_mean_beta*scale,15,clr(2,:),'x');
scatter(3-offs-offs/3,CP_GPP_casia_annual_mean_beta*scale,15,clr(2,:),'x');
scatter(5-offs-offs/3,CP_GPP_eastus_annual_mean_beta*scale,15,clr(2,:),'x');
scatter(6-offs-offs/3,CP_GPP_europe_annual_mean_beta*scale,15,clr(2,:),'x');
scatter(7-offs-offs/3,CP_GPP_sahel_annual_mean_beta*scale,15,clr(2,:),'x');
scatter(8-offs-offs/3,CP_GPP_westna_annual_mean_beta*scale,15,clr(2,:),'x');

scatter(1+offs-offs/3,EP_GPP_amazon_annual_mean_beta*scale,15,clr(10,:),'x');
scatter(2+offs-offs/3,EP_GPP_africa_annual_mean_beta*scale,15,clr(10,:),'x');
scatter(4+offs-offs/3,EP_GPP_austr_annual_mean_beta*scale,15,clr(10,:),'x');
scatter(3+offs-offs/3,EP_GPP_casia_annual_mean_beta*scale,15,clr(10,:),'x');
scatter(5+offs-offs/3,EP_GPP_eastus_annual_mean_beta*scale,15,clr(10,:),'x');
scatter(6+offs-offs/3,EP_GPP_europe_annual_mean_beta*scale,15,clr(10,:),'x');
scatter(7+offs-offs/3,EP_GPP_sahel_annual_mean_beta*scale,15,clr(10,:),'x');
scatter(8+offs-offs/3,EP_GPP_westna_annual_mean_beta*scale,15,clr(10,:),'x');

set(gca, 'XLim',[0.5 8.5], 'XTick',1:8, 'TickDir','out',...
    'TickLength',[0.02 0.04], 'XTickLabels',{'Trop. S. Am.','Trop. & S. Africa',...
    'Trop. Asia','Australia','Eastern U.S.', 'Europe', 'The Sahel', 'W. N. America'},...
    'FontSize',7, 'YLim',[-0.5 0.3]); xtickangle(20)
box off;
text(0.6, 0.3, 'A', 'FontSize',12, 'VerticalAlignment','top')

ylabel('GPP response (Pg C yr^{-1} SD^{-1})','FontSize',9)

scatter(6, -0.32, 20, 'k', 'filled','^');
text(6.1, -0.32, 'MsTMIP', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)
scatter(6, -0.4, 20, 'k', 'x');
text(6.1, -0.4, 'LUE', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)

legend([b1 b2], 'CP','EP', 'Location','southeast', 'FontSize',9);
legend('boxoff');


%% NEP
clear *_GPP_*;
subplot(2,1,2)
plot([0 9], [0 0], 'k-')
hold on;

% MsTMIP
load ./data/cp_ep_nep_mstmip_regional.mat;
clear *_tropical_* *_extratropical* *_grass_* *_semiarid_* *_tundra_*;

% Ensemble
dberg_box(1-offs,CP_NEP_amazon_annual_beta*scale,clr(4,:),'w',wdth,15);
dberg_box(2-offs,CP_NEP_africa_annual_beta*scale,clr(4,:),'w',wdth,15);
dberg_box(4-offs,CP_NEP_austr_annual_beta*scale,clr(4,:),'w',wdth,15);
dberg_box(3-offs,CP_NEP_casia_annual_beta*scale,clr(4,:),'w',wdth,15);
dberg_box(5-offs,CP_NEP_eastus_annual_beta*scale,clr(4,:),'w',wdth,15);
dberg_box(6-offs,CP_NEP_europe_annual_beta*scale,clr(4,:),'w',wdth,15);
dberg_box(7-offs,CP_NEP_sahel_annual_beta*scale,clr(4,:),'w',wdth,15);
dberg_box(8-offs,CP_NEP_westna_annual_beta*scale,clr(4,:),'w',wdth,15);

dberg_box(1+offs,EP_NEP_amazon_annual_beta*scale,clr(8,:),'w',wdth,15);
dberg_box(2+offs,EP_NEP_africa_annual_beta*scale,clr(8,:),'w',wdth,15);
dberg_box(4+offs,EP_NEP_austr_annual_beta*scale,clr(8,:),'w',wdth,15);
dberg_box(3+offs,EP_NEP_casia_annual_beta*scale,clr(8,:),'w',wdth,15);
dberg_box(5+offs,EP_NEP_eastus_annual_beta*scale,clr(8,:),'w',wdth,15);
dberg_box(6+offs,EP_NEP_europe_annual_beta*scale,clr(8,:),'w',wdth,15);
dberg_box(7+offs,EP_NEP_sahel_annual_beta*scale,clr(8,:),'w',wdth,15);
dberg_box(8+offs,EP_NEP_westna_annual_beta*scale,clr(8,:),'w',wdth,15);

% Means and CIs
plot([1-offs+offs/3 1-offs+offs/3],...
    [CP_NEP_amazon_annual_mean_beta(1)-CP_NEP_amazon_annual_mean_beta_CI(1) CP_NEP_amazon_annual_mean_beta(1)+CP_NEP_amazon_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([2-offs+offs/3 2-offs+offs/3],...
    [CP_NEP_africa_annual_mean_beta(1)-CP_NEP_africa_annual_mean_beta_CI(1) CP_NEP_africa_annual_mean_beta(1)+CP_NEP_africa_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([4-offs+offs/3 4-offs+offs/3],...
    [CP_NEP_austr_annual_mean_beta(1)-CP_NEP_austr_annual_mean_beta_CI(1) CP_NEP_austr_annual_mean_beta(1)+CP_NEP_austr_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([3-offs+offs/3 3-offs+offs/3],...
    [CP_NEP_casia_annual_mean_beta(1)-CP_NEP_casia_annual_mean_beta_CI(1) CP_NEP_casia_annual_mean_beta(1)+CP_NEP_casia_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([5-offs+offs/3 5-offs+offs/3],...
    [CP_NEP_eastus_annual_mean_beta(1)-CP_NEP_eastus_annual_mean_beta_CI(1) CP_NEP_eastus_annual_mean_beta(1)+CP_NEP_eastus_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([6-offs+offs/3 6-offs+offs/3],...
    [CP_NEP_europe_annual_mean_beta(1)-CP_NEP_europe_annual_mean_beta_CI(1) CP_NEP_europe_annual_mean_beta(1)+CP_NEP_europe_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([7-offs+offs/3 7-offs+offs/3],...
    [CP_NEP_sahel_annual_mean_beta(1)-CP_NEP_sahel_annual_mean_beta_CI(1) CP_NEP_sahel_annual_mean_beta(1)+CP_NEP_sahel_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([8-offs+offs/3 8-offs+offs/3],...
    [CP_NEP_westna_annual_mean_beta(1)-CP_NEP_westna_annual_mean_beta_CI(1) CP_NEP_westna_annual_mean_beta(1)+CP_NEP_westna_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));

plot([1+offs+offs/3 1+offs+offs/3],...
    [EP_NEP_amazon_annual_mean_beta(1)-EP_NEP_amazon_annual_mean_beta_CI(1) EP_NEP_amazon_annual_mean_beta(1)+EP_NEP_amazon_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([2+offs+offs/3 2+offs+offs/3],...
    [EP_NEP_africa_annual_mean_beta(1)-EP_NEP_africa_annual_mean_beta_CI(1) EP_NEP_africa_annual_mean_beta(1)+EP_NEP_africa_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([4+offs+offs/3 4+offs+offs/3],...
    [EP_NEP_austr_annual_mean_beta(1)-EP_NEP_austr_annual_mean_beta_CI(1) EP_NEP_austr_annual_mean_beta(1)+EP_NEP_austr_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([3+offs+offs/3 3+offs+offs/3],...
    [EP_NEP_casia_annual_mean_beta(1)-EP_NEP_casia_annual_mean_beta_CI(1) EP_NEP_casia_annual_mean_beta(1)+EP_NEP_casia_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([5+offs+offs/3 5+offs+offs/3],...
    [EP_NEP_eastus_annual_mean_beta(1)-EP_NEP_eastus_annual_mean_beta_CI(1) EP_NEP_eastus_annual_mean_beta(1)+EP_NEP_eastus_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([6+offs+offs/3 6+offs+offs/3],...
    [EP_NEP_europe_annual_mean_beta(1)-EP_NEP_europe_annual_mean_beta_CI(1) EP_NEP_europe_annual_mean_beta(1)+EP_NEP_europe_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([7+offs+offs/3 7+offs+offs/3],...
    [EP_NEP_sahel_annual_mean_beta(1)-EP_NEP_sahel_annual_mean_beta_CI(1) EP_NEP_sahel_annual_mean_beta(1)+EP_NEP_sahel_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([8+offs+offs/3 8+offs+offs/3],...
    [EP_NEP_westna_annual_mean_beta(1)-EP_NEP_westna_annual_mean_beta_CI(1) EP_NEP_westna_annual_mean_beta(1)+EP_NEP_westna_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));

scatter(1-offs+offs/3,CP_NEP_amazon_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(2-offs+offs/3,CP_NEP_africa_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(4-offs+offs/3,CP_NEP_austr_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(3-offs+offs/3,CP_NEP_casia_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(5-offs+offs/3,CP_NEP_eastus_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(6-offs+offs/3,CP_NEP_europe_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(7-offs+offs/3,CP_NEP_sahel_annual_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(8-offs+offs/3,CP_NEP_westna_annual_mean_beta*scale,15,clr(2,:),'filled','^');

scatter(1+offs+offs/3,EP_NEP_amazon_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(2+offs+offs/3,EP_NEP_africa_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(4+offs+offs/3,EP_NEP_austr_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(3+offs+offs/3,EP_NEP_casia_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(5+offs+offs/3,EP_NEP_eastus_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(6+offs+offs/3,EP_NEP_europe_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(7+offs+offs/3,EP_NEP_sahel_annual_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(8+offs+offs/3,EP_NEP_westna_annual_mean_beta*scale,15,clr(10,:),'filled','^');

% Inversions
clear CP* EP* ep_* cp_* models;
load ./data/cp_ep_nep_inversions_regional.mat;
clear *_tropical_* *_extratropical* *_grass_* *_semiarid_* *_tundra_*;

% Means and CIs
plot([1-offs-offs/3 1-offs-offs/3],...
    [CP_NEP_amazon_annual_mean_beta(1)-CP_NEP_amazon_annual_mean_beta_CI(1) CP_NEP_amazon_annual_mean_beta(1)+CP_NEP_amazon_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([2-offs-offs/3 2-offs-offs/3],...
    [CP_NEP_africa_annual_mean_beta(1)-CP_NEP_africa_annual_mean_beta_CI(1) CP_NEP_africa_annual_mean_beta(1)+CP_NEP_africa_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([4-offs-offs/3 4-offs-offs/3],...
    [CP_NEP_austr_annual_mean_beta(1)-CP_NEP_austr_annual_mean_beta_CI(1) CP_NEP_austr_annual_mean_beta(1)+CP_NEP_austr_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([3-offs-offs/3 3-offs-offs/3],...
    [CP_NEP_casia_annual_mean_beta(1)-CP_NEP_casia_annual_mean_beta_CI(1) CP_NEP_casia_annual_mean_beta(1)+CP_NEP_casia_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([5-offs-offs/3 5-offs-offs/3],...
    [CP_NEP_eastus_annual_mean_beta(1)-CP_NEP_eastus_annual_mean_beta_CI(1) CP_NEP_eastus_annual_mean_beta(1)+CP_NEP_eastus_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([6-offs-offs/3 6-offs-offs/3],...
    [CP_NEP_europe_annual_mean_beta(1)-CP_NEP_europe_annual_mean_beta_CI(1) CP_NEP_europe_annual_mean_beta(1)+CP_NEP_europe_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([7-offs-offs/3 7-offs-offs/3],...
    [CP_NEP_sahel_annual_mean_beta(1)-CP_NEP_sahel_annual_mean_beta_CI(1) CP_NEP_sahel_annual_mean_beta(1)+CP_NEP_sahel_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([8-offs-offs/3 8-offs-offs/3],...
    [CP_NEP_westna_annual_mean_beta(1)-CP_NEP_westna_annual_mean_beta_CI(1) CP_NEP_westna_annual_mean_beta(1)+CP_NEP_westna_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));

plot([1+offs-offs/3 1+offs-offs/3],...
    [EP_NEP_amazon_annual_mean_beta(1)-EP_NEP_amazon_annual_mean_beta_CI(1) EP_NEP_amazon_annual_mean_beta(1)+EP_NEP_amazon_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([2+offs-offs/3 2+offs-offs/3],...
    [EP_NEP_africa_annual_mean_beta(1)-EP_NEP_africa_annual_mean_beta_CI(1) EP_NEP_africa_annual_mean_beta(1)+EP_NEP_africa_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([4+offs-offs/3 4+offs-offs/3],...
    [EP_NEP_austr_annual_mean_beta(1)-EP_NEP_austr_annual_mean_beta_CI(1) EP_NEP_austr_annual_mean_beta(1)+EP_NEP_austr_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([3+offs-offs/3 3+offs-offs/3],...
    [EP_NEP_casia_annual_mean_beta(1)-EP_NEP_casia_annual_mean_beta_CI(1) EP_NEP_casia_annual_mean_beta(1)+EP_NEP_casia_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([5+offs-offs/3 5+offs-offs/3],...
    [EP_NEP_eastus_annual_mean_beta(1)-EP_NEP_eastus_annual_mean_beta_CI(1) EP_NEP_eastus_annual_mean_beta(1)+EP_NEP_eastus_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([6+offs-offs/3 6+offs-offs/3],...
    [EP_NEP_europe_annual_mean_beta(1)-EP_NEP_europe_annual_mean_beta_CI(1) EP_NEP_europe_annual_mean_beta(1)+EP_NEP_europe_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([7+offs-offs/3 7+offs-offs/3],...
    [EP_NEP_sahel_annual_mean_beta(1)-EP_NEP_sahel_annual_mean_beta_CI(1) EP_NEP_sahel_annual_mean_beta(1)+EP_NEP_sahel_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([8+offs-offs/3 8+offs-offs/3],...
    [EP_NEP_westna_annual_mean_beta(1)-EP_NEP_westna_annual_mean_beta_CI(1) EP_NEP_westna_annual_mean_beta(1)+EP_NEP_westna_annual_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));

scatter(1-offs-offs/3,CP_NEP_amazon_annual_mean_beta*scale,15,clr(2,:),'filled','s');
scatter(2-offs-offs/3,CP_NEP_africa_annual_mean_beta*scale,15,clr(2,:),'filled','s');
scatter(4-offs-offs/3,CP_NEP_austr_annual_mean_beta*scale,15,clr(2,:),'filled','s');
scatter(3-offs-offs/3,CP_NEP_casia_annual_mean_beta*scale,15,clr(2,:),'filled','s');
scatter(5-offs-offs/3,CP_NEP_eastus_annual_mean_beta*scale,15,clr(2,:),'filled','s');
scatter(6-offs-offs/3,CP_NEP_europe_annual_mean_beta*scale,15,clr(2,:),'filled','s');
scatter(7-offs-offs/3,CP_NEP_sahel_annual_mean_beta*scale,15,clr(2,:),'filled','s');
scatter(8-offs-offs/3,CP_NEP_westna_annual_mean_beta*scale,15,clr(2,:),'filled','s');

scatter(1+offs-offs/3,EP_NEP_amazon_annual_mean_beta*scale,15,clr(10,:),'filled','s');
scatter(2+offs-offs/3,EP_NEP_africa_annual_mean_beta*scale,15,clr(10,:),'filled','s');
scatter(4+offs-offs/3,EP_NEP_austr_annual_mean_beta*scale,15,clr(10,:),'filled','s');
scatter(3+offs-offs/3,EP_NEP_casia_annual_mean_beta*scale,15,clr(10,:),'filled','s');
scatter(5+offs-offs/3,EP_NEP_eastus_annual_mean_beta*scale,15,clr(10,:),'filled','s');
scatter(6+offs-offs/3,EP_NEP_europe_annual_mean_beta*scale,15,clr(10,:),'filled','s');
scatter(7+offs-offs/3,EP_NEP_sahel_annual_mean_beta*scale,15,clr(10,:),'filled','s');
scatter(8+offs-offs/3,EP_NEP_westna_annual_mean_beta*scale,15,clr(10,:),'filled','s');


set(gca, 'XLim',[0.5 8.5], 'XTick',1:8, 'TickDir','out',...
    'TickLength',[0.02 0.04], 'XTickLabels',{'Trop. S. Am.','Trop. & S. Africa',...
    'Trop. Asia','Australia','Eastern U.S.', 'Europe', 'The Sahel', 'W. N. America'},...
    'FontSize',7, 'YLim',[-0.33 0.23]); xtickangle(20);
box off;
text(0.6, 0.23, 'B', 'FontSize',12, 'VerticalAlignment','top')

scatter(6, -0.22, 20, 'k', 'filled','s');
text(6.1, -0.22, 'CO_{2} Inversions', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)
scatter(6, -0.28, 20, 'k', 'filled','^');
text(6.1, -0.28, 'MsTMIP', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)

ylabel('NEP response (Pg C yr^{-1} SD^{-1})','FontSize',9)


set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/epi-cpi-gpp-nee-bySpecificRegion.tif')
close all;

