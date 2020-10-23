% Make plots of Ahlstrom regional GPP/NEE response to EP/CP ENSO

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
h.Position = [1 1 5 5];


%% GPP
subplot(2,1,1)
plot([0 6], [0 0], 'k-')
hold on;

% MsTMIP
load ./data/cp_ep_gpp_mstmip_regional.mat;

% Ensemble
b1 = dberg_box(1-offs,CP_GPP_tropical_shyear_beta*scale,clr(4,:),'w',wdth,10);
dberg_box(2-offs,CP_GPP_semiarid_shyear_beta*scale,clr(4,:),'w',wdth,10);
dberg_box(4-offs,CP_GPP_grass_shyear_beta*scale,clr(4,:),'w',wdth,10);
dberg_box(3-offs,CP_GPP_extratropical_shyear_beta*scale,clr(4,:),'w',wdth,10);
dberg_box(5-offs,CP_GPP_tundra_shyear_beta*scale,clr(4,:),'w',wdth,10);

b2 = dberg_box(1+offs,EP_GPP_tropical_shyear_beta*scale,clr(8,:),'w',wdth,10);
dberg_box(2+offs,EP_GPP_semiarid_shyear_beta*scale,clr(8,:),'w',wdth,10);
dberg_box(4+offs,EP_GPP_grass_shyear_beta*scale,clr(8,:),'w',wdth,10);
dberg_box(3+offs,EP_GPP_extratropical_shyear_beta*scale,clr(8,:),'w',wdth,10);
dberg_box(5+offs,EP_GPP_tundra_shyear_beta*scale,clr(8,:),'w',wdth,10);

% Means and CIs
plot([1-offs+offs/3 1-offs+offs/3],...
    [CP_GPP_tropical_shyear_mean_beta(1)-CP_GPP_tropical_shyear_mean_beta_CI(1) CP_GPP_tropical_shyear_mean_beta(1)+CP_GPP_tropical_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([2-offs+offs/3 2-offs+offs/3],...
    [CP_GPP_semiarid_shyear_mean_beta(1)-CP_GPP_semiarid_shyear_mean_beta_CI(1) CP_GPP_semiarid_shyear_mean_beta(1)+CP_GPP_semiarid_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([4-offs+offs/3 4-offs+offs/3],...
    [CP_GPP_grass_shyear_mean_beta(1)-CP_GPP_grass_shyear_mean_beta_CI(1) CP_GPP_grass_shyear_mean_beta(1)+CP_GPP_grass_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([3-offs+offs/3 3-offs+offs/3],...
    [CP_GPP_extratropical_shyear_mean_beta(1)-CP_GPP_extratropical_shyear_mean_beta_CI(1) CP_GPP_extratropical_shyear_mean_beta(1)+CP_GPP_extratropical_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([5-offs+offs/3 5-offs+offs/3],...
    [CP_GPP_tundra_shyear_mean_beta(1)-CP_GPP_tundra_shyear_mean_beta_CI(1) CP_GPP_tundra_shyear_mean_beta(1)+CP_GPP_tundra_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));

plot([1+offs+offs/3 1+offs+offs/3],...
    [EP_GPP_tropical_shyear_mean_beta(1)-EP_GPP_tropical_shyear_mean_beta_CI(1) EP_GPP_tropical_shyear_mean_beta(1)+EP_GPP_tropical_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([2+offs+offs/3 2+offs+offs/3],...
    [EP_GPP_semiarid_shyear_mean_beta(1)-EP_GPP_semiarid_shyear_mean_beta_CI(1) EP_GPP_semiarid_shyear_mean_beta(1)+EP_GPP_semiarid_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([4+offs+offs/3 4+offs+offs/3],...
    [EP_GPP_grass_shyear_mean_beta(1)-EP_GPP_grass_shyear_mean_beta_CI(1) EP_GPP_grass_shyear_mean_beta(1)+EP_GPP_grass_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([3+offs+offs/3 3+offs+offs/3],...
    [EP_GPP_extratropical_shyear_mean_beta(1)-EP_GPP_extratropical_shyear_mean_beta_CI(1) EP_GPP_extratropical_shyear_mean_beta(1)+EP_GPP_extratropical_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([5+offs+offs/3 5+offs+offs/3],...
    [EP_GPP_tundra_shyear_mean_beta(1)-EP_GPP_tundra_shyear_mean_beta_CI(1) EP_GPP_tundra_shyear_mean_beta(1)+EP_GPP_tundra_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));

scatter(1-offs+offs/3,CP_GPP_tropical_shyear_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(2-offs+offs/3,CP_GPP_semiarid_shyear_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(4-offs+offs/3,CP_GPP_grass_shyear_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(3-offs+offs/3,CP_GPP_extratropical_shyear_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(5-offs+offs/3,CP_GPP_tundra_shyear_mean_beta*scale,15,clr(2,:),'filled','^');

scatter(1+offs+offs/3,EP_GPP_tropical_shyear_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(2+offs+offs/3,EP_GPP_semiarid_shyear_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(4+offs+offs/3,EP_GPP_grass_shyear_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(3+offs+offs/3,EP_GPP_extratropical_shyear_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(5+offs+offs/3,EP_GPP_tundra_shyear_mean_beta*scale,15,clr(10,:),'filled','^');

% LUE
clear CP* EP* ep_* cp_* models;
load ./data/cp_ep_gpp_lue_regional.mat;

plot([1-offs-offs/3 1-offs-offs/3],...
    [CP_GPP_tropical_shyear_mean_beta-CP_GPP_tropical_shyear_mean_beta_CI CP_GPP_tropical_shyear_mean_beta+CP_GPP_tropical_shyear_mean_beta_CI]/1000,...
    '-','Color',clr(2,:));
plot([2-offs-offs/3 2-offs-offs/3],...
    [CP_GPP_semiarid_shyear_mean_beta-CP_GPP_semiarid_shyear_mean_beta_CI CP_GPP_semiarid_shyear_mean_beta+CP_GPP_semiarid_shyear_mean_beta_CI]/1000,...
    '-','Color',clr(2,:));
plot([4-offs-offs/3 4-offs-offs/3],...
    [CP_GPP_grass_shyear_mean_beta-CP_GPP_grass_shyear_mean_beta_CI CP_GPP_grass_shyear_mean_beta+CP_GPP_grass_shyear_mean_beta_CI]/1000,...
    '-','Color',clr(2,:));
plot([3-offs-offs/3 3-offs-offs/3],...
    [CP_GPP_extratropical_shyear_mean_beta-CP_GPP_extratropical_shyear_mean_beta_CI CP_GPP_extratropical_shyear_mean_beta+CP_GPP_extratropical_shyear_mean_beta_CI]/1000,...
    '-','Color',clr(2,:));
plot([5-offs-offs/3 5-offs-offs/3],...
    [CP_GPP_tundra_shyear_mean_beta-CP_GPP_tundra_shyear_mean_beta_CI CP_GPP_tundra_shyear_mean_beta+CP_GPP_tundra_shyear_mean_beta_CI]/1000,...
    '-','Color',clr(2,:));

plot([1+offs-offs/3 1+offs-offs/3],...
    [EP_GPP_tropical_shyear_mean_beta-EP_GPP_tropical_shyear_mean_beta_CI EP_GPP_tropical_shyear_mean_beta+EP_GPP_tropical_shyear_mean_beta_CI]/1000,...
    '-','Color',clr(10,:));
plot([2+offs-offs/3 2+offs-offs/3],...
    [EP_GPP_semiarid_shyear_mean_beta-EP_GPP_semiarid_shyear_mean_beta_CI EP_GPP_semiarid_shyear_mean_beta+EP_GPP_semiarid_shyear_mean_beta_CI]/1000,...
    '-','Color',clr(10,:));
plot([4+offs-offs/3 4+offs-offs/3],...
    [EP_GPP_grass_shyear_mean_beta-EP_GPP_grass_shyear_mean_beta_CI EP_GPP_grass_shyear_mean_beta+EP_GPP_grass_shyear_mean_beta_CI]/1000,...
    '-','Color',clr(10,:));
plot([3+offs-offs/3 3+offs-offs/3],...
    [EP_GPP_extratropical_shyear_mean_beta-EP_GPP_extratropical_shyear_mean_beta_CI EP_GPP_extratropical_shyear_mean_beta+EP_GPP_extratropical_shyear_mean_beta_CI]/1000,...
    '-','Color',clr(10,:));
plot([5+offs-offs/3 5+offs-offs/3],...
    [EP_GPP_tundra_shyear_mean_beta-EP_GPP_tundra_shyear_mean_beta_CI EP_GPP_tundra_shyear_mean_beta+EP_GPP_tundra_shyear_mean_beta_CI]/1000,...
    '-','Color',clr(10,:));

scatter(1-offs-offs/3,CP_GPP_tropical_shyear_mean_beta*scale,15,clr(2,:),'x');
scatter(2-offs-offs/3,CP_GPP_semiarid_shyear_mean_beta*scale,15,clr(2,:),'x');
scatter(4-offs-offs/3,CP_GPP_grass_shyear_mean_beta*scale,15,clr(2,:),'x');
scatter(3-offs-offs/3,CP_GPP_extratropical_shyear_mean_beta*scale,15,clr(2,:),'x');
scatter(5-offs-offs/3,CP_GPP_tundra_shyear_mean_beta*scale,15,clr(2,:),'x');

scatter(1+offs-offs/3,EP_GPP_tropical_shyear_mean_beta*scale,15,clr(10,:),'x');
scatter(2+offs-offs/3,EP_GPP_semiarid_shyear_mean_beta*scale,15,clr(10,:),'x');
scatter(4+offs-offs/3,EP_GPP_grass_shyear_mean_beta*scale,15,clr(10,:),'x');
scatter(3+offs-offs/3,EP_GPP_extratropical_shyear_mean_beta*scale,15,clr(10,:),'x');
scatter(5+offs-offs/3,EP_GPP_tundra_shyear_mean_beta*scale,15,clr(10,:),'x');

set(gca, 'XLim',[0.5 5.5], 'XTick',1:5, 'TickDir','out',...
    'TickLength',[0.02 0.04], 'XTickLabels',{'Tropical Forest','Semiarid',...
    'Extratropical Forest','Grass/Crop','Tundra/Arctic Shrub'},...
    'FontSize',7, 'YLim',[-0.7 0.4], 'YTick',-0.8:0.2:0.4)
ylim = get(gca, 'YLim');
box off;
text(0.6, ylim(2), 'A', 'FontSize',12, 'VerticalAlignment','top')

ylabel('GPP response (Pg C yr^{-1} K^{-1})','FontSize',9)

scatter(3.5, -0.46, 20, 'k', 'filled','^');
text(3.6, -0.46, 'MsTMIP', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)
scatter(3.5, -0.56, 20, 'k', 'x');
text(3.6, -0.56, 'LUE', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)

legend([b1 b2], 'CP','EP', 'Location','southeast', 'FontSize',9);
legend('boxoff');


%% NEP
clear *_GPP_*;
subplot(2,1,2)
plot([0 6], [0 0], 'k-')
hold on;

% MsTMIP
load ./data/cp_ep_nep_mstmip_regional.mat;

% Ensemble
dberg_box(1-offs,CP_NEP_tropical_shyear_beta*scale,clr(4,:),'w',wdth,15);
dberg_box(2-offs,CP_NEP_semiarid_shyear_beta*scale,clr(4,:),'w',wdth,15);
dberg_box(4-offs,CP_NEP_grass_shyear_beta*scale,clr(4,:),'w',wdth,15);
dberg_box(3-offs,CP_NEP_extratropical_shyear_beta*scale,clr(4,:),'w',wdth,15);
dberg_box(5-offs,CP_NEP_tundra_shyear_beta*scale,clr(4,:),'w',wdth,15);

dberg_box(1+offs,EP_NEP_tropical_shyear_beta*scale,clr(8,:),'w',wdth,15);
dberg_box(2+offs,EP_NEP_semiarid_shyear_beta*scale,clr(8,:),'w',wdth,15);
dberg_box(4+offs,EP_NEP_grass_shyear_beta*scale,clr(8,:),'w',wdth,15);
dberg_box(3+offs,EP_NEP_extratropical_shyear_beta*scale,clr(8,:),'w',wdth,15);
dberg_box(5+offs,EP_NEP_tundra_shyear_beta*scale,clr(8,:),'w',wdth,15);

% Means and CIs
plot([1-offs 1-offs],...
    [CP_NEP_tropical_shyear_mean_beta(1)-CP_NEP_tropical_shyear_mean_beta_CI(1) CP_NEP_tropical_shyear_mean_beta(1)+CP_NEP_tropical_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([2-offs 2-offs],...
    [CP_NEP_semiarid_shyear_mean_beta(1)-CP_NEP_semiarid_shyear_mean_beta_CI(1) CP_NEP_semiarid_shyear_mean_beta(1)+CP_NEP_semiarid_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([4-offs 4-offs],...
    [CP_NEP_grass_shyear_mean_beta(1)-CP_NEP_grass_shyear_mean_beta_CI(1) CP_NEP_grass_shyear_mean_beta(1)+CP_NEP_grass_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([3-offs 3-offs],...
    [CP_NEP_extratropical_shyear_mean_beta(1)-CP_NEP_extratropical_shyear_mean_beta_CI(1) CP_NEP_extratropical_shyear_mean_beta(1)+CP_NEP_extratropical_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));
plot([5-offs 5-offs],...
    [CP_NEP_tundra_shyear_mean_beta(1)-CP_NEP_tundra_shyear_mean_beta_CI(1) CP_NEP_tundra_shyear_mean_beta(1)+CP_NEP_tundra_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(2,:));

plot([1+offs 1+offs],...
    [EP_NEP_tropical_shyear_mean_beta(1)-EP_NEP_tropical_shyear_mean_beta_CI(1) EP_NEP_tropical_shyear_mean_beta(1)+EP_NEP_tropical_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([2+offs 2+offs],...
    [EP_NEP_semiarid_shyear_mean_beta(1)-EP_NEP_semiarid_shyear_mean_beta_CI(1) EP_NEP_semiarid_shyear_mean_beta(1)+EP_NEP_semiarid_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([4+offs 4+offs],...
    [EP_NEP_grass_shyear_mean_beta(1)-EP_NEP_grass_shyear_mean_beta_CI(1) EP_NEP_grass_shyear_mean_beta(1)+EP_NEP_grass_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([3+offs 3+offs],...
    [EP_NEP_extratropical_shyear_mean_beta(1)-EP_NEP_extratropical_shyear_mean_beta_CI(1) EP_NEP_extratropical_shyear_mean_beta(1)+EP_NEP_extratropical_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));
plot([5+offs 5+offs],...
    [EP_NEP_tundra_shyear_mean_beta(1)-EP_NEP_tundra_shyear_mean_beta_CI(1) EP_NEP_tundra_shyear_mean_beta(1)+EP_NEP_tundra_shyear_mean_beta_CI(1)]/1000,...
    '-','Color',clr(10,:));

scatter(1-offs,CP_NEP_tropical_shyear_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(2-offs,CP_NEP_semiarid_shyear_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(4-offs,CP_NEP_grass_shyear_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(3-offs,CP_NEP_extratropical_shyear_mean_beta*scale,15,clr(2,:),'filled','^');
scatter(5-offs,CP_NEP_tundra_shyear_mean_beta*scale,15,clr(2,:),'filled','^');

scatter(1+offs,EP_NEP_tropical_shyear_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(2+offs,EP_NEP_semiarid_shyear_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(4+offs,EP_NEP_grass_shyear_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(3+offs,EP_NEP_extratropical_shyear_mean_beta*scale,15,clr(10,:),'filled','^');
scatter(5+offs,EP_NEP_tundra_shyear_mean_beta*scale,15,clr(10,:),'filled','^');


set(gca, 'XLim',[0.5 5.5], 'XTick',1:5, 'TickDir','out',...
    'TickLength',[0.02 0.04], 'XTickLabels',{'Tropical Forest','Semiarid',...
    'Extratropical Forest','Grass/Crop','Tundra/Arctic Shrub'},...
    'FontSize',7, 'YLim',[-0.6 0.2], 'YTick',-0.6:0.2:0.2)
ylim = get(gca, 'YLim');
box off;
text(0.6, ylim(2), 'B', 'FontSize',12, 'VerticalAlignment','top')

ylabel('NEP response (Pg C yr^{-1} K^{-1})','FontSize',9)


set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-epi-cpi-gpp-nee-byAhlstromRegion-shyear.tif')
close all;

