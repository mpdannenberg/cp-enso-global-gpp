% Make some maps and stuff

clr = cbrewer('div','RdBu',12);

latlim = [-80 80];
lonlim = [-180 180];
worldland = shaperead('landareas','UseGeoCoords', true);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 6];

set(h, 'defaultAxesColorOrder',[0 0 0; 0 0 0]/255)

%% Map MsTMIP Beta
load('./data/cp_ep_nep_mstmip.mat');

subplot(3,2,1)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, EP_NEP_annual_beta)
caxis([-0.03 0.03])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
text(-2.2,1.3,'A', 'FontSize',12);
ttl = title('Eastern Pacific', 'FontSize',12);
ttl.Position(2) = 1.75;

subplot(3,2,2)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, CP_NEP_annual_beta)
caxis([-0.03 0.03])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
set(gca, 'Position',pos);
text(-2.2,1.3,'B', 'FontSize',12);
ttl = title('Central Pacific', 'FontSize',12);
ttl.Position(2) = 1.75;

cb = colorbar('eastoutside');
cb.Position = [0.52 0.71 0.04 0.21];
cb.Ticks = -0.03:0.005:0.03;
cb.TickLength = 0.23;
cb.TickLabels = {'-30','','-20','','-10','','0','','10','','20','','30'};
% cb.YDir = 'reverse';
ylb = ylabel(cb, {'Mean NEP response'; '(g C m^{-2} yr^{-1} K^{-1})'});
ylb.Position = [-0.75 0.0000 0];

%% Plot MsTMIP beta through time
wdth = 0.32;
offs = 0.19;
scale = 0.001;
clr2 = [103,0,31
    178,24,43
    214,96,77
    244,165,130
    253,219,199
    224,224,224
    186,186,186
    135,135,135
    77,77,77
    26,26,26]/255;

ax1 = subplot(3,2,[3 4]);
plot([0 22], [0 0], 'k-')
set(ax1, 'Position', [0.09 0.41 0.835 0.22])

% Monthly
yyaxis left;
hold on;
for i = 1:18
    
    dberg_box(i-offs,CP_NEP_global_monthly_beta(i,:), clr2(4,:), 'w', wdth, 10);
    dberg_box(i+offs,EP_NEP_global_monthly_beta(i,:), clr2(8,:), 'w', wdth, 10);
    
    plot([i-offs+offs/3 i-offs+offs/3],...
        [CP_NEP_global_monthly_mean_beta(i)-CP_NEP_global_monthly_mean_beta_CI(i) CP_NEP_global_monthly_mean_beta(i)+CP_NEP_global_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(2,:));
    scatter(i-offs+offs/3,CP_NEP_global_monthly_mean_beta(i),10,clr2(2,:),'filled','^');
    
    plot([i+offs+offs/3 i+offs+offs/3],...
        [EP_NEP_global_monthly_mean_beta(i)-EP_NEP_global_monthly_mean_beta_CI(i) EP_NEP_global_monthly_mean_beta(i)+EP_NEP_global_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(10,:));
    scatter(i+offs+offs/3,EP_NEP_global_monthly_mean_beta(i),10,clr2(10,:),'filled','^');
    
    
end
set(gca, 'XLim',[0.5 21.5], 'YLim',[-6 4], 'XTick',[1:18 20 21], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.01 0.05], 'XTickLabels','');
ylb = ylabel('Monthly NEP response (Tg C day^{-1} K^{-1})', 'FontSize',7);
text(0.9, 3.5, 'C) Global', 'FontSize',12);

% Prior Jul-Jun
yyaxis right;
b1 = dberg_box(20-offs,CP_NEP_global_shyear_beta*scale, clr2(4,:), 'w', wdth, 10);
b2 = dberg_box(20+offs,EP_NEP_global_shyear_beta*scale, clr2(8,:), 'w', wdth, 10);
plot([20-offs+offs/3 20-offs+offs/3],...
    [CP_NEP_global_shyear_mean_beta*scale-CP_NEP_global_shyear_mean_beta_CI*scale CP_NEP_global_shyear_mean_beta*scale+CP_NEP_global_shyear_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(20-offs+offs/3,CP_NEP_global_shyear_mean_beta*scale,10,clr2(2,:),'filled','^');
plot([20+offs+offs/3 20+offs+offs/3],...
    [EP_NEP_global_shyear_mean_beta*scale-EP_NEP_global_shyear_mean_beta_CI*scale EP_NEP_global_shyear_mean_beta*scale+EP_NEP_global_shyear_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(20+offs+offs/3,EP_NEP_global_shyear_mean_beta*scale,10,clr2(10,:),'filled','^');

% Calendar year
b1 = dberg_box(21-offs,CP_NEP_global_annual_beta*scale, clr2(4,:), 'w', wdth, 10);
b2 = dberg_box(21+offs,EP_NEP_global_annual_beta*scale, clr2(8,:), 'w', wdth, 10);
plot([21-offs+offs/3 21-offs+offs/3],...
    [CP_NEP_global_annual_mean_beta*scale-CP_NEP_global_annual_mean_beta_CI*scale CP_NEP_global_annual_mean_beta*scale+CP_NEP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(21-offs+offs/3,CP_NEP_global_annual_mean_beta*scale,10,clr2(2,:),'filled','^');
plot([21+offs+offs/3 21+offs+offs/3],...
    [EP_NEP_global_annual_mean_beta*scale-EP_NEP_global_annual_mean_beta_CI*scale EP_NEP_global_annual_mean_beta*scale+EP_NEP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(21+offs+offs/3,EP_NEP_global_annual_mean_beta*scale,10,clr2(10,:),'filled','^');

set(gca, 'YLim',[-1500 1000]/1000);
ylb = ylabel('Annual NEP response (Pg C yr^{-1} K^{-1})', 'FontSize',7);

hold off;
box off;


%% Plot MsTMIP beta through time
wdth = 0.32;
offs = 0.19;
scale = 0.001;
clr2 = [103,0,31
    178,24,43
    214,96,77
    244,165,130
    253,219,199
    224,224,224
    186,186,186
    135,135,135
    77,77,77
    26,26,26]/255;

ax2 = subplot(3,2,[5 6]);
plot([0 22], [0 0], 'k-')
set(ax2, 'Position', [0.09 0.08 0.835 0.22])

% Monthly
yyaxis left;
hold on;
for i = 1:18
    
    dberg_box(i-offs,CP_NEP_tropics_monthly_beta(i,:), clr2(4,:), 'w', wdth, 10);
    dberg_box(i+offs,EP_NEP_tropics_monthly_beta(i,:), clr2(8,:), 'w', wdth, 10);
    
    plot([i-offs+offs/3 i-offs+offs/3],...
        [CP_NEP_tropics_monthly_mean_beta(i)-CP_NEP_tropics_monthly_mean_beta_CI(i) CP_NEP_tropics_monthly_mean_beta(i)+CP_NEP_tropics_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(2,:));
    scatter(i-offs+offs/3,CP_NEP_tropics_monthly_mean_beta(i),10,clr2(2,:),'filled','^');
    
    plot([i+offs+offs/3 i+offs+offs/3],...
        [EP_NEP_tropics_monthly_mean_beta(i)-EP_NEP_tropics_monthly_mean_beta_CI(i) EP_NEP_tropics_monthly_mean_beta(i)+EP_NEP_tropics_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(10,:));
    scatter(i+offs+offs/3,EP_NEP_tropics_monthly_mean_beta(i),10,clr2(10,:),'filled','^');
    
    
end
set(gca, 'XLim',[0.5 21.5], 'YLim',[-6 4], 'XTick',[1:18 20 21], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.01 0.05], 'XTickLabels',...
    {'Jul*','Aug*','Sep*','Oct*','Nov*','Dec*','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec','Jul*-Jun','Jan-Dec'});
xtickangle(-45)
ylb = ylabel('Monthly NEP response (Tg C day^{-1} K^{-1})', 'FontSize',7);
text(0.9, 3.5, 'D) Tropics', 'FontSize',12);

% prior Jul - Jun
yyaxis right;
b1 = dberg_box(20-offs,CP_NEP_tropics_shyear_beta*scale, clr2(4,:), 'w', wdth, 10);
b2 = dberg_box(20+offs,EP_NEP_tropics_shyear_beta*scale, clr2(8,:), 'w', wdth, 10);
plot([20-offs+offs/3 20-offs+offs/3],...
    [CP_NEP_tropics_shyear_mean_beta*scale-CP_NEP_tropics_shyear_mean_beta_CI*scale CP_NEP_tropics_shyear_mean_beta*scale+CP_NEP_tropics_shyear_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(20-offs+offs/3,CP_NEP_tropics_shyear_mean_beta*scale,10,clr2(2,:),'filled','^');
plot([20+offs+offs/3 20+offs+offs/3],...
    [EP_NEP_tropics_shyear_mean_beta*scale-EP_NEP_tropics_shyear_mean_beta_CI*scale EP_NEP_tropics_shyear_mean_beta*scale+EP_NEP_tropics_shyear_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(20+offs+offs/3,EP_NEP_tropics_shyear_mean_beta*scale,10,clr2(10,:),'filled','^');

% Calendar year
b1 = dberg_box(21-offs,CP_NEP_tropics_annual_beta*scale, clr2(4,:), 'w', wdth, 10);
b2 = dberg_box(21+offs,EP_NEP_tropics_annual_beta*scale, clr2(8,:), 'w', wdth, 10);
plot([21-offs+offs/3 21-offs+offs/3],...
    [CP_NEP_tropics_annual_mean_beta*scale-CP_NEP_tropics_annual_mean_beta_CI*scale CP_NEP_tropics_annual_mean_beta*scale+CP_NEP_tropics_annual_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(21-offs+offs/3,CP_NEP_tropics_annual_mean_beta*scale,10,clr2(2,:),'filled','^');
plot([21+offs+offs/3 21+offs+offs/3],...
    [EP_NEP_tropics_annual_mean_beta*scale-EP_NEP_tropics_annual_mean_beta_CI*scale EP_NEP_tropics_annual_mean_beta*scale+EP_NEP_tropics_annual_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(21+offs+offs/3,EP_NEP_tropics_annual_mean_beta*scale,10,clr2(10,:),'filled','^');

set(gca, 'YLim',[-1500 1000]/1000, 'YTick',-1.5:0.5:1);
ylb = ylabel('Annual NEP response (Pg C yr^{-1} K^{-1})', 'FontSize',7);

hold off;
box off;


%% Plot inversion beta through time: global
load('./data/cp_ep_nep_inversions.mat');
set(h, 'currentaxes',ax1);
yyaxis left;
hold on;
for i = 1:18
    
    plot([i-offs-offs/3 i-offs-offs/3],...
        [CP_NEP_global_monthly_mean_beta(i)-CP_NEP_global_monthly_mean_beta_CI(i) CP_NEP_global_monthly_mean_beta(i)+CP_NEP_global_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(2,:));
    scatter(i-offs-offs/3,CP_NEP_global_monthly_mean_beta(i),12,clr2(2,:),'filled','s');
    
    plot([i+offs-offs/3 i+offs-offs/3],...
        [EP_NEP_global_monthly_mean_beta(i)-EP_NEP_global_monthly_mean_beta_CI(i) EP_NEP_global_monthly_mean_beta(i)+EP_NEP_global_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(10,:));
    scatter(i+offs-offs/3,EP_NEP_global_monthly_mean_beta(i),12,clr2(10,:),'filled','s');
    
    
end

% prior Jul-Jun
yyaxis right;
plot([20-offs-offs/3 20-offs-offs/3],...
    [CP_NEP_global_shyear_mean_beta*scale-CP_NEP_global_shyear_mean_beta_CI*scale CP_NEP_global_shyear_mean_beta*scale+CP_NEP_global_shyear_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(20-offs-offs/3,CP_NEP_global_shyear_mean_beta*scale,12,clr2(2,:),'filled','s');
plot([20+offs-offs/3 20+offs-offs/3],...
    [EP_NEP_global_shyear_mean_beta*scale-EP_NEP_global_shyear_mean_beta_CI*scale EP_NEP_global_shyear_mean_beta*scale+EP_NEP_global_shyear_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(20+offs-offs/3,EP_NEP_global_shyear_mean_beta*scale,12,clr2(10,:),'filled','s');

% Calendar year
yyaxis right;
plot([21-offs-offs/3 21-offs-offs/3],...
    [CP_NEP_global_annual_mean_beta*scale-CP_NEP_global_annual_mean_beta_CI*scale CP_NEP_global_annual_mean_beta*scale+CP_NEP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(21-offs-offs/3,CP_NEP_global_annual_mean_beta*scale,12,clr2(2,:),'filled','s');
plot([21+offs-offs/3 21+offs-offs/3],...
    [EP_NEP_global_annual_mean_beta*scale-EP_NEP_global_annual_mean_beta_CI*scale EP_NEP_global_annual_mean_beta*scale+EP_NEP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(21+offs-offs/3,EP_NEP_global_annual_mean_beta*scale,12,clr2(10,:),'filled','s');

%% Plot inversion beta through time: tropics
load('./data/cp_ep_nep_inversions.mat');
set(h, 'currentaxes',ax2);
yyaxis left;
hold on;
for i = 1:18
    
    plot([i-offs-offs/3 i-offs-offs/3],...
        [CP_NEP_tropics_monthly_mean_beta(i)-CP_NEP_tropics_monthly_mean_beta_CI(i) CP_NEP_tropics_monthly_mean_beta(i)+CP_NEP_tropics_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(2,:));
    scatter(i-offs-offs/3,CP_NEP_tropics_monthly_mean_beta(i),12,clr2(2,:),'filled','s');
    
    plot([i+offs-offs/3 i+offs-offs/3],...
        [EP_NEP_tropics_monthly_mean_beta(i)-EP_NEP_tropics_monthly_mean_beta_CI(i) EP_NEP_tropics_monthly_mean_beta(i)+EP_NEP_tropics_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(10,:));
    scatter(i+offs-offs/3,EP_NEP_tropics_monthly_mean_beta(i),12,clr2(10,:),'filled','s');
    
    
end

% prior Jul-Jun
yyaxis right;
plot([20-offs-offs/3 20-offs-offs/3],...
    [CP_NEP_tropics_shyear_mean_beta*scale-CP_NEP_tropics_shyear_mean_beta_CI*scale CP_NEP_tropics_shyear_mean_beta*scale+CP_NEP_tropics_shyear_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(20-offs-offs/3,CP_NEP_tropics_shyear_mean_beta*scale,12,clr2(2,:),'filled','s');

plot([20+offs-offs/3 20+offs-offs/3],...
    [EP_NEP_tropics_shyear_mean_beta*scale-EP_NEP_tropics_shyear_mean_beta_CI*scale EP_NEP_tropics_shyear_mean_beta*scale+EP_NEP_tropics_shyear_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(20+offs-offs/3,EP_NEP_tropics_shyear_mean_beta*scale,12,clr2(10,:),'filled','s');

% Calendar year
yyaxis right;
plot([21-offs-offs/3 21-offs-offs/3],...
    [CP_NEP_tropics_annual_mean_beta*scale-CP_NEP_tropics_annual_mean_beta_CI*scale CP_NEP_tropics_annual_mean_beta*scale+CP_NEP_tropics_annual_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(21-offs-offs/3,CP_NEP_tropics_annual_mean_beta*scale,12,clr2(2,:),'filled','s');

plot([21+offs-offs/3 21+offs-offs/3],...
    [EP_NEP_tropics_annual_mean_beta*scale-EP_NEP_tropics_annual_mean_beta_CI*scale EP_NEP_tropics_annual_mean_beta*scale+EP_NEP_tropics_annual_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(21+offs-offs/3,EP_NEP_tropics_annual_mean_beta*scale,12,clr2(10,:),'filled','s');

%% Legend
scatter(12.5+0.1, -1.3, 20, 'k', 'filled','s');
text(12.6+0.1, -1.3, 'CAMS', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)
scatter(12.5+0.1, -1, 20, 'k', 'filled','^');
text(12.6+0.1, -1, 'MsTMIP', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)

lgd = legend([b1 b2], 'CP','EP', 'Location','southeast', 'FontSize',9);
legend('boxoff');
lgd.Position = [0.66 0.075 0.1012 0.07];


set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-epi-cpi-nep.tif')
close all;

