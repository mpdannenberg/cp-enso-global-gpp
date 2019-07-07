% Make some maps and stuff

clr = [103,0,31
178,24,43
214,96,77
244,165,130
253,219,199
209,229,240
146,197,222
67,147,195
33,102,172
5,48,97]/255;

latlim = [-80 80];
lonlim = [-180 180];
worldland = shaperead('landareas','UseGeoCoords', true);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 4];

set(h, 'defaultAxesColorOrder',[0 0 0; 0 0 0]/255)

%% Map MsTMIP Beta
load('./data/cp_ep_nep_mstmip.mat');

subplot(2,2,1)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, EP_NEP_annual_beta)
caxis([-0.025 0.025])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
text(-2.2,1.3,'A', 'FontSize',12);
ttl = title('Eastern Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 1.75;

subplot(2,2,2)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, CP_NEP_annual_beta)
caxis([-0.025 0.025])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
set(gca, 'Position',pos);
text(-2.2,1.3,'B', 'FontSize',12);
ttl = title('Central Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 1.75;

cb = colorbar('eastoutside');
cb.Position = [0.52 0.6 0.04 0.3];
cb.Ticks = -0.025:0.005:0.025;
cb.TickLength = 0.23;
cb.TickLabels = {'-25','','','','','0','','','','','25'};
% cb.YDir = 'reverse';
ylb = ylabel(cb, {'Mean NEP response'; '(g C m^{-2} yr^{-1} SD^{-1})'});
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

ax1 = subplot(2,2,[3 4]);
plot([0 15], [0 0], 'k-')
set(ax1, 'Position', [0.09 0.07 0.835 0.45])

% Monthly
yyaxis left;
hold on;
for i = 1:12
    
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
set(gca, 'XLim',[0.5 14], 'YLim',[-5 3], 'XTick',[1:12 13.5], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.01 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
ylb = ylabel('Monthly NEP response (Tg C day^{-1} SD^{-1})', 'FontSize',7);
ylb.Position = [-0.1 -1 -1];
text(0.9, 2.75, 'C', 'FontSize',12);

% Annual
yyaxis right;
b1 = dberg_box(13.5-offs,CP_NEP_global_annual_beta*scale, clr2(4,:), 'w', wdth, 10);
b2 = dberg_box(13.5+offs,EP_NEP_global_annual_beta*scale, clr2(8,:), 'w', wdth, 10);

plot([13.5-offs+offs/3 13.5-offs+offs/3],...
    [CP_NEP_global_annual_mean_beta*scale-CP_NEP_global_annual_mean_beta_CI*scale CP_NEP_global_annual_mean_beta*scale+CP_NEP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(13.5-offs+offs/3,CP_NEP_global_annual_mean_beta*scale,10,clr2(2,:),'filled','^');

plot([13.5+offs+offs/3 13.5+offs+offs/3],...
    [EP_NEP_global_annual_mean_beta*scale-EP_NEP_global_annual_mean_beta_CI*scale EP_NEP_global_annual_mean_beta*scale+EP_NEP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(13.5+offs+offs/3,EP_NEP_global_annual_mean_beta*scale,10,clr2(10,:),'filled','^');
set(gca, 'YLim',[-1250 750]/1000);
ylb = ylabel('Annual NEP response (Pg C yr^{-1} SD^{-1})', 'FontSize',7);
ylb.Position = [14.7    -0.25   -1.0000];

hold off;
box off;


%% Plot inversion beta through time
load('./data/cp_ep_nep_inversions.mat');
yyaxis left;
hold on;
for i = 1:12
    
    plot([i-offs-offs/3 i-offs-offs/3],...
        [CP_NEP_global_monthly_mean_beta(i)-CP_NEP_global_monthly_mean_beta_CI(i) CP_NEP_global_monthly_mean_beta(i)+CP_NEP_global_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(2,:));
    scatter(i-offs-offs/3,CP_NEP_global_monthly_mean_beta(i),12,clr2(2,:),'filled','s');
    
    plot([i+offs-offs/3 i+offs-offs/3],...
        [EP_NEP_global_monthly_mean_beta(i)-EP_NEP_global_monthly_mean_beta_CI(i) EP_NEP_global_monthly_mean_beta(i)+EP_NEP_global_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(10,:));
    scatter(i+offs-offs/3,EP_NEP_global_monthly_mean_beta(i),12,clr2(10,:),'filled','s');
    
    
end

% Annual
yyaxis right;
plot([13.5-offs-offs/3 13.5-offs-offs/3],...
    [CP_NEP_global_annual_mean_beta*scale-CP_NEP_global_annual_mean_beta_CI*scale CP_NEP_global_annual_mean_beta*scale+CP_NEP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(13.5-offs-offs/3,CP_NEP_global_annual_mean_beta*scale,12,clr2(2,:),'filled','s');

plot([13.5+offs-offs/3 13.5+offs-offs/3],...
    [EP_NEP_global_annual_mean_beta*scale-EP_NEP_global_annual_mean_beta_CI*scale EP_NEP_global_annual_mean_beta*scale+EP_NEP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(13.5+offs-offs/3,EP_NEP_global_annual_mean_beta*scale,12,clr2(10,:),'filled','s');

%% Legend
scatter(11, -1.05, 20, 'k', 'filled','s');
text(11.1, -1.05, 'Inversions', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)
scatter(11, -0.85, 20, 'k', 'filled','^');
text(11.1, -0.85, 'MsTMIP', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)

lgd = legend([b1 b2], 'CP','EP', 'Location','southeast', 'FontSize',9);
legend('boxoff');
lgd.Position = [0.6 0.09 0.1012 0.0951];


set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/epi-cpi-nep.tif')
close all;

