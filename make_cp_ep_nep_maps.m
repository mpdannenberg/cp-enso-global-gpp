% Make some maps and stuff

clr = [84,48,5
    140,81,10
    191,129,45
    223,194,125
    246,232,195
    199,234,229
    128,205,193
    53,151,143
    1,102,94
    0,60,48]/255;

latlim = [-80 80];
lonlim = [-180 180];
worldland = shaperead('landareas','UseGeoCoords', true);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 4];

set(h, 'defaultAxesColorOrder',[0 0 0; 140,81,10]/255)


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

cb = colorbar('eastoutside');
cb.Position = [0.5    0.6    0.04    0.3];
cb.Ticks = -0.025:0.005:0.025;
cb.TickLength = 0.24;
% cb.TickLabels = {'-0.025','','','','','0','','','','','0.025'};
% ylb = ylabel(cb, 'kg C m^{-2} yr^{-1}');
cb.TickLabels = {'-25','','','','','0','','','','','25'};
ylb = ylabel(cb, 'g C m^{-2} yr^{-1}');
ylb.Position = [-0.75 0.0000 0];

%% Plot MsTMIP beta through time
ax1 = subplot(2,2,3);
yyaxis left;
for i = 1:12
    lower = min(EP_NEP_global_monthly_beta(i, :));
    upper = max(EP_NEP_global_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [EP_NEP_global_monthly_mean_beta(i) EP_NEP_global_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
text(1, 3, 'C', 'FontSize',12);
yyaxis right;
lower = min(EP_NEP_global_annual_beta)/1000;
upper = max(EP_NEP_global_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [223,194,125]/255, 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [EP_NEP_global_annual_mean_beta EP_NEP_global_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color',[140,81,10]/255)
set(gca, 'YLim',[-1000 1000]/1000, 'FontSize',8);
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

ax2 = subplot(2,2,4);
yyaxis left;
for i = 1:12
    lower = min(CP_NEP_global_monthly_beta(i, :));
    upper = max(CP_NEP_global_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [CP_NEP_global_monthly_mean_beta(i) CP_NEP_global_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-4 4], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
set(gca, 'Position',pos);
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
text(1, 3, 'D', 'FontSize',12);
yyaxis right;
lower = min(CP_NEP_global_annual_beta)/1000;
upper = max(CP_NEP_global_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [223,194,125]/255, 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [CP_NEP_global_annual_mean_beta CP_NEP_global_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color',[140,81,10]/255)
set(gca, 'YLim',[-1000 1000]/1000);
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

%% Plot inversion beta through time
load('./data/cp_ep_nep_inversions.mat');

clr2 = [166,206,227
31,120,180
178,223,138
51,160,44]/255;

set(h, 'currentaxes',ax1);
yyaxis left;
for i = 1:12
    for j = 1:4
        hold on;
        plot([i-0.4 i+0.4], [EP_NEP_global_monthly_beta(i,j) EP_NEP_global_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [EP_NEP_global_annual_beta(1) EP_NEP_global_annual_beta(1)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [EP_NEP_global_annual_beta(2) EP_NEP_global_annual_beta(2)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(2,:))

% Manual legend
i = 4;
plot([i-0.4 i+0.4], [0.85 0.85],...
    '-', 'LineWidth',3, 'Color',clr2(1,:));
text(i+0.5, 0.85, 'CAMS v15r2', 'Color',clr2(1,:), 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')
plot([i-0.4 i+0.4], [0.65 0.65],...
    '-', 'LineWidth',3, 'Color',clr2(2,:))
text(i+0.5, 0.65, 'CAMS v15r4', 'Color',clr2(2,:), 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')
i = 9.5;
plot([i-0.4 i+0.4], [0.85 0.85],...
    '-', 'LineWidth',3, 'Color',clr2(3,:));
text(i+0.5, 0.85, 'MACC v13r1', 'Color',clr2(3,:), 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')
plot([i-0.4 i+0.4], [0.65 0.65],...
    '-', 'LineWidth',3, 'Color',clr2(4,:))
text(i+0.5, 0.65, 'MACC v14r2', 'Color',clr2(4,:), 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')

set(h, 'currentaxes',ax2); 
yyaxis left;
for i = 1:12
    for j = 1:4
        hold on;
        plot([i-0.4 i+0.4], [CP_NEP_global_monthly_beta(i,j) CP_NEP_global_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [CP_NEP_global_annual_beta(1) CP_NEP_global_annual_beta(1)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [CP_NEP_global_annual_beta(2) CP_NEP_global_annual_beta(2)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(2,:))

% Manual legend
i = 9.2;
fill([i-0.3 i+0.3 i+0.3 i-0.3], [0.45 0.45 0.95 0.95],...
    [0.8 0.8 0.8], 'EdgeColor','none');
plot([i-0.3 i+0.3], [0.7 0.7],...
    'k-', 'LineWidth',3)
text(i+0.5, 0.7, 'MsTMIP mean', 'Color','k', 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')
text(i+0.5, 0.9, 'MsTMIP range', 'Color',[0.7 0.7 0.7], 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/epi-cpi-gpp.tif')
close all;
