%% Make some maps and stuff
load('./data/cp_ep_gpp_mstmip.mat');

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
h.Position = [1 1 7 6];

set(h, 'defaultAxesColorOrder',[0 0 0; 140,81,10]/255)

% Correlations
subplot(3,2,1)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, EP_GPP_annual_r)
caxis([-1 1])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
ttl = title('Eastern Pacific Index','FontSize',14);
ttl.Position = [0.0000    1.9    0.0000];
text(-2.2,1.3,'A', 'FontSize',12);

subplot(3,2,2)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, CP_GPP_annual_r)
caxis([-1 1])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
set(gca, 'Position',pos);
ttl = title('Central Pacific Index','FontSize',14);
ttl.Position = [0.0000    1.9    0.0000];
text(-2.2,1.3,'B', 'FontSize',12);

cb = colorbar('eastoutside');
cb.Position = [0.5    0.71    0.04    0.22];
cb.Ticks = -1:0.2:1;
cb.TickLength = 0.21;
cb.TickLabels = {'-1','','','','','0','','','','','1'};
ylabel(cb, 'Pearson''s correlation (R)');

% Beta
subplot(3,2,3)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, EP_GPP_annual_beta)
caxis([-0.05 0.05])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
text(-2.2,1.3,'C', 'FontSize',12);

subplot(3,2,4)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, CP_GPP_annual_beta)
caxis([-0.05 0.05])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
set(gca, 'Position',pos);
text(-2.2,1.3,'D', 'FontSize',12);

cb = colorbar('eastoutside');
cb.Position = [0.5    0.41    0.04    0.22];
cb.Ticks = -0.05:0.01:0.05;
cb.TickLength = 0.21;
cb.TickLabels = {'-0.05','','','','','0','','','','','0.05'};
ylb = ylabel(cb, 'kg C m^{-2} yr^{-1}');
ylb.Position = [-0.75 0.0000 0];

% Plot beta through time
subplot(3,2,5)
yyaxis left;
for i = 1:12
    lower = min(EP_GPP_global_monthly_beta(i, :));
    upper = max(EP_GPP_global_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [EP_GPP_global_monthly_mean_beta(i) EP_GPP_global_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
text(1, 3, 'E', 'FontSize',12);
yyaxis right;
lower = min(EP_GPP_global_annual_beta)/1000;
upper = max(EP_GPP_global_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [223,194,125]/255, 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [EP_GPP_global_annual_mean_beta EP_GPP_global_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color',[140,81,10]/255)
set(gca, 'YLim',[-1000 1000]/1000, 'FontSize',8);
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

subplot(3,2,6)
yyaxis left;
for i = 1:12
    lower = min(CP_GPP_global_monthly_beta(i, :));
    upper = max(CP_GPP_global_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [CP_GPP_global_monthly_mean_beta(i) CP_GPP_global_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-4 4], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
set(gca, 'Position',pos);
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
text(1, 3, 'F', 'FontSize',12);
yyaxis right;
lower = min(CP_GPP_global_annual_beta)/1000;
upper = max(CP_GPP_global_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [223,194,125]/255, 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [CP_GPP_global_annual_mean_beta CP_GPP_global_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color',[140,81,10]/255)
set(gca, 'YLim',[-1000 1000]/1000);
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/epi-cpi-mstmip.tif')
close all;

