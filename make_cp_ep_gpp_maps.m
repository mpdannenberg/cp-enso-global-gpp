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
h.Position = [1 1 7 6];

set(h, 'defaultAxesColorOrder',[0 0 0; 0,0,0]/255)

%% Map MsTMIP Beta
load('./data/cp_ep_gpp_mstmip.mat');

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
colormap(clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
text(-2.2,1.3,'C', 'FontSize',12);
text(-3.2,0,'MsTMIP', 'FontSize',12, 'FontWeight','bold','Rotation',90,...
    'HorizontalAlignment','center','VerticalAlignment','middle');

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
colormap(clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
set(gca, 'Position',pos);
text(-2.2,1.3,'D', 'FontSize',12);

cb = colorbar('eastoutside');
cb.Position = [0.51    0.41    0.04    0.52];
cb.Ticks = -0.05:0.01:0.05;
cb.TickLength = 0.09;
% cb.TickLabels = {'-0.05','','','','','0','','','','','0.05'};
% ylb = ylabel(cb, 'kg C m^{-2} yr^{-1}');
% cb.TickLabels = {'-50','','','','','0','','','','','50'};
cb.TickLabels = {'-50','-40','-30','-20','-10','0','10','20','30','40','50'};
ylb = ylabel(cb, 'Mean GPP response (g C m^{-2} yr^{-1} SD^{-1})', 'FontSize',9);
ylb.Position = [-0.95 0.0000 0];

%% Plot MsTMIP beta through time
% ax1 = subplot(3,2,5);
% yyaxis left;
% for i = 1:12
%     lower = min(EP_GPP_global_monthly_beta(i, :));
%     upper = max(EP_GPP_global_monthly_beta(i, :));
%     fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
%         [0.8 0.8 0.8], 'EdgeColor','none');
%     hold on;
%     plot([i-0.4 i+0.4], [EP_GPP_global_monthly_mean_beta(i) EP_GPP_global_monthly_mean_beta(i)],...
%         'k-', 'LineWidth',3)
% end
% set(gca, 'XLim',[0 15], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
%     'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
% pos = get(gca, 'Position');
% pos(1) = pos(1)-0.04;
% set(gca, 'Position',pos);
% hold off;
% ylabel('Tg C day^{-1}', 'FontSize',8);
% text(1, 3, 'E', 'FontSize',12);
% yyaxis right;
% lower = min(EP_GPP_global_annual_beta)/1000;
% upper = max(EP_GPP_global_annual_beta)/1000;
% i = 14;
% fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
%     [0.8 0.8 0.8], 'EdgeColor','none');
% hold on;
% plot([i-0.4 i+0.4], [EP_GPP_global_annual_mean_beta EP_GPP_global_annual_mean_beta]/1000,...
%     '-', 'LineWidth',3, 'Color','k')
% set(gca, 'YLim',[-1000 1000]/1000, 'FontSize',8);
% plot([0 15],[0 0],'k-')
% ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
% ylb.Position = [16.8    0.0000   -1.0000];
% 
% ax2 = subplot(3,2,6);
% yyaxis left;
% for i = 1:12
%     lower = min(CP_GPP_global_monthly_beta(i, :));
%     upper = max(CP_GPP_global_monthly_beta(i, :));
%     fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
%         [0.8 0.8 0.8], 'EdgeColor','none');
%     hold on;
%     plot([i-0.4 i+0.4], [CP_GPP_global_monthly_mean_beta(i) CP_GPP_global_monthly_mean_beta(i)],...
%         'k-', 'LineWidth',3)
% end
% set(gca, 'XLim',[0 15], 'YLim',[-4 4], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
%     'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
% pos = get(gca, 'Position');
% pos(1) = pos(1)+0.02;
% set(gca, 'Position',pos);
% hold off;
% ylabel('Tg C day^{-1}', 'FontSize',8);
% text(1, 3, 'F', 'FontSize',12);
% yyaxis right;
% lower = min(CP_GPP_global_annual_beta)/1000;
% upper = max(CP_GPP_global_annual_beta)/1000;
% i = 14;
% fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
%     [0.8 0.8 0.8], 'EdgeColor','none');
% hold on;
% plot([i-0.4 i+0.4], [CP_GPP_global_annual_mean_beta CP_GPP_global_annual_mean_beta]/1000,...
%     '-', 'LineWidth',3, 'Color','k')
% set(gca, 'YLim',[-1000 1000]/1000);
% plot([0 15],[0 0],'k-')
% ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
% ylb.Position = [16.8    0.0000   -1.0000];

% Parameters
wdth = 0.22;
offs = 0.15;
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

ax1 = subplot(3,2,[5 6]);
plot([0 15], [0 0], 'k-')
set(ax1, 'Position', [0.09 0.07 0.835 0.3])

% Monthly
yyaxis left;
hold on;
for i = 1:12
    
    dberg_box(i-offs,CP_GPP_global_monthly_beta(i,:), clr2(4,:), 'w', wdth);
    dberg_box(i+offs,EP_GPP_global_monthly_beta(i,:), clr2(8,:), 'w', wdth);
    
    plot([i-offs i-offs],...
        [CP_GPP_global_monthly_mean_beta(i)-CP_GPP_global_monthly_mean_beta_CI(i) CP_GPP_global_monthly_mean_beta(i)+CP_GPP_global_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(2,:));
    scatter(i-offs,CP_GPP_global_monthly_mean_beta(i),10,clr2(2,:),'filled','^');
    
    plot([i+offs i+offs],...
        [EP_GPP_global_monthly_mean_beta(i)-EP_GPP_global_monthly_mean_beta_CI(i) EP_GPP_global_monthly_mean_beta(i)+EP_GPP_global_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(10,:));
    scatter(i+offs,EP_GPP_global_monthly_mean_beta(i),10,clr2(10,:),'filled','^');
    
    
end
set(gca, 'XLim',[0.5 14.5], 'YLim',[-5 5], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.01 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
ylabel('Monthly GPP response (Tg C day^{-1} SD^{-1})', 'FontSize',7);
text(0.9, 4.5, 'E', 'FontSize',12);

% Annual
yyaxis right;
b1 = dberg_box(14-offs,CP_GPP_global_annual_beta*scale, clr2(4,:), 'w', wdth);
b2 = dberg_box(14+offs,EP_GPP_global_annual_beta*scale, clr2(8,:), 'w', wdth);

plot([14-offs 14-offs],...
    [CP_GPP_global_annual_mean_beta*scale-CP_GPP_global_annual_mean_beta_CI*scale CP_GPP_global_annual_mean_beta*scale+CP_GPP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(14-offs,CP_GPP_global_annual_mean_beta*scale,10,clr2(2,:),'filled','^');

plot([14+offs 14+offs],...
    [EP_GPP_global_annual_mean_beta*scale-EP_GPP_global_annual_mean_beta_CI*scale EP_GPP_global_annual_mean_beta*scale+EP_GPP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(14+offs,EP_GPP_global_annual_mean_beta*scale,10,clr2(10,:),'filled','^');
set(gca, 'YLim',[-1000 1000]/1000);
ylb = ylabel('Annual GPP response (Pg C yr^{-1} SD^{-1})', 'FontSize',7);
ylb.Position = [15.2    0.0000   -1.0000];

hold off;
box off;

%% Map LUE Beta
load('./data/cp_ep_gpp_lue.mat');
CP_GPP_annual_beta(CP_GPP_annual_beta==0) = NaN;
EP_GPP_annual_beta(EP_GPP_annual_beta==0) = NaN;

subplot(3,2,1)
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
text(-2.2,1.3,'A', 'FontSize',12);
ttl = title('Eastern Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 1.75;
text(-3.2,0,'LUE', 'FontSize',12, 'FontWeight','bold','Rotation',90,...
    'HorizontalAlignment','center','VerticalAlignment','middle');

subplot(3,2,2)
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
text(-2.2,1.3,'B', 'FontSize',12);
ttl = title('Central Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 1.75;

%% Plot LUE beta through time

% clr2 = [55,126,184
% 77,175,74]/255;
% 
% set(h, 'currentaxes',ax1);
% yyaxis left;
% for i = 1:12
%     for j = 1:2
%         hold on;
%         plot([i-0.4 i+0.4], [EP_GPP_global_monthly_beta(i,j) EP_GPP_global_monthly_beta(i,j)],...
%             'k-', 'LineWidth',3, 'Color',clr2(j,:))
%     end
% end
% yyaxis right;
% i = 14;
% plot([i-0.4 i+0.4], [EP_GPP_global_annual_beta(1) EP_GPP_global_annual_beta(1)]/1000,...
%     '-', 'LineWidth',3, 'Color',clr2(1,:))
% plot([i-0.4 i+0.4], [EP_GPP_global_annual_beta(2) EP_GPP_global_annual_beta(2)]/1000,...
%     '-', 'LineWidth',3, 'Color',clr2(2,:))
% 
% set(h, 'currentaxes',ax2); 
% yyaxis left;
% for i = 1:12
%     for j = 1:2
%         hold on;
%         plot([i-0.4 i+0.4], [CP_GPP_global_monthly_beta(i,j) CP_GPP_global_monthly_beta(i,j)],...
%             'k-', 'LineWidth',3, 'Color',clr2(j,:))
%     end
% end
% yyaxis right;
% i = 14;
% plot([i-0.4 i+0.4], [CP_GPP_global_annual_beta(1) CP_GPP_global_annual_beta(1)]/1000,...
%     '-', 'LineWidth',3, 'Color',clr2(1,:))
% plot([i-0.4 i+0.4], [CP_GPP_global_annual_beta(2) CP_GPP_global_annual_beta(2)]/1000,...
%     '-', 'LineWidth',3, 'Color',clr2(2,:))
% 
% % Manual legend
% i = 3.5;
% plot([i-0.4 i+0.4], [0.85 0.85],...
%     '-', 'LineWidth',3, 'Color',clr2(1,:));
% text(i+0.5, 0.85, 'CCW', 'Color',clr2(1,:), 'FontWeight','bold', 'FontSize',8,...
%     'VerticalAlignment','middle', 'HorizontalAlignment','left')
% plot([i-0.4 i+0.4], [0.65 0.65],...
%     '-', 'LineWidth',3, 'Color',clr2(2,:))
% text(i+0.5, 0.65, 'MOD17', 'Color',clr2(2,:), 'FontWeight','bold', 'FontSize',8,...
%     'VerticalAlignment','middle', 'HorizontalAlignment','left')
% 
% i = 8.5;
% fill([i-0.3 i+0.3 i+0.3 i-0.3], [0.45 0.45 0.95 0.95],...
%     [0.8 0.8 0.8], 'EdgeColor','none');
% plot([i-0.3 i+0.3], [0.7 0.7],...
%     'k-', 'LineWidth',3)
% text(i+0.5, 0.7, 'MsTMIP mean', 'Color','k', 'FontWeight','bold', 'FontSize',8,...
%     'VerticalAlignment','middle', 'HorizontalAlignment','left')
% text(i+0.5, 0.9, 'MsTMIP range', 'Color',[0.7 0.7 0.7], 'FontWeight','bold', 'FontSize',8,...
%     'VerticalAlignment','middle', 'HorizontalAlignment','left')

set(h, 'currentaxes',ax1);

% CCW
% Monthly
yyaxis left;
hold on;
for i = 1:12
    
    plot([i-offs+offs/2 i-offs+offs/2],...
        [CP_GPP_global_monthly_beta(i,1)-CP_GPP_global_monthly_beta_CI(i,1) CP_GPP_global_monthly_beta(i,1)+CP_GPP_global_monthly_beta_CI(i,1)],...
        '-','Color',clr2(2,:));
    scatter(i-offs+offs/2,CP_GPP_global_monthly_beta(i,1),15,clr2(2,:),'x');
    
    plot([i+offs+offs/2 i+offs+offs/2],...
        [EP_GPP_global_monthly_beta(i,1)-EP_GPP_global_monthly_beta_CI(i,1) EP_GPP_global_monthly_beta(i,1)+EP_GPP_global_monthly_beta_CI(i,1)],...
        '-','Color',clr2(10,:));
    scatter(i+offs+offs/2,EP_GPP_global_monthly_beta(i,1),15,clr2(10,:),'x');
    
    
end

% Annual
yyaxis right;
plot([14-offs+offs/2 14-offs+offs/2],...
    [CP_GPP_global_annual_beta(1)*scale-CP_GPP_global_annual_beta_CI(1)*scale CP_GPP_global_annual_beta(1)*scale+CP_GPP_global_annual_beta_CI(1)*scale],...
    '-','Color',clr2(2,:));
scatter(14-offs+offs/2,CP_GPP_global_annual_beta(1)*scale,15,clr2(2,:),'x');

plot([14+offs+offs/2 14+offs+offs/2],...
    [EP_GPP_global_annual_beta(1)*scale-EP_GPP_global_annual_beta_CI(1)*scale EP_GPP_global_annual_beta(1)*scale+EP_GPP_global_annual_beta_CI(1)*scale],...
    '-','Color',clr2(10,:));
scatter(14+offs+offs/2,EP_GPP_global_annual_beta(1)*scale,15,clr2(10,:),'x');

% MOD17
% Monthly
yyaxis left;
hold on;
for i = 1:12
    
    plot([i-offs-offs/2 i-offs-offs/2],...
        [CP_GPP_global_monthly_beta(i,2)-CP_GPP_global_monthly_beta_CI(i,2) CP_GPP_global_monthly_beta(i,2)+CP_GPP_global_monthly_beta_CI(i,2)],...
        '-','Color',clr2(2,:));
    scatter(i-offs-offs/2,CP_GPP_global_monthly_beta(i,2),12,clr2(2,:),'filled','s');
    
    plot([i+offs-offs/2 i+offs-offs/2],...
        [EP_GPP_global_monthly_beta(i,2)-EP_GPP_global_monthly_beta_CI(i,2) EP_GPP_global_monthly_beta(i,2)+EP_GPP_global_monthly_beta_CI(i,2)],...
        '-','Color',clr2(10,:));
    scatter(i+offs-offs/2,EP_GPP_global_monthly_beta(i,2),12,clr2(10,:),'filled','s');
    
    
end

% Annual
yyaxis right;
plot([14-offs-offs/2 14-offs-offs/2],...
    [CP_GPP_global_annual_beta(2)*scale-CP_GPP_global_annual_beta_CI(2)*scale CP_GPP_global_annual_beta(2)*scale+CP_GPP_global_annual_beta_CI(2)*scale],...
    '-','Color',clr2(2,:));
scatter(14-offs-offs/2,CP_GPP_global_annual_beta(2)*scale,12,clr2(2,:),'filled','s');

plot([14+offs-offs/2 14+offs-offs/2],...
    [EP_GPP_global_annual_beta(2)*scale-EP_GPP_global_annual_beta_CI(2)*scale EP_GPP_global_annual_beta(2)*scale+EP_GPP_global_annual_beta_CI(2)*scale],...
    '-','Color',clr2(10,:));
scatter(14+offs-offs/2,EP_GPP_global_annual_beta(2)*scale,12,clr2(10,:),'filled','s');

%% Legend
scatter(2, 0.9, 20, 'k', 'filled','s');
text(2.1, 0.9, 'MOD17', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)
scatter(2, 0.75, 20, 'k', 'filled','^');
text(2.1, 0.75, 'MsTMIP', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)
scatter(2, 0.6, 20, 'k', 'x');
text(2.1, 0.6, 'CCW', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)

lgd = legend([b1 b2], 'CP','EP', 'Location','northwest', 'FontSize',9);
legend('boxoff');
lgd.Position = [0.27    0.3    0.1012    0.0634];

%% Save
set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/epi-cpi-gpp.tif')
close all;

