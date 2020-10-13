% Make some maps and stuff

clr = cbrewer('div','RdBu',12);

latlim = [-80 80];
lonlim = [-180 180];
worldland = shaperead('landareas','UseGeoCoords', true);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 0 7 8];

set(h, 'defaultAxesColorOrder',[0 0 0; 0,0,0]/255)

%% Map MsTMIP Beta
load('./data/cp_ep_gpp_mstmip.mat');

subplot(4,2,3)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, EP_GPP_annual_beta)
caxis([-0.06 0.06])
colormap(clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
text(-2.2,1.3,'C', 'FontSize',12);
text(-3.2,0,'MsTMIP', 'FontSize',12, 'FontWeight','bold','Rotation',90,...
    'HorizontalAlignment','center','VerticalAlignment','middle');

subplot(4,2,4)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(lat, lon, CP_GPP_annual_beta)
caxis([-0.06 0.06])
colormap(clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
set(gca, 'Position',pos);
text(-2.2,1.3,'D', 'FontSize',12);

cb = colorbar('eastoutside');
cb.Position = [0.51    0.57    0.04    0.38];
cb.Ticks = -0.06:0.01:0.06;
cb.TickLength = 0.09;
cb.TickLabels = {'-60','-50','-40','-30','-20','-10','0','10','20','30','40','50','60'};
ylb = ylabel(cb, 'Mean GPP response (g C m^{-2} yr^{-1} K^{-1})', 'FontSize',9);
ylb.Position = [-0.95 0.0000 0];

%% Plot MsTMIP beta through time
% Parameters
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

ax1 = subplot(4,2,[5 6]);
plot([0 15], [0 0], 'k-')
set(ax1, 'Position', [0.09 0.32 0.835 0.2])

% Monthly
yyaxis left;
hold on;
for i = 1:18
    
    dberg_box(i-offs,CP_GPP_global_monthly_beta(i,:), clr2(4,:), 'w', wdth, 10);
    dberg_box(i+offs,EP_GPP_global_monthly_beta(i,:), clr2(8,:), 'w', wdth, 10);
    
    plot([i-offs+offs/3 i-offs+offs/3],...
        [CP_GPP_global_monthly_mean_beta(i)-CP_GPP_global_monthly_mean_beta_CI(i) CP_GPP_global_monthly_mean_beta(i)+CP_GPP_global_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(2,:));
    scatter(i-offs+offs/3,CP_GPP_global_monthly_mean_beta(i),10,clr2(2,:),'filled','^');
    
    plot([i+offs+offs/3 i+offs+offs/3],...
        [EP_GPP_global_monthly_mean_beta(i)-EP_GPP_global_monthly_mean_beta_CI(i) EP_GPP_global_monthly_mean_beta(i)+EP_GPP_global_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(10,:));
    scatter(i+offs+offs/3,EP_GPP_global_monthly_mean_beta(i),10,clr2(10,:),'filled','^');
    
    
end
set(gca, 'XLim',[0.5 14], 'YLim',[-5 5], 'XTick',[1:12 13.5], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.01 0.05], 'XTickLabels','');
ylabel('Monthly GPP response (Tg C day^{-1} SD^{-1})', 'FontSize',7);
text(0.9, 4.5, 'E) Global', 'FontSize',12); 

% Annual
yyaxis right;
b1 = dberg_box(13.5-offs,CP_GPP_global_annual_beta*scale, clr2(4,:), 'w', wdth, 10);
b2 = dberg_box(13.5+offs,EP_GPP_global_annual_beta*scale, clr2(8,:), 'w', wdth, 10);

plot([13.5-offs+offs/3 13.5-offs+offs/3],...
    [CP_GPP_global_annual_mean_beta*scale-CP_GPP_global_annual_mean_beta_CI*scale CP_GPP_global_annual_mean_beta*scale+CP_GPP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(13.5-offs+offs/3,CP_GPP_global_annual_mean_beta*scale,10,clr2(2,:),'filled','^');

plot([13.5+offs+offs/3 13.5+offs+offs/3],...
    [EP_GPP_global_annual_mean_beta*scale-EP_GPP_global_annual_mean_beta_CI*scale EP_GPP_global_annual_mean_beta*scale+EP_GPP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(13.5+offs+offs/3,EP_GPP_global_annual_mean_beta*scale,10,clr2(10,:),'filled','^');
set(gca, 'YLim',[-1000 1000]/1000);
ylb = ylabel('Annual GPP response (Pg C yr^{-1} SD^{-1})', 'FontSize',7);
ylb.Position = [14.7    0.0000   -1.0000];

hold off;
box off;

%% Plot MsTMIP TROPICAL beta through time
% Parameters
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

ax2 = subplot(4,2,[7 8]);
plot([0 15], [0 0], 'k-')
set(ax2, 'Position', [0.09 0.07 0.835 0.2])

% Monthly
yyaxis left;
hold on;
for i = 1:12
    
    dberg_box(i-offs,CP_GPP_tropics_monthly_beta(i,:), clr2(4,:), 'w', wdth, 10);
    dberg_box(i+offs,EP_GPP_tropics_monthly_beta(i,:), clr2(8,:), 'w', wdth, 10);
    
    plot([i-offs+offs/3 i-offs+offs/3],...
        [CP_GPP_tropics_monthly_mean_beta(i)-CP_GPP_tropics_monthly_mean_beta_CI(i) CP_GPP_tropics_monthly_mean_beta(i)+CP_GPP_tropics_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(2,:));
    scatter(i-offs+offs/3,CP_GPP_tropics_monthly_mean_beta(i),10,clr2(2,:),'filled','^');
    
    plot([i+offs+offs/3 i+offs+offs/3],...
        [EP_GPP_tropics_monthly_mean_beta(i)-EP_GPP_tropics_monthly_mean_beta_CI(i) EP_GPP_tropics_monthly_mean_beta(i)+EP_GPP_tropics_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(10,:));
    scatter(i+offs+offs/3,EP_GPP_tropics_monthly_mean_beta(i),10,clr2(10,:),'filled','^');
    
    
end
set(gca, 'XLim',[0.5 14], 'YLim',[-4 4], 'XTick',[1:12 13.5], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.01 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
ylabel('Monthly GPP response (Tg C day^{-1} SD^{-1})', 'FontSize',7);
text(0.9, 3.5, 'F) Tropics', 'FontSize',12);

% Annual
yyaxis right;
b1 = dberg_box(13.5-offs,CP_GPP_tropics_annual_beta*scale, clr2(4,:), 'w', wdth, 10);
b2 = dberg_box(13.5+offs,EP_GPP_tropics_annual_beta*scale, clr2(8,:), 'w', wdth, 10);

plot([13.5-offs+offs/3 13.5-offs+offs/3],...
    [CP_GPP_tropics_annual_mean_beta*scale-CP_GPP_tropics_annual_mean_beta_CI*scale CP_GPP_tropics_annual_mean_beta*scale+CP_GPP_tropics_annual_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(13.5-offs+offs/3,CP_GPP_tropics_annual_mean_beta*scale,10,clr2(2,:),'filled','^');

plot([13.5+offs+offs/3 13.5+offs+offs/3],...
    [EP_GPP_tropics_annual_mean_beta*scale-EP_GPP_tropics_annual_mean_beta_CI*scale EP_GPP_tropics_annual_mean_beta*scale+EP_GPP_tropics_annual_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(13.5+offs+offs/3,EP_GPP_tropics_annual_mean_beta*scale,10,clr2(10,:),'filled','^');
set(gca, 'YLim',[-800 800]/1000, 'YTick',-0.8:0.4:0.8);
ylb = ylabel('Annual GPP response (Pg C yr^{-1} SD^{-1})', 'FontSize',7);
ylb.Position = [14.7    0.0000   -1.0000];

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
caxis([-0.06 0.06])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
pos(2) = pos(2)+0.05;
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
caxis([-0.06 0.06])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
pos(2) = pos(2)+0.05;
set(gca, 'Position',pos);
text(-2.2,1.3,'B', 'FontSize',12);
ttl = title('Central Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 1.75;

%% Plot LUE beta through time
set(h, 'currentaxes',ax1);

% CCW
% Monthly
yyaxis left;
hold on;
for i = 1:12
    
    plot([i-offs-offs/3 i-offs-offs/3],...
        [CP_GPP_global_monthly_mean_beta(i)-CP_GPP_global_monthly_mean_beta_CI(i) CP_GPP_global_monthly_mean_beta(i)+CP_GPP_global_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(2,:));
    scatter(i-offs-offs/3,CP_GPP_global_monthly_mean_beta(i),15,clr2(2,:),'x');
    
    plot([i+offs-offs/3 i+offs-offs/3],...
        [EP_GPP_global_monthly_mean_beta(i)-EP_GPP_global_monthly_mean_beta_CI(i) EP_GPP_global_monthly_mean_beta(i)+EP_GPP_global_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(10,:));
    scatter(i+offs-offs/3,EP_GPP_global_monthly_mean_beta(i),15,clr2(10,:),'x');
    
    
end

% Annual
yyaxis right;
plot([13.5-offs-offs/3 13.5-offs-offs/3],...
    [CP_GPP_global_annual_mean_beta*scale-CP_GPP_global_annual_mean_beta_CI*scale CP_GPP_global_annual_mean_beta*scale+CP_GPP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(13.5-offs-offs/3,CP_GPP_global_annual_mean_beta*scale,15,clr2(2,:),'x');

plot([13.5+offs-offs/3 13.5+offs-offs/3],...
    [EP_GPP_global_annual_mean_beta*scale-EP_GPP_global_annual_mean_beta_CI*scale EP_GPP_global_annual_mean_beta*scale+EP_GPP_global_annual_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(13.5+offs-offs/3,EP_GPP_global_annual_mean_beta*scale,15,clr2(10,:),'x');

%% Plot TROPICAL LUE beta through time
set(h, 'currentaxes',ax2);

% CCW
% Monthly
yyaxis left;
hold on;
for i = 1:12
    
    plot([i-offs-offs/3 i-offs-offs/3],...
        [CP_GPP_tropics_monthly_mean_beta(i)-CP_GPP_tropics_monthly_mean_beta_CI(i) CP_GPP_tropics_monthly_mean_beta(i)+CP_GPP_tropics_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(2,:));
    scatter(i-offs-offs/3,CP_GPP_tropics_monthly_mean_beta(i),15,clr2(2,:),'x');
    
    plot([i+offs-offs/3 i+offs-offs/3],...
        [EP_GPP_tropics_monthly_mean_beta(i)-EP_GPP_tropics_monthly_mean_beta_CI(i) EP_GPP_tropics_monthly_mean_beta(i)+EP_GPP_tropics_monthly_mean_beta_CI(i)],...
        '-','Color',clr2(10,:));
    scatter(i+offs-offs/3,EP_GPP_tropics_monthly_mean_beta(i),15,clr2(10,:),'x');
    
    
end

% Annual
yyaxis right;
plot([13.5-offs-offs/3 13.5-offs-offs/3],...
    [CP_GPP_tropics_annual_mean_beta*scale-CP_GPP_tropics_annual_mean_beta_CI*scale CP_GPP_tropics_annual_mean_beta*scale+CP_GPP_tropics_annual_mean_beta_CI*scale],...
    '-','Color',clr2(2,:));
scatter(13.5-offs-offs/3,CP_GPP_tropics_annual_mean_beta*scale,15,clr2(2,:),'x');

plot([13.5+offs-offs/3 13.5+offs-offs/3],...
    [EP_GPP_tropics_annual_mean_beta*scale-EP_GPP_tropics_annual_mean_beta_CI*scale EP_GPP_tropics_annual_mean_beta*scale+EP_GPP_tropics_annual_mean_beta_CI*scale],...
    '-','Color',clr2(10,:));
scatter(13.5+offs-offs/3,EP_GPP_tropics_annual_mean_beta*scale,15,clr2(10,:),'x');


%% Legend
scatter(11-0.3, 0.78, 20, 'k', 'filled','^');
text(11.1-0.3, 0.78, 'MsTMIP', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)
scatter(11-0.3, 0.58, 20, 'k', 'x');
text(11.1-0.3, 0.58, 'LUE', 'HorizontalAlignment','left', 'VerticalAlignment','middle', 'FontSize',9)

lgd = legend([b1 b2], 'CP','EP', 'Location','northwest', 'FontSize',9);
legend('boxoff');
lgd.Position = [0.81    0.22    0.1012    0.0634];

%% Save
set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/epi-cpi-gpp.tif')
close all;

