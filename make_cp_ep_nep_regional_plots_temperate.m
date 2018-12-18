% Make some plots

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 0 7 8];

set(h, 'defaultAxesColorOrder',[0 0 0; 0,0,0]/255)


%% Plot MsTMIP beta through time
% Central Asia
load('./data/cp_ep_nep_mstmip_regional.mat');
ax1 = subplot(4,2,1);
yyaxis left;
for i = 1:12
    lower = min(EP_NEP_casia_monthly_beta(i, :));
    upper = max(EP_NEP_casia_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [EP_NEP_casia_monthly_mean_beta(i) EP_NEP_casia_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
pos(2) = pos(2)+0.02;
set(gca, 'Position',pos, 'YDir','reverse');
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
text(1, 0.75, 'Central Asia', 'FontSize',10);
yyaxis right;
lower = min(EP_NEP_casia_annual_beta)/1000;
upper = max(EP_NEP_casia_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [0.8 0.8 0.8], 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [EP_NEP_casia_annual_mean_beta EP_NEP_casia_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color','k')
set(gca, 'YLim',[-200 200]/1000, 'FontSize',8, 'YDir','reverse');
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];
ttl = title('Eastern Pacific ENSO', 'FontSize',11);
ttl.Position(2) = -1.25;

ax2 = subplot(4,2,2);
yyaxis left;
for i = 1:12
    lower = min(CP_NEP_casia_monthly_beta(i, :));
    upper = max(CP_NEP_casia_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [CP_NEP_casia_monthly_mean_beta(i) CP_NEP_casia_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
pos(2) = pos(2)+0.02;
set(gca, 'Position',pos, 'YDir','reverse');
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
% text(1, 3, 'F', 'FontSize',12);
yyaxis right;
lower = min(CP_NEP_casia_annual_beta)/1000;
upper = max(CP_NEP_casia_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [0.8 0.8 0.8], 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [CP_NEP_casia_annual_mean_beta CP_NEP_casia_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color','k')
set(gca, 'YLim',[-200 200]/1000, 'YDir','reverse');
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];
ttl = title('Central Pacific ENSO', 'FontSize',11);
ttl.Position(2) = -1.25;

% Eastern U.S.
load('./data/cp_ep_nep_mstmip_regional.mat');
ax5 = subplot(4,2,3);
yyaxis left;
for i = 1:12
    lower = min(EP_NEP_eastus_monthly_beta(i, :));
    upper = max(EP_NEP_eastus_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [EP_NEP_eastus_monthly_mean_beta(i) EP_NEP_eastus_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
pos(2) = pos(2)-0.0;
set(gca, 'Position',pos, 'YDir','reverse');
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
text(1, 0.75, 'Eastern U.S.', 'FontSize',10);
yyaxis right;
lower = min(EP_NEP_eastus_annual_beta)/1000;
upper = max(EP_NEP_eastus_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [0.8 0.8 0.8], 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [EP_NEP_eastus_annual_mean_beta EP_NEP_eastus_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color','k')
set(gca, 'YLim',[-200 200]/1000, 'FontSize',8, 'YDir','reverse');
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

ax6 = subplot(4,2,4);
yyaxis left;
for i = 1:12
    lower = min(CP_NEP_eastus_monthly_beta(i, :));
    upper = max(CP_NEP_eastus_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [CP_NEP_eastus_monthly_mean_beta(i) CP_NEP_eastus_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
pos(2) = pos(2)-0.0;
set(gca, 'Position',pos, 'YDir','reverse');
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
% text(1, 3, 'F', 'FontSize',12);
yyaxis right;
lower = min(CP_NEP_eastus_annual_beta)/1000;
upper = max(CP_NEP_eastus_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [0.8 0.8 0.8], 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [CP_NEP_eastus_annual_mean_beta CP_NEP_eastus_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color','k')
set(gca, 'YLim',[-200 200]/1000, 'YDir','reverse');
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

% Western U.S.
load('./data/cp_ep_nep_mstmip_regional.mat');
ax7 = subplot(4,2,5);
yyaxis left;
for i = 1:12
    lower = min(EP_NEP_westna_monthly_beta(i, :));
    upper = max(EP_NEP_westna_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [EP_NEP_westna_monthly_mean_beta(i) EP_NEP_westna_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
pos(2) = pos(2)-0.02;
set(gca, 'Position',pos, 'YDir','reverse');
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
text(1, 0.75, 'Western North America', 'FontSize',10);
yyaxis right;
lower = min(EP_NEP_westna_annual_beta)/1000;
upper = max(EP_NEP_westna_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [0.8 0.8 0.8], 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [EP_NEP_westna_annual_mean_beta EP_NEP_westna_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color','k')
set(gca, 'YLim',[-200 200]/1000, 'FontSize',8, 'YDir','reverse');
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

ax8 = subplot(4,2,6);
yyaxis left;
for i = 1:12
    lower = min(CP_NEP_westna_monthly_beta(i, :));
    upper = max(CP_NEP_westna_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [CP_NEP_westna_monthly_mean_beta(i) CP_NEP_westna_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
pos(2) = pos(2)-0.02;
set(gca, 'Position',pos, 'YDir','reverse');
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
% text(1, 3, 'F', 'FontSize',12);
yyaxis right;
lower = min(CP_NEP_westna_annual_beta)/1000;
upper = max(CP_NEP_westna_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [0.8 0.8 0.8], 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [CP_NEP_westna_annual_mean_beta CP_NEP_westna_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color','k')
set(gca, 'YLim',[-200 200]/1000, 'YDir','reverse');
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

% Europe
load('./data/cp_ep_nep_mstmip_regional.mat');
ax3 = subplot(4,2,7);
yyaxis left;
for i = 1:12
    lower = min(EP_NEP_europe_monthly_beta(i, :));
    upper = max(EP_NEP_europe_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [EP_NEP_europe_monthly_mean_beta(i) EP_NEP_europe_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
pos(2) = pos(2)-0.04;
set(gca, 'Position',pos, 'YDir','reverse');
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
text(1, 0.75, 'Europe', 'FontSize',10);
yyaxis right;
lower = min(EP_NEP_europe_annual_beta)/1000;
upper = max(EP_NEP_europe_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [0.8 0.8 0.8], 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [EP_NEP_europe_annual_mean_beta EP_NEP_europe_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color','k')
set(gca, 'YLim',[-200 200]/1000, 'FontSize',8, 'YDir','reverse');
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

ax4 = subplot(4,2,8);
yyaxis left;
for i = 1:12
    lower = min(CP_NEP_europe_monthly_beta(i, :));
    upper = max(CP_NEP_europe_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [CP_NEP_europe_monthly_mean_beta(i) CP_NEP_europe_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
pos(2) = pos(2)-0.04;
set(gca, 'Position',pos, 'YDir','reverse');
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
yyaxis right;
lower = min(CP_NEP_europe_annual_beta)/1000;
upper = max(CP_NEP_europe_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [0.8 0.8 0.8], 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [CP_NEP_europe_annual_mean_beta CP_NEP_europe_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color','k')
set(gca, 'YLim',[-200 200]/1000, 'YDir','reverse');
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];


%% Plot LUE beta through time
load('./data/cp_ep_nep_inversions_regional.mat');

clr2 = [166,206,227
31,120,180
178,223,138
51,160,44]/255;

% Central Asia
set(h, 'currentaxes',ax1);
yyaxis left;
for i = 1:12
    for j = 1:4
        hold on;
        plot([i-0.4 i+0.4], [EP_NEP_casia_monthly_beta(i,j) EP_NEP_casia_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [EP_NEP_casia_annual_beta(1) EP_NEP_casia_annual_beta(1)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [EP_NEP_casia_annual_beta(2) EP_NEP_casia_annual_beta(2)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(2,:))
plot([i-0.4 i+0.4], [EP_NEP_casia_annual_beta(3) EP_NEP_casia_annual_beta(3)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(3,:))
plot([i-0.4 i+0.4], [EP_NEP_casia_annual_beta(4) EP_NEP_casia_annual_beta(4)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(4,:))

set(h, 'currentaxes',ax2); 
yyaxis left;
for i = 1:12
    for j = 1:4
        hold on;
        plot([i-0.4 i+0.4], [CP_NEP_casia_monthly_beta(i,j) CP_NEP_casia_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [CP_NEP_casia_annual_beta(1) CP_NEP_casia_annual_beta(1)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [CP_NEP_casia_annual_beta(2) CP_NEP_casia_annual_beta(2)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(2,:))
plot([i-0.4 i+0.4], [CP_NEP_casia_annual_beta(3) CP_NEP_casia_annual_beta(3)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(3,:))
plot([i-0.4 i+0.4], [CP_NEP_casia_annual_beta(4) CP_NEP_casia_annual_beta(4)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(4,:))



% Eastern U.S.
set(h, 'currentaxes',ax5);
yyaxis left;
for i = 1:12
    for j = 1:4
        hold on;
        plot([i-0.4 i+0.4], [EP_NEP_eastus_monthly_beta(i,j) EP_NEP_eastus_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [EP_NEP_eastus_annual_beta(1) EP_NEP_eastus_annual_beta(1)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [EP_NEP_eastus_annual_beta(2) EP_NEP_eastus_annual_beta(2)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(2,:))
plot([i-0.4 i+0.4], [EP_NEP_eastus_annual_beta(3) EP_NEP_eastus_annual_beta(3)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(3,:))
plot([i-0.4 i+0.4], [EP_NEP_eastus_annual_beta(4) EP_NEP_eastus_annual_beta(4)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(4,:))

set(h, 'currentaxes',ax6); 
yyaxis left;
for i = 1:12
    for j = 1:4
        hold on;
        plot([i-0.4 i+0.4], [CP_NEP_eastus_monthly_beta(i,j) CP_NEP_eastus_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [CP_NEP_eastus_annual_beta(1) CP_NEP_eastus_annual_beta(1)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [CP_NEP_eastus_annual_beta(2) CP_NEP_eastus_annual_beta(2)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(2,:))
plot([i-0.4 i+0.4], [CP_NEP_eastus_annual_beta(3) CP_NEP_eastus_annual_beta(3)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(3,:))
plot([i-0.4 i+0.4], [CP_NEP_eastus_annual_beta(4) CP_NEP_eastus_annual_beta(4)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(4,:))


% western North America
set(h, 'currentaxes',ax7);
yyaxis left;
for i = 1:12
    for j = 1:4
        hold on;
        plot([i-0.4 i+0.4], [EP_NEP_westna_monthly_beta(i,j) EP_NEP_westna_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [EP_NEP_westna_annual_beta(1) EP_NEP_westna_annual_beta(1)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [EP_NEP_westna_annual_beta(2) EP_NEP_westna_annual_beta(2)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(2,:))
plot([i-0.4 i+0.4], [EP_NEP_westna_annual_beta(3) EP_NEP_westna_annual_beta(3)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(3,:))
plot([i-0.4 i+0.4], [EP_NEP_westna_annual_beta(4) EP_NEP_westna_annual_beta(4)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(4,:))

set(h, 'currentaxes',ax8); 
yyaxis left;
for i = 1:12
    for j = 1:4
        hold on;
        plot([i-0.4 i+0.4], [CP_NEP_westna_monthly_beta(i,j) CP_NEP_westna_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [CP_NEP_westna_annual_beta(1) CP_NEP_westna_annual_beta(1)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [CP_NEP_westna_annual_beta(2) CP_NEP_westna_annual_beta(2)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(2,:))
plot([i-0.4 i+0.4], [CP_NEP_westna_annual_beta(3) CP_NEP_westna_annual_beta(3)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(3,:))
plot([i-0.4 i+0.4], [CP_NEP_westna_annual_beta(4) CP_NEP_westna_annual_beta(4)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(4,:))


% europe
set(h, 'currentaxes',ax3);
yyaxis left;
for i = 1:12
    for j = 1:4
        hold on;
        plot([i-0.4 i+0.4], [EP_NEP_europe_monthly_beta(i,j) EP_NEP_europe_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [EP_NEP_europe_annual_beta(1) EP_NEP_europe_annual_beta(1)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [EP_NEP_europe_annual_beta(2) EP_NEP_europe_annual_beta(2)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(2,:))
plot([i-0.4 i+0.4], [EP_NEP_europe_annual_beta(3) EP_NEP_europe_annual_beta(3)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(3,:))
plot([i-0.4 i+0.4], [EP_NEP_europe_annual_beta(4) EP_NEP_europe_annual_beta(4)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(4,:))

set(h, 'currentaxes',ax4); 
yyaxis left;
for i = 1:12
    for j = 1:4
        hold on;
        plot([i-0.4 i+0.4], [CP_NEP_europe_monthly_beta(i,j) CP_NEP_europe_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [CP_NEP_europe_annual_beta(1) CP_NEP_europe_annual_beta(1)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [CP_NEP_europe_annual_beta(2) CP_NEP_europe_annual_beta(2)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(2,:))
plot([i-0.4 i+0.4], [CP_NEP_europe_annual_beta(3) CP_NEP_europe_annual_beta(3)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(3,:))
plot([i-0.4 i+0.4], [CP_NEP_europe_annual_beta(4) CP_NEP_europe_annual_beta(4)]/1000,...
    '-', 'LineWidth',2, 'Color',clr2(4,:))




%% Manual legend
set(h, 'currentaxes',ax3);
i = 4;
plot([i-0.4 i+0.4], -[0.16 0.16],...
    '-', 'LineWidth',2, 'Color',clr2(1,:));
text(i+0.5, -0.16, 'CAMS v15r2', 'Color',clr2(1,:), 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')
plot([i-0.4 i+0.4], -[0.12 0.12],...
    '-', 'LineWidth',2, 'Color',clr2(2,:))
text(i+0.5, -0.12, 'CAMS v15r4', 'Color',clr2(2,:), 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')
i = 9.5;
plot([i-0.4 i+0.4], -[0.16 0.16],...
    '-', 'LineWidth',2, 'Color',clr2(3,:));
text(i+0.5, -0.16, 'MACC v13r1', 'Color',clr2(3,:), 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')
plot([i-0.4 i+0.4], -[0.12 0.12],...
    '-', 'LineWidth',2, 'Color',clr2(4,:))
text(i+0.5, -0.12, 'MACC v14r2', 'Color',clr2(4,:), 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')

set(h, 'currentaxes',ax4);
i = 8.5;
fill([i-0.3 i+0.3 i+0.3 i-0.3], -1*[0.08 0.08 0.18 0.18],...
    [0.8 0.8 0.8], 'EdgeColor','none');
plot([i-0.3 i+0.3], -1*[0.13 0.13],...
    'k-', 'LineWidth',3)
text(i+0.5, -1*0.13, 'MsTMIP mean', 'Color','k', 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')
text(i+0.5, -1*0.17, 'MsTMIP range', 'Color',[0.7 0.7 0.7], 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/epi-cpi-nep-temperate.tif')
close all;

