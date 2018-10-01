% Make some plots

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 0 7 8];

set(h, 'defaultAxesColorOrder',[0 0 0; 140,81,10]/255)


%% Plot MsTMIP beta through time
% Central Asia
load('./data/cp_ep_gpp_mstmip_regional.mat');
ax1 = subplot(4,2,1);
yyaxis left;
for i = 1:12
    lower = min(EP_GPP_casia_monthly_beta(i, :));
    upper = max(EP_GPP_casia_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [EP_GPP_casia_monthly_mean_beta(i) EP_GPP_casia_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
pos(2) = pos(2)+0.02;
set(gca, 'Position',pos);
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
text(1, 0.75, 'Central Asia', 'FontSize',10);
yyaxis right;
lower = min(EP_GPP_casia_annual_beta)/1000;
upper = max(EP_GPP_casia_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [223,194,125]/255, 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [EP_GPP_casia_annual_mean_beta EP_GPP_casia_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color',[140,81,10]/255)
set(gca, 'YLim',[-200 200]/1000, 'FontSize',8);
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];
ttl = title('Eastern Pacific ENSO', 'FontSize',11);
ttl.Position(2) = 1.25;

ax2 = subplot(4,2,2);
yyaxis left;
for i = 1:12
    lower = min(CP_GPP_casia_monthly_beta(i, :));
    upper = max(CP_GPP_casia_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [CP_GPP_casia_monthly_mean_beta(i) CP_GPP_casia_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
pos(2) = pos(2)+0.02;
set(gca, 'Position',pos);
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
% text(1, 3, 'F', 'FontSize',12);
yyaxis right;
lower = min(CP_GPP_casia_annual_beta)/1000;
upper = max(CP_GPP_casia_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [223,194,125]/255, 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [CP_GPP_casia_annual_mean_beta CP_GPP_casia_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color',[140,81,10]/255)
set(gca, 'YLim',[-200 200]/1000);
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];
ttl = title('Central Pacific ENSO', 'FontSize',11);
ttl.Position(2) = 1.25;

% Eastern U.S.
load('./data/cp_ep_gpp_mstmip_regional.mat');
ax5 = subplot(4,2,3);
yyaxis left;
for i = 1:12
    lower = min(EP_GPP_eastus_monthly_beta(i, :));
    upper = max(EP_GPP_eastus_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [EP_GPP_eastus_monthly_mean_beta(i) EP_GPP_eastus_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
pos(2) = pos(2)-0.0;
set(gca, 'Position',pos);
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
text(1, 0.75, 'Eastern U.S.', 'FontSize',10);
yyaxis right;
lower = min(EP_GPP_eastus_annual_beta)/1000;
upper = max(EP_GPP_eastus_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [223,194,125]/255, 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [EP_GPP_eastus_annual_mean_beta EP_GPP_eastus_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color',[140,81,10]/255)
set(gca, 'YLim',[-200 200]/1000, 'FontSize',8);
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

ax6 = subplot(4,2,4);
yyaxis left;
for i = 1:12
    lower = min(CP_GPP_eastus_monthly_beta(i, :));
    upper = max(CP_GPP_eastus_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [CP_GPP_eastus_monthly_mean_beta(i) CP_GPP_eastus_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
pos(2) = pos(2)-0.0;
set(gca, 'Position',pos);
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
% text(1, 3, 'F', 'FontSize',12);
yyaxis right;
lower = min(CP_GPP_eastus_annual_beta)/1000;
upper = max(CP_GPP_eastus_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [223,194,125]/255, 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [CP_GPP_eastus_annual_mean_beta CP_GPP_eastus_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color',[140,81,10]/255)
set(gca, 'YLim',[-200 200]/1000);
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

% Western U.S.
load('./data/cp_ep_gpp_mstmip_regional.mat');
ax7 = subplot(4,2,5);
yyaxis left;
for i = 1:12
    lower = min(EP_GPP_westna_monthly_beta(i, :));
    upper = max(EP_GPP_westna_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [EP_GPP_westna_monthly_mean_beta(i) EP_GPP_westna_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
pos(2) = pos(2)-0.02;
set(gca, 'Position',pos);
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
text(1, 0.75, 'Western North America', 'FontSize',10);
yyaxis right;
lower = min(EP_GPP_westna_annual_beta)/1000;
upper = max(EP_GPP_westna_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [223,194,125]/255, 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [EP_GPP_westna_annual_mean_beta EP_GPP_westna_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color',[140,81,10]/255)
set(gca, 'YLim',[-200 200]/1000, 'FontSize',8);
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

ax8 = subplot(4,2,6);
yyaxis left;
for i = 1:12
    lower = min(CP_GPP_westna_monthly_beta(i, :));
    upper = max(CP_GPP_westna_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [CP_GPP_westna_monthly_mean_beta(i) CP_GPP_westna_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
pos(2) = pos(2)-0.02;
set(gca, 'Position',pos);
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
% text(1, 3, 'F', 'FontSize',12);
yyaxis right;
lower = min(CP_GPP_westna_annual_beta)/1000;
upper = max(CP_GPP_westna_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [223,194,125]/255, 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [CP_GPP_westna_annual_mean_beta CP_GPP_westna_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color',[140,81,10]/255)
set(gca, 'YLim',[-200 200]/1000);
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

% Europe
load('./data/cp_ep_gpp_mstmip_regional.mat');
ax3 = subplot(4,2,7);
yyaxis left;
for i = 1:12
    lower = min(EP_GPP_europe_monthly_beta(i, :));
    upper = max(EP_GPP_europe_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [EP_GPP_europe_monthly_mean_beta(i) EP_GPP_europe_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
pos(2) = pos(2)-0.04;
set(gca, 'Position',pos);
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
text(1, 0.75, 'Europe', 'FontSize',10);
yyaxis right;
lower = min(EP_GPP_europe_annual_beta)/1000;
upper = max(EP_GPP_europe_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [223,194,125]/255, 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [EP_GPP_europe_annual_mean_beta EP_GPP_europe_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color',[140,81,10]/255)
set(gca, 'YLim',[-200 200]/1000, 'FontSize',8);
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];

ax4 = subplot(4,2,8);
yyaxis left;
for i = 1:12
    lower = min(CP_GPP_europe_monthly_beta(i, :));
    upper = max(CP_GPP_europe_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [CP_GPP_europe_monthly_mean_beta(i) CP_GPP_europe_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end
set(gca, 'XLim',[0 15], 'YLim',[-1 1], 'XTick',[1:12 14], 'TickDir','out', 'FontSize',8,...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Annual'});
pos = get(gca, 'Position');
pos(1) = pos(1)+0.02;
pos(2) = pos(2)-0.04;
set(gca, 'Position',pos);
hold off;
ylabel('Tg C day^{-1}', 'FontSize',8);
yyaxis right;
lower = min(CP_GPP_europe_annual_beta)/1000;
upper = max(CP_GPP_europe_annual_beta)/1000;
i = 14;
fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
    [223,194,125]/255, 'EdgeColor','none');
hold on;
plot([i-0.4 i+0.4], [CP_GPP_europe_annual_mean_beta CP_GPP_europe_annual_mean_beta]/1000,...
    '-', 'LineWidth',3, 'Color',[140,81,10]/255)
set(gca, 'YLim',[-200 200]/1000);
plot([0 15],[0 0],'k-')
ylb = ylabel('Pg C yr^{-1}', 'FontSize',8);
ylb.Position = [16.8    0.0000   -1.0000];


%% Plot LUE beta through time
load('./data/cp_ep_gpp_lue_regional.mat');

clr2 = [55,126,184
77,175,74]/255;

% Central Asia
set(h, 'currentaxes',ax1);
yyaxis left;
for i = 1:12
    for j = 1:2
        hold on;
        plot([i-0.4 i+0.4], [EP_GPP_casia_monthly_beta(i,j) EP_GPP_casia_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [EP_GPP_casia_annual_beta(1) EP_GPP_casia_annual_beta(1)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [EP_GPP_casia_annual_beta(2) EP_GPP_casia_annual_beta(2)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(2,:))

set(h, 'currentaxes',ax2); 
yyaxis left;
for i = 1:12
    for j = 1:2
        hold on;
        plot([i-0.4 i+0.4], [CP_GPP_casia_monthly_beta(i,j) CP_GPP_casia_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [CP_GPP_casia_annual_beta(1) CP_GPP_casia_annual_beta(1)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [CP_GPP_casia_annual_beta(2) CP_GPP_casia_annual_beta(2)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(2,:))



% Eastern U.S.
set(h, 'currentaxes',ax5);
yyaxis left;
for i = 1:12
    for j = 1:2
        hold on;
        plot([i-0.4 i+0.4], [EP_GPP_eastus_monthly_beta(i,j) EP_GPP_eastus_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [EP_GPP_eastus_annual_beta(1) EP_GPP_eastus_annual_beta(1)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [EP_GPP_eastus_annual_beta(2) EP_GPP_eastus_annual_beta(2)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(2,:))

set(h, 'currentaxes',ax6); 
yyaxis left;
for i = 1:12
    for j = 1:2
        hold on;
        plot([i-0.4 i+0.4], [CP_GPP_eastus_monthly_beta(i,j) CP_GPP_eastus_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [CP_GPP_eastus_annual_beta(1) CP_GPP_eastus_annual_beta(1)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [CP_GPP_eastus_annual_beta(2) CP_GPP_eastus_annual_beta(2)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(2,:))


% western North America
set(h, 'currentaxes',ax7);
yyaxis left;
for i = 1:12
    for j = 1:2
        hold on;
        plot([i-0.4 i+0.4], [EP_GPP_westna_monthly_beta(i,j) EP_GPP_westna_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [EP_GPP_westna_annual_beta(1) EP_GPP_westna_annual_beta(1)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [EP_GPP_westna_annual_beta(2) EP_GPP_westna_annual_beta(2)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(2,:))

set(h, 'currentaxes',ax8); 
yyaxis left;
for i = 1:12
    for j = 1:2
        hold on;
        plot([i-0.4 i+0.4], [CP_GPP_westna_monthly_beta(i,j) CP_GPP_westna_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [CP_GPP_westna_annual_beta(1) CP_GPP_westna_annual_beta(1)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [CP_GPP_westna_annual_beta(2) CP_GPP_westna_annual_beta(2)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(2,:))


% europe
set(h, 'currentaxes',ax3);
yyaxis left;
for i = 1:12
    for j = 1:2
        hold on;
        plot([i-0.4 i+0.4], [EP_GPP_europe_monthly_beta(i,j) EP_GPP_europe_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [EP_GPP_europe_annual_beta(1) EP_GPP_europe_annual_beta(1)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [EP_GPP_europe_annual_beta(2) EP_GPP_europe_annual_beta(2)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(2,:))

set(h, 'currentaxes',ax4); 
yyaxis left;
for i = 1:12
    for j = 1:2
        hold on;
        plot([i-0.4 i+0.4], [CP_GPP_europe_monthly_beta(i,j) CP_GPP_europe_monthly_beta(i,j)],...
            'k-', 'LineWidth',3, 'Color',clr2(j,:))
    end
end
yyaxis right;
i = 14;
plot([i-0.4 i+0.4], [CP_GPP_europe_annual_beta(1) CP_GPP_europe_annual_beta(1)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(1,:))
plot([i-0.4 i+0.4], [CP_GPP_europe_annual_beta(2) CP_GPP_europe_annual_beta(2)]/1000,...
    '-', 'LineWidth',3, 'Color',clr2(2,:))




%% Manual legend
i = 4.5;
plot([i-0.4 i+0.4], [0.16 0.16],...
    '-', 'LineWidth',3, 'Color',clr2(1,:));
text(i+0.5, 0.16, 'CCW', 'Color',clr2(1,:), 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')
plot([i-0.4 i+0.4], [0.12 0.12],...
    '-', 'LineWidth',3, 'Color',clr2(2,:))
text(i+0.5, 0.12, 'MOD17', 'Color',clr2(2,:), 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')

i = 8.5;
fill([i-0.3 i+0.3 i+0.3 i-0.3], [0.08 0.08 0.18 0.18],...
    [0.8 0.8 0.8], 'EdgeColor','none');
plot([i-0.3 i+0.3], [0.13 0.13],...
    'k-', 'LineWidth',3)
text(i+0.5, 0.13, 'MsTMIP mean', 'Color','k', 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')
text(i+0.5, 0.17, 'MsTMIP range', 'Color',[0.7 0.7 0.7], 'FontWeight','bold', 'FontSize',8,...
    'VerticalAlignment','middle', 'HorizontalAlignment','left')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/epi-cpi-gpp-temperate.tif')
close all;

