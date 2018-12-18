% Correlations between tree-ring widths and EPI/CPI

syear = 1959;
eyear = 2016;

%% load and process Mauna Loa data
% CO2 growth rate anomaly
T = readtable('./data/co2_gr_mlo.txt', 'ReadVariableNames',true, 'HeaderLines',59);
mdl = fitlm(T.year, T.anninc);
gr = mdl.Residuals.Raw(T.year>=syear & T.year<=eyear);

% CO2 seasonal amplitude
T = readtable('./data/co2_mm_mlo.txt', 'ReadVariableNames',true, 'HeaderLines',72);
am = NaN(size(gr));
for i = syear:eyear
    co2 = T.Interpolated(T.year == i);
    am(i-syear+1) = max(co2)-min(co2);
end

clear T co2 mdl i;

%% load EPI and CPI
load ./data/cpi_epi_1951-2016;
[~,idx,~] = intersect(yr, syear:eyear);
epi = epi(idx); cpi = cpi(idx); yr = yr(idx); clear idx;

%% Plot time series
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 4];

ha = tight_subplot(3, 1, 0, [0.1 0.1], [0.1 0.1]);

axes(ha(1))
plot(yr, gr, 'k-', 'LineWidth',1.5);
set(gca, 'XLim',[min(yr) max(yr)], 'XTickLabels','','XColor','w', 'TickDir','out','TickLength',[0.01 0.02]);
ylabel({'CO_{2} growth rate', 'anomaly (ppm)'},'FontSize',8)
text(1961, 1.5, 'A', 'FontSize',12)
box off;

axes(ha(2))
plot(yr, epi, 'k-', 'LineWidth',1.5);
set(gca, 'XLim',[min(yr) max(yr)], 'YAxisLocation','right', 'XTickLabels','', 'YLim',[-3 3],'XColor','w', 'TickDir','out','TickLength',[0.01 0.02]);
ylabel({'Eastern Pacific', 'ENSO index'},'FontSize',8)
text(1961, 2.25, 'B', 'FontSize',12)
box off;

axes(ha(3))
plot(yr, cpi, 'k-', 'LineWidth',1.5);
set(gca, 'XLim',[min(yr) max(yr)], 'YLim',[-3 3], 'TickDir','out','TickLength',[0.01 0.02]);
ylabel({'Central Pacific', 'ENSO index'},'FontSize',8)
text(1961, 2.25, 'C', 'FontSize',12)
box off;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-mauna-loa-time-series.tif')
close all;

%% Scatterplots
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 3];

ha = tight_subplot(1, 2, 0, [0.15 0.05], [0.1 0.1]);

% Growth rate anomalies
axes(ha(1))
plot(epi, gr, 'ko', 'LineWidth',1.5);
set(gca, 'XLim',[-3.5 3.5], 'YLim',[-1.75 1.75],'TickLength',[0 0]);
ylabel({'CO_{2} growth rate anomaly (ppm)'},'FontSize',9)
grid on;
[r,p] = corr(epi, gr);
text(-3,-1.1,['R = ',num2str(round(r,2))],'FontSize',10)
text(-3,-1.4,['{\itp} = ',num2str(round(p,2))],'FontSize',10)
mdl = fitlm(epi, gr);
hold on;
plot([min(epi) max(epi)],...
    mdl.Coefficients.Estimate(2)*[min(epi) max(epi)]+mdl.Coefficients.Estimate(1),...
    '-', 'Color',[0.4 0.4 0.4], 'LineWidth',2)
hold off;
text(-2.9, 1.4, 'A', 'FontSize',14);
xlabel('Eastern Pacific ENSO index','FontSize',9)

axes(ha(2))
plot(cpi, gr, 'ko', 'LineWidth',1.5);
set(gca, 'XLim',[-3.5 3.5], 'YLim',[-1.75 1.75],'TickLength',[0 0],'YAxisLocation','right');
ylabel({'CO_{2} growth rate anomaly (ppm)'},'FontSize',9)
grid on;
[r,p] = corr(cpi, gr);
text(2,1.4,['R = ',num2str(round(r,2))],'FontSize',10)
text(2,1.1,['{\itp} = ',num2str(round(p,2))],'FontSize',10)
text(-2.9, 1.4, 'B', 'FontSize',14);
xlabel('Central Pacific ENSO index','FontSize',9)

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-mauna-loa-scatterplot.tif')
close all;


