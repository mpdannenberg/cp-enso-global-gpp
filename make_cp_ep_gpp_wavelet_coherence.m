% Wavelet coherence and cross spectrum
syear = 1982;
eyear = 2016;

%% Load monthly CP/EP indices
load ./data/cpi_epi_monthly.mat;
enso_yr = yr;
enso_mo = mo;

%% Figure setup
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 7];

cb = cbrewer('seq','Blues',10);

%% Load CCW GPP data
load ./data/gpp_ccw.mat;
clear GPP_annual GPP_global_annual GPP_monthly temp lat lon;
gpp_yr = yr; clear yr years;
gpp_mo = mo; clear mo;

GPP_global_monthly_anom = reshape([GPP_global_monthly - repmat(mean(GPP_global_monthly), length(syear:eyear), 1)]', length(gpp_yr), []);

%% Wavelet coherence of CCW GPP and CP/EP indices
subplot(3,2,1)
[epwcoh, epwcs, f, coi]=wcoherence(ep_idx(enso_yr>=syear & enso_yr<=eyear),...
    GPP_global_monthly_anom(gpp_yr>=syear & gpp_yr<=eyear),...
    years(1/12),'PhaseDisplayThreshold',0.7);
wcoherence(ep_idx(enso_yr>=syear & enso_yr<=eyear),...
    GPP_global_monthly_anom(gpp_yr>=syear & gpp_yr<=eyear),...
    years(1/12),'PhaseDisplayThreshold',1);
colormap(cb);
colorbar('off')
ttl = title('Eastern Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 5.5;
xlabel('')
ax = gca;
ax.Position(1) = 0.19;
set(ax, 'XTick',3:10:40, 'XTickLabel',{'1985','1995','2005','2015'},...
    'TickDir','out', 'TickLength',[0.02 0.03])
text(1, 4, 'A', 'FontSize',12)
text(-14, 1, 'CCW', 'FontWeight','bold','FontSize',12,'HorizontalAlignment','center','Rotation',90)

subplot(3,2,2)
[cpwcoh, cpwcs]=wcoherence(cp_idx(enso_yr>=syear & enso_yr<=eyear),...
    GPP_global_monthly_anom(gpp_yr>=syear & gpp_yr<=eyear),...
    years(1/12),'PhaseDisplayThreshold',0.7);
wcoherence(cp_idx(enso_yr>=syear & enso_yr<=eyear),...
    GPP_global_monthly_anom(gpp_yr>=syear & gpp_yr<=eyear),...
    years(1/12),'PhaseDisplayThreshold',1);
colormap(cb);
colorbar('off')
ttl = title('Central Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 5.5;
xlabel('')
ylabel('')
ax = gca;
ax.Position(1) = 0.62;
set(ax, 'XTick',3:10:40, 'XTickLabel',{'1985','1995','2005','2015'},...
    'TickDir','out', 'TickLength',[0.02 0.03])
text(1, 4, 'B', 'FontSize',12)

%% Load MOD17 GPP data
load ./data/gpp_mod17.mat;
clear GPP_annual GPP_global_annual GPP_monthly temp lat lon;
gpp_yr = yr; clear yr years;
gpp_mo = mo; clear mo;

GPP_global_monthly_anom = reshape([GPP_global_monthly - repmat(mean(GPP_global_monthly), length(syear:eyear), 1)]', length(gpp_yr), []);

%% Wavelet coherence of MOD17 GPP and CP/EP indices
subplot(3,2,3)
[epwcoh, epwcs, f, coi]=wcoherence(ep_idx(enso_yr>=syear & enso_yr<=eyear),...
    GPP_global_monthly_anom(gpp_yr>=syear & gpp_yr<=eyear),...
    years(1/12),'PhaseDisplayThreshold',0.7);
wcoherence(ep_idx(enso_yr>=syear & enso_yr<=eyear),...
    GPP_global_monthly_anom(gpp_yr>=syear & gpp_yr<=eyear),...
    years(1/12),'PhaseDisplayThreshold',1);
colormap(cb);
colorbar('off')
title('');
xlabel('')
ax = gca;
ax.Position(1) = 0.19;
ax.Position(2) = 0.43;
set(ax, 'XTick',3:10:40, 'XTickLabel',{'1985','1995','2005','2015'},...
    'TickDir','out', 'TickLength',[0.02 0.03])
text(1, 4, 'C', 'FontSize',12)
text(-14, 1, 'MOD17', 'FontWeight','bold','FontSize',12,'HorizontalAlignment','center','Rotation',90)

subplot(3,2,4)
[cpwcoh, cpwcs]=wcoherence(cp_idx(enso_yr>=syear & enso_yr<=eyear),...
    GPP_global_monthly_anom(gpp_yr>=syear & gpp_yr<=eyear),...
    years(1/12),'PhaseDisplayThreshold',0.7);
wcoherence(cp_idx(enso_yr>=syear & enso_yr<=eyear),...
    GPP_global_monthly_anom(gpp_yr>=syear & gpp_yr<=eyear),...
    years(1/12),'PhaseDisplayThreshold',1);
colormap(cb);
colorbar('off')
title('');
xlabel('')
ylabel('')
ax = gca;
ax.Position(1) = 0.62;
ax.Position(2) = 0.43;
set(ax, 'XTick',3:10:40, 'XTickLabel',{'1985','1995','2005','2015'},...
    'TickDir','out', 'TickLength',[0.02 0.03])
text(1, 4, 'D', 'FontSize',12)

%% Load MsTMIP GPP data
syear = 1951;
eyear = 2010;

load ./data/gpp_mstmip.mat;
clear GPP_annual_mean GPP_global_annual* GPP_monthly GPP_global_monthly temp lat lon;
gpp_yr = reshape(repmat(syear:eyear, 12, 1), 1, [])'; clear years;
gpp_mo = repmat(1:12, 1, length(syear:eyear))';

GPP_global_monthly_anom = reshape([GPP_global_monthly_mean - repmat(mean(GPP_global_monthly_mean), length(syear:eyear), 1)]', length(gpp_yr), []);

%% Wavelet coherence of MsTMIP GPP and CP/EP indices
subplot(3,2,5)
[epwcoh, epwcs, f, coi]=wcoherence(ep_idx(enso_yr>=syear & enso_yr<=eyear),...
    GPP_global_monthly_anom(gpp_yr>=syear & gpp_yr<=eyear),...
    years(1/12),'PhaseDisplayThreshold',0.7);
wcoherence(ep_idx(enso_yr>=syear & enso_yr<=eyear),...
    GPP_global_monthly_anom(gpp_yr>=syear & gpp_yr<=eyear),...
    years(1/12),'PhaseDisplayThreshold',1);
colormap(cb);
colorbar('off')
title('');
xlabel('')
ax = gca;
ax.Position(1) = 0.19;
ax.Position(2) = 0.15;
ax.Position(3) = 0.3347;
ax.Position(4) = 0.2157;
set(ax, 'XTick',4:10:60, 'XTickLabel',{'1955','1965','1975','1985','1995','2005'},...
    'TickDir','out', 'TickLength',[0.02 0.03])
text(2, 5, 'E', 'FontSize',12)
text(-24, 1.5, 'MsTMIP', 'FontWeight','bold','FontSize',12,'HorizontalAlignment','center','Rotation',90)

subplot(3,2,6)
[cpwcoh, cpwcs]=wcoherence(cp_idx(enso_yr>=syear & enso_yr<=eyear),...
    GPP_global_monthly_anom(gpp_yr>=syear & gpp_yr<=eyear),...
    years(1/12),'PhaseDisplayThreshold',0.7);
wcoherence(cp_idx(enso_yr>=syear & enso_yr<=eyear),...
    GPP_global_monthly_anom(gpp_yr>=syear & gpp_yr<=eyear),...
    years(1/12),'PhaseDisplayThreshold',1);
colormap(cb);
colorbar('off')
title('');
xlabel('')
ylabel('')
ax = gca;
ax.Position(1) = 0.62;
ax.Position(2) = 0.15;
ax.Position(3) = 0.3347;
ax.Position(4) = 0.2157;
set(ax, 'XTick',4:10:60, 'XTickLabel',{'1955','1965','1975','1985','1995','2005'},...
    'TickDir','out', 'TickLength',[0.02 0.03])
text(2, 5, 'F', 'FontSize',12)
cb = colorbar('southoutside');
cb.Position = [0.19    0.08    0.765    0.02];
cb.TickLength = 0.027;
xlabel(cb, 'Magnitude-squared coherence', 'FontSize',11);

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/ep-cp-gpp-wavelet-coherence.tif')
close all;

