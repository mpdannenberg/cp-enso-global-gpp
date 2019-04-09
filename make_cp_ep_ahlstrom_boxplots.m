% Make plots of Ahlstrom regional GPP/NEE response to EP/CP ENSO

% Parameters
wdth = 0.25;
offs = 0.15;
scale = 0.001;

load ./data/cp_ep_gpp_mstmip_regional.mat;

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 4];

plot([0 6], [0 0], 'k-')
hold on;
dberg_box(1-offs,CP_GPP_tropical_annual_beta*scale,'r','w',wdth);
dberg_box(2-offs,CP_GPP_semiarid_annual_beta*scale,'r','w',wdth);
dberg_box(4-offs,CP_GPP_grass_annual_beta*scale,'r','w',wdth);
dberg_box(3-offs,CP_GPP_extratropical_annual_beta*scale,'r','w',wdth);
dberg_box(5-offs,CP_GPP_tundra_annual_beta*scale,'r','w',wdth);

dberg_box(1+offs,EP_GPP_tropical_annual_beta*scale,'k','w',wdth);
hold on;
dberg_box(2+offs,EP_GPP_semiarid_annual_beta*scale,'k','w',wdth);
dberg_box(4+offs,EP_GPP_grass_annual_beta*scale,'k','w',wdth);
dberg_box(3+offs,EP_GPP_extratropical_annual_beta*scale,'k','w',wdth);
dberg_box(5+offs,EP_GPP_tundra_annual_beta*scale,'k','w',wdth);

% Means and CIs
plot(1-offs/2,CP_GPP_tropical_annual_mean_beta*scale,'^r');
plot(2-offs/2,CP_GPP_semiarid_annual_mean_beta*scale,'^r');
plot(4-offs/2,CP_GPP_grass_annual_mean_beta*scale,'^r');
plot(3-offs/2,CP_GPP_extratropical_annual_mean_beta*scale,'^r');
plot(5-offs/2,CP_GPP_tundra_annual_mean_beta*scale,'^r');

plot(1+offs/2,EP_GPP_tropical_annual_mean_beta*scale,'^k');
plot(2+offs/2,EP_GPP_semiarid_annual_mean_beta*scale,'^k');
plot(4+offs/2,EP_GPP_grass_annual_mean_beta*scale,'^k');
plot(3+offs/2,EP_GPP_extratropical_annual_mean_beta*scale,'^k');
plot(5+offs/2,EP_GPP_tundra_annual_mean_beta*scale,'^k');


set(gca, 'XLim',[0.5 5.5], 'XTick',1:5, 'TickDir','out',...
    'TickLength',[0.02 0.04], 'XTickLabels',{'Tropical Forest','Semiarid',...
    'Extratropical Forest','Grass/Crop','Tundra/Arctic Shrub'})

ylabel('ENSO GPP response (Pg C yr^{-1} SD^{-1})')




