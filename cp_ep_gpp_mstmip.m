% Examine relationships between modeled MsTMIP GPP and EPI/CPI over 
% 1951-2010 period

%% Units for final GPP data (original: kgC m-2 s-1)
% monthly gridded GPP: kgC m-2 day-1
% monthly global GPP: TgC day-1
% annual gridded GPP: kgC m-2 year-1
% annual global GPP: TgC year-1

%% Setup
syear = 1951; % First year of analysis
scale = 10^-9; % kg --> Tg
models = {'BIOME-BGC','CLASS-CTEM-N','CLM4','CLM4VIC','DLEM','GTEC',...
    'ISAM','LPJ-wsl','ORCHIDEE-LSCE','SiB3','SiBCASA','TEM6','VEGAS2.1',...
    'VISIT'};

%% Load MsTMIP GPP data
cd('C:\Users\dannenberg\Documents\Data_Analysis\MsTMIP');
lat = ncread('BIOME-BGC_SG1_Monthly_GPP.nc4','lat');
lon = ncread('BIOME-BGC_SG1_Monthly_GPP.nc4','lon');
ny = length(lat); nx = length(lon);

yr = reshape(repmat(1901:2010, 12, 1), [], 1);
idx = yr >= syear;

yr = yr(idx); nt = length(yr);
mo = repmat(1:12, 1, length(syear:2010))';
ndys = repmat([31,28,31,30,31,30,31,31,30,31,30,31], 1, length(syear:2010))'; % Number of days, excluding leap years

GPP = NaN(ny, nx, nt, length(models));
for i = 1:length(models)
    
    gpp = ncread([models{i},'_SG1_Monthly_GPP.nc4'],'GPP') * 60 * 60 * 24; % from kgC m-2 s-1 --> kgC m-2 day-1
    gpp(gpp==-9999) = NaN;
    
    GPP(:, :, :, i) = permute(gpp(:, :, idx), [2 1 3]);
    
end
clear i gpp idx;

%% Aggregate to monthly and annual scales
% Multimodel monthly gridded mean
GPP_monthly = bsxfun(@times,GPP,reshape(ndys,1,1,[])); % daily average --> monthly total GPP

% Multimodel annual gridded mean
windowSize = 12;
b = ones(1,windowSize);
a = 1;
GPP_annual = filter(b, a, GPP_monthly, [], 3); % 12-month running sums (kgC m-2 yr-1)
GPP_annual = GPP_annual(:, :, mo==12, :); % Get calendar year sum
GPP_annual_mean = nanmean(GPP_annual, 4);
clear GPP_monthly a b windowSize ndys;

%% Calculate global GPP at monthly and annual scale
[LON, LAT] = meshgrid(lon, lat);
e = referenceEllipsoid('World Geodetic System 1984');
area = areaquad(reshape(LAT-0.25,[],1),reshape(LON-0.25,[],1),reshape(LAT+0.25,[],1),reshape(LON+0.25,[],1),e);
area = reshape(area, ny, nx); 
clear LON LAT e;

yrs = syear:2010;
GPP_global_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP(:, :, yr==yrs(i) & mo==j, k);
            GPP_global_monthly(i,j,k) = nansum(nansum( gpp.*area )) * scale; % TgC day-1
        end
    end
end
GPP_global_monthly_mean = nanmean(GPP_global_monthly, 3);

GPP_global_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_global_annual(i, k) = nansum(nansum( gpp.*area )) * scale; % TgC yr-1
    end
end
GPP_global_annual_mean = nanmean(GPP_global_annual, 2);

clear i j k gpp yrs area;


%% Regress monthly & annual gpp against CPI and EPI indices
cd('C:\Users\dannenberg\Documents\Publications\Dannenberg_et_al_CPElNinoGlobalGPP\cp-enso-global-gpp');
load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr<=2010);
epi = epi(yr<=2010);
yr = yr(yr<=2010);

EP_GPP_global_monthly_beta = NaN(12, length(models));
EP_GPP_global_monthly_mean_beta = NaN(12, 1);
EP_GPP_global_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_global_annual_mean);
EP_GPP_global_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(epi, GPP_global_monthly_mean(:, i));
    EP_GPP_global_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_global_monthly(:, i, j));
        EP_GPP_global_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_global_annual(:, j));
            EP_GPP_global_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_GPP_global_monthly_beta = NaN(12, length(models));
CP_GPP_global_monthly_mean_beta = NaN(12, 1);
CP_GPP_global_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_global_annual_mean);
CP_GPP_global_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(cpi, GPP_global_monthly_mean(:, i));
    CP_GPP_global_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_global_monthly(:, i, j));
        CP_GPP_global_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_global_annual(:, j));
            CP_GPP_global_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

%% Regress EPI and CPI vs. gridded GPP
EP_GPP_annual_r = NaN(ny, nx);
EP_GPP_annual_p = NaN(ny, nx);
EP_GPP_annual_beta = NaN(ny, nx);
CP_GPP_annual_r = NaN(ny, nx);
CP_GPP_annual_p = NaN(ny, nx);
CP_GPP_annual_beta = NaN(ny, nx);

for i = 1:ny
    for j = 1:nx
        ts = squeeze(GPP_annual_mean(i,j,:));
        if sum(isnan(ts))==0
            [r,p] = corr(epi, ts);
            mdl = fitlm(epi, ts);
            EP_GPP_annual_r(i,j) = r;
            EP_GPP_annual_p(i,j) = p;
            EP_GPP_annual_beta(i,j) = mdl.Coefficients.Estimate(2);
            
            [r,p] = corr(cpi, ts);
            mdl = fitlm(cpi, ts);
            CP_GPP_annual_r(i,j) = r;
            CP_GPP_annual_p(i,j) = p;
            CP_GPP_annual_beta(i,j) = mdl.Coefficients.Estimate(2);
        end
    end
end

clear i j ts r p mdl nt nx ny scale syear GPP mo;

save('./data/cp_ep_gpp_mstmip.mat');

%% Make some maps and stuff
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
pos(1) = pos(1)-0.08;
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
pos(1) = pos(1)-0.1;
set(gca, 'Position',pos);
ttl = title('Central Pacific Index','FontSize',14);
ttl.Position = [0.0000    1.9    0.0000];
text(-2.2,1.3,'B', 'FontSize',12);

cb = colorbar('eastoutside');
cb.Position = [0.84    0.74    0.04    0.22];
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
pos(1) = pos(1)-0.08;
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
pos(1) = pos(1)-0.1;
set(gca, 'Position',pos);
text(-2.2,1.3,'D', 'FontSize',12);

cb = colorbar('eastoutside');
cb.Position = [0.84    0.41    0.04    0.22];
cb.Ticks = -0.05:0.01:0.05;
cb.TickLength = 0.21;
cb.TickLabels = {'-0.05','','','','','0','','','','','0.05'};
ylabel(cb, '\beta (kg C m^{-2} yr^{-1})');

% Plot beta through time
subplot(3,2,5)
for i = 1:12
    lower = min(EP_GPP_global_monthly_beta(i, :));
    upper = max(EP_GPP_global_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [EP_GPP_global_monthly_mean_beta(i) EP_GPP_global_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end

set(gca, 'XLim',[0 15], 'XTick',[1:12 14], 'TickDir','out',...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Year'});
pos = get(gca, 'Position');
pos(1) = pos(1)-0.08;
set(gca, 'Position',pos);

subplot(3,2,6)
for i = 1:12
    lower = min(CP_GPP_global_monthly_beta(i, :));
    upper = max(CP_GPP_global_monthly_beta(i, :));
    fill([i-0.4 i+0.4 i+0.4 i-0.4], [lower lower upper upper],...
        [0.8 0.8 0.8], 'EdgeColor','none');
    hold on;
    plot([i-0.4 i+0.4], [CP_GPP_global_monthly_mean_beta(i) CP_GPP_global_monthly_mean_beta(i)],...
        'k-', 'LineWidth',3)
end

set(gca, 'XLim',[0 15], 'YLim',[-4 4], 'XTick',[1:12 14], 'TickDir','out',...
    'TickLength',[0.025 0.05], 'XTickLabels',{'J','F','M','A','M','J','J','A','S','O','N','D','Year'});
pos = get(gca, 'Position');
pos(1) = pos(1)-0.1;
set(gca, 'Position',pos);




