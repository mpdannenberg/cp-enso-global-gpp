% Create EP/CP Niño dataset from COBE-2 SST
t0 = datenum(1891, 1, 1);
latlim = [-20 20];
lonlim = [120 290];
syear = 1951;

warning('off');

windowSize = 6; % Number of months over which to average SSTs
endMonth = 2; % End month of integration period

%% Load COBE-2 SST data
% info = ncinfo('./data/sst.mon.mean.nc');
lat = double(ncread('./data/sst.mon.mean.nc','lat'));
lon = double(ncread('./data/sst.mon.mean.nc','lon'));
t = ncread('./data/sst.mon.mean.nc','time');
sst = double(ncread('./data/sst.mon.mean.nc','sst'));
sst = permute(sst, [2 1 3]);
sst(sst>40) = NaN;
[yr, mo] = datevec(t0 + t); clear t0 t;

%% Need to do monthly anomalies
ssta = NaN(size(sst));
for i = 1:12
    
    x = sst(:, :, mo==i);
    xclim = sst(:, :, mo==i & yr>=1981 & yr<=2010);
    xm = mean(xclim, 3);
    
    ssta(:,:,mo==i) = x - repmat(xm, 1,1,size(x, 3));
    
end
clear x xclim xm i sst;

lat_weight = repmat(sqrt(cosd(lat)), 1, length(lon));

%% Regress SST data onto Nino-4 and Nino-1+2 indices to get independent CP- and EP- SST anomalies
nino4 = ssta(lat>=-5 & lat<=5, lon>=160 & lon<=210, :);
nino4 = reshape(permute(nino4, [3 1 2]), length(yr), []);
w = repmat(reshape(lat_weight(lat>=-5 & lat<=5, lon>=160 & lon<=210), 1, []), length(yr), 1);
% nino4 = mean(nino4, 2);
nino4 = sum(nino4.*w, 2) ./ sum(w, 2);
clear w;

nino12 = ssta(lat>=-10 & lat<=0, lon>=270 & lon<=280, :);
nino12 = reshape(permute(nino12, [3 1 2]), length(yr), []);
idx = sum(isnan(nino12))==0;
w = repmat(reshape(lat_weight(lat>=-10 & lat<=0, lon>=270 & lon<=280), 1, []), length(yr), 1);
% nino12 = nanmean(nino12, 2);
nino12 = sum(nino12(:,idx).*w(:,idx), 2) ./ sum(w(:,idx), 2);
clear w idx;

sst_tp = ssta(lat>=min(latlim) & lat<=max(latlim), lon>=min(lonlim) & lon<=max(lonlim), :);
[ny, nx, ~] = size(sst_tp);
sst_tp = reshape(permute(sst_tp, [3 1 2]), length(yr), []);
w_tp = reshape(lat_weight(lat>=min(latlim) & lat<=max(latlim), lon>=min(lonlim) & lon<=max(lonlim)), 1, []);
sst_idx = sum(isnan(sst_tp)) == 0;
sst_tp = sst_tp(:, sst_idx);
w_tp = w_tp(sst_idx);
sst_ep = NaN(size(sst_tp));
sst_cp = NaN(size(sst_tp));
for i = 1:size(sst_tp, 2)
    y = sst_tp(:, i);
    lm_ep = fitlm(nino4, y);
    sst_ep(:, i) = lm_ep.Residuals.Raw;
    lm_cp = fitlm(nino12, y);
    sst_cp(:, i) = lm_cp.Residuals.Raw;
end
clear y lm_ep lm_cp i sst_tp;

[ep_coef, ep_pc, ep_latent] = pca(sst_ep, 'VariableWeights',w_tp);
[cp_coef, cp_pc, cp_latent] = pca(sst_cp, 'VariableWeights',w_tp);
ep_eof = NaN(size(sst_idx));
cp_eof = NaN(size(sst_idx));
ep_eof(sst_idx) = ep_coef(:, 1);
cp_eof(sst_idx) = cp_coef(:, 1);
ep_eof = reshape(ep_eof, ny, nx);
cp_eof = reshape(cp_eof, ny, nx);
ep_idx = ep_pc(:, 1);
cp_idx = cp_pc(:, 1);
clear sst_idx w_tp ny nx;

%% Get spatial patterns by correlating SSTa against CPI and EPI
[ny, nx, nt] = size(ssta);
cpi_r = NaN(ny, nx);
epi_r = NaN(ny, nx);
for i = 1:ny
    for j = 1:nx
        y = squeeze(ssta(i, j, :));
        cpi_r(i,j) = corr(y, cp_idx);
        epi_r(i,j) = corr(y, ep_idx);
    end
end
clear i j y;


%% Filter to 6-month mean CPI and EPI
b = ones(1,windowSize)/windowSize;
a = 1;
cpi = filter(b, a, cp_idx, [], 1);
epi = filter(b, a, ep_idx, [], 1);
clear b a windowSize;

%% Include only time period of interest
idx = yr >= syear & mo == endMonth; % Sep-Feb 1951-2016
cpi = cpi(idx);
epi = epi(idx);
yr = yr(idx);

%% Standardize CPI and EPI
cpi = (cpi - mean(cpi(yr>=1981 & yr<=2010))) / std(cpi(yr>=1981 & yr<=2010));
epi = (epi - mean(epi(yr>=1981 & yr<=2010))) / std(epi(yr>=1981 & yr<=2010));
clear cp_coef cp_eof cp_idx cp_latent cp_pc ep_coef ep_eof ep_idx ep_latent ep_pc;

save('./data/cpi_epi_1951-2016.mat', 'cpi','epi','yr', '-v7.3');

%% Figure
% ENSO regions
nino3_lat = [-5 5 5 -5 -5]; nino3_lon = [-150 -150 -90 -90 -150];
nino12_lat = [-10 0 0 -10 -10]; nino12_lon = [-90 -90 -80 -80 -90];
nino4_lat = [-5 5 5 -5 -5]; nino4_lon = [160 160 -150 -150 160];
nino34_lat = [-5 5 5 -5 -5]; nino34_lon = [-170 -170 -120 -120 -170];

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 6.4];

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

latlim = [-25 25];

ax = tight_subplot(3,1,0,0.1,[0.06 0.14]);

axes(ax(1))
axesm('mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
            'on','PLineLocation',10,'MLineLocation',30,'MeridianLabel','on',...
            'ParallelLabel','on','GLineWidth',0.5,'Frame','on','FFaceColor',...
            [0.7 0.7 0.7], 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
            'FLineWidth',1, 'FontColor',[0.2 0.2 0.2], 'MLabelParallel',25,...
            'FEdgeColor','w','FontSize',7);
axis off;
axis image;
contourm(lat, lon, epi_r, 'Fill','on', 'LevelList',-1:0.2:1);
colormap(gca, flipud(clr));
caxis([-1 1]);
plotm(nino3_lat, nino3_lon, '-', 'Color',[0.3 0.3 0.3], 'LineWidth',1.5);
plotm(nino4_lat, nino4_lon, '-', 'Color',[0.3 0.3 0.3], 'LineWidth',1.5);
plotm(nino12_lat, nino12_lon, '-', 'Color',[0.3 0.3 0.3], 'LineWidth',1.5);
plotm(nino34_lat, nino34_lon, '-', 'Color','k', 'LineWidth',2);
textm(0,175, 'Niño 4', 'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.3 0.3 0.3], 'FontWeight','bold')
textm(0,-105, 'Niño 3', 'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.3 0.3 0.3], 'FontWeight','bold')
textm(-10,-85, 'Niño 1+2', 'HorizontalAlignment','center','VerticalAlignment','top','Color',[0.3 0.3 0.3], 'FontWeight','bold')
textm(5,-145, 'Niño 3.4', 'HorizontalAlignment','center','VerticalAlignment','bottom','Color','k', 'FontWeight','bold')
textm(20,125,'A. Eastern Pacific ENSO Index', 'FontSize',12);

axes(ax(2))
axesm('mercator','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
            'on','PLineLocation',10,'MLineLocation',30,'MeridianLabel','off',...
            'ParallelLabel','on','GLineWidth',0.5,'Frame','on','FFaceColor',...
            [0.7 0.7 0.7], 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
            'FLineWidth',1, 'FontColor',[0.2 0.2 0.2], 'MLabelParallel',25,...
            'FEdgeColor','w','FontSize',7);
axis off;
axis image;
contourm(lat, lon, cpi_r, 'Fill','on', 'LevelList',-1:0.2:1);
colormap(gca, flipud(clr));
caxis([-1 1]);
plotm(nino3_lat, nino3_lon, '-', 'Color',[0.3 0.3 0.3], 'LineWidth',1.5);
plotm(nino4_lat, nino4_lon, '-', 'Color',[0.3 0.3 0.3], 'LineWidth',1.5);
plotm(nino12_lat, nino12_lon, '-', 'Color',[0.3 0.3 0.3], 'LineWidth',1.5);
plotm(nino34_lat, nino34_lon, '-', 'Color','k', 'LineWidth',2);
textm(0,175, 'Niño 4', 'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.3 0.3 0.3], 'FontWeight','bold')
textm(0,-105, 'Niño 3', 'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.3 0.3 0.3], 'FontWeight','bold')
textm(-10,-85, 'Niño 1+2', 'HorizontalAlignment','center','VerticalAlignment','top','Color',[0.3 0.3 0.3], 'FontWeight','bold')
textm(5,-145, 'Niño 3.4', 'HorizontalAlignment','center','VerticalAlignment','bottom','Color','k', 'FontWeight','bold')
textm(20,125,'B. Central Pacific ENSO Index', 'FontSize',12);

cb = colorbar('eastoutside');
cb.Position = [0.87 0.367 0.04 0.535];
cb.TickLength = 0.08;
xlabel(cb, 'Pearson''s R', 'FontSize',11);

axes(ax(3))
plot(yr, epi, 'k-', 'LineWidth',2);
hold on;
plot(yr, cpi, '-', 'Color',[0.5 0.5 0.5], 'LineWidth',2);
set(gca, 'YLim',[-3.5 3.5],'YTick',-3:3, 'XLim',[1950 2016], 'TickLength',[0 0],...
    'XColor',[0.2 0.2 0.2], 'YColor',[0.2 0.2 0.2])
grid on;
hold off;
xlabel('Year')
text(1952,2.5, 'C. Time series', 'FontSize',12)
lgd = legend('Eastern Pacific ENSO Index (EPI)','Central Pacific ENSO Index (CPI)','Location','southwest');
legend('boxoff');
lgd.FontSize = 11;
lgd.Position = [0.22 0.1 0.1161 0.0757];

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/epi-cpi-1951-2016.tif')
close all;
