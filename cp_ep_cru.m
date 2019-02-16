% Examine relationships between CRU climate and CP/EP ENSO index over
% 1951-2016 period

syear = 1951;
eyear = 2016;
latlim = [-80 80];
lonlim = [-180 180];
worldland = shaperead('landareas','UseGeoCoords', true);
clr = cbrewer('div','RdBu',12);

%% EP and CP ENSO indices
load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

%% CRU climate data
tmx = ncread('C:\Users\dannenberg\Documents\Data_Analysis\CRU\cru_ts4.01.1901.2016.tmx.dat.nc', 'tmx');
tmx = permute(tmx, [3 2 1]);
tmx(tmx==9.9692e+36) = NaN;
tmn = ncread('C:\Users\dannenberg\Documents\Data_Analysis\CRU\cru_ts4.01.1901.2016.tmn.dat.nc', 'tmn');
tmn = permute(tmn, [3 2 1]);
tmn(tmn==9.9692e+36) = NaN;
pre = ncread('C:\Users\dannenberg\Documents\Data_Analysis\CRU\cru_ts4.01.1901.2016.pre.dat.nc', 'pre');
pre = permute(pre, [3 2 1]);
pre(pre==9.9692e+36) = NaN;
vap = ncread('C:\Users\dannenberg\Documents\Data_Analysis\CRU\cru_ts4.01.1901.2016.vap.dat.nc', 'vap');
vap = permute(vap, [3 2 1]);
vap(vap==9.9692e+36) = NaN;

cru_yr = reshape(repmat(1901:2016, 12, 1), [], 1);
cru_mo = repmat(1:12, 1, length(1901:2016))';
cru_lat = double(ncread('C:\Users\dannenberg\Documents\Data_Analysis\CRU\cru_ts4.01.1901.2016.vap.dat.nc', 'lat'));
cru_lon = double(ncread('C:\Users\dannenberg\Documents\Data_Analysis\CRU\cru_ts4.01.1901.2016.vap.dat.nc', 'lon'));

estmx = 10 * 0.611 * exp(17.502*tmx ./ (tmx+240.97)); % e(s) (hPa) at Tmax
estmn = 10 * 0.611 * exp(17.502*tmn ./ (tmn+240.97)); % e(s) (hPa) at Tmin
es = (estmx + estmn) / 2;

vpd = es - vap;

clear estmx estmn es;
[nt, ny, nx] = size(vpd);

%% CRU to three month means/sums
% 3 month running sums of precip
windowSize = 3;
b = ones(1,windowSize);
a = 1;
pre = filter(b, a, pre, [], 1);
clear b a windowSize;

% 3 month running means of precip
windowSize = 3;
b = ones(1,windowSize) / windowSize;
a = 1;
tmn = filter(b, a, tmn, [], 1);
tmx = filter(b, a, tmx, [], 1);
vpd = filter(b, a, vpd, [], 1);
clear b a windowSize;

%% Minimum temperature ~ EPI/CPI
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 5.75];
set(h, 'defaultAxesColorOrder',[0 0 0; 0,0,0]/255)

% DJF
clim = tmn(cru_yr>=syear & cru_yr<=eyear & cru_mo==2, :, :);

[r, p] = corr(epi, reshape(clim, [], nx*ny));
subplot(3,2,1)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
text(-2.2,1.3,'A', 'FontSize',12);
text(-3.5,0,'Dec-Feb', 'FontSize',12, 'Rotation',90, 'HorizontalAlignment',...
    'center', 'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold');
ttl = title('Eastern Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 1.75;

[r, p] = corr(cpi, reshape(clim, [], nx*ny));
subplot(3,2,2)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
text(-2.2,1.3,'B', 'FontSize',12);
ttl = title('Central Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 1.75;

% MAM
clim = tmn(cru_yr>=syear & cru_yr<=eyear & cru_mo==5, :, :);

[r, p] = corr(epi, reshape(clim, [], nx*ny));
subplot(3,2,3)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.03;
set(gca, 'Position',pos);
text(-2.2,1.3,'C', 'FontSize',12);
text(-3.5,0,'Mar-May', 'FontSize',12, 'Rotation',90, 'HorizontalAlignment',...
    'center', 'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold');

[r, p] = corr(cpi, reshape(clim, [], nx*ny));
subplot(3,2,4)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.03;
set(gca, 'Position',pos);
text(-2.2,1.3,'D', 'FontSize',12);

% JJA
clim = tmn(cru_yr>=syear & cru_yr<=eyear & cru_mo==8, :, :);

[r, p] = corr(epi, reshape(clim, [], nx*ny));
subplot(3,2,5)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.06;
set(gca, 'Position',pos);
text(-2.2,1.3,'E', 'FontSize',12);
text(-3.5,0,'Jun-Aug', 'FontSize',12, 'Rotation',90, 'HorizontalAlignment',...
    'center', 'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold');

[r, p] = corr(cpi, reshape(clim, [], nx*ny));
subplot(3,2,6)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.06;
set(gca, 'Position',pos);
text(-2.2,1.3,'F', 'FontSize',12);

cb = colorbar('southoutside');
cb.Position = [0.15    0.1    0.7    0.04];
cb.Ticks = -0.6:0.1:0.6;
cb.TickLength = 0.054;
set(cb, 'FontSize',8)
ylb = ylabel(cb, 'Pearson''s \itR', 'FontSize',11);
cbarrow;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/epi-cpi-tmn.tif')
close all;

%% Maximum temperature ~ EPI/CPI
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 5.75];
set(h, 'defaultAxesColorOrder',[0 0 0; 0,0,0]/255)

% DJF
clim = tmx(cru_yr>=syear & cru_yr<=eyear & cru_mo==2, :, :);

[r] = corr(epi, reshape(clim, [], nx*ny));
subplot(3,2,1)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
text(-2.2,1.3,'A', 'FontSize',12);
text(-3.5,0,'Dec-Feb', 'FontSize',12, 'Rotation',90, 'HorizontalAlignment',...
    'center', 'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold');
ttl = title('Eastern Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 1.75;

[r] = corr(cpi, reshape(clim, [], nx*ny));
subplot(3,2,2)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
text(-2.2,1.3,'B', 'FontSize',12);
ttl = title('Central Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 1.75;

% MAM
clim = tmx(cru_yr>=syear & cru_yr<=eyear & cru_mo==5, :, :);

[r] = corr(epi, reshape(clim, [], nx*ny));
subplot(3,2,3)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.03;
set(gca, 'Position',pos);
text(-2.2,1.3,'C', 'FontSize',12);
text(-3.5,0,'Mar-May', 'FontSize',12, 'Rotation',90, 'HorizontalAlignment',...
    'center', 'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold');

[r] = corr(cpi, reshape(clim, [], nx*ny));
subplot(3,2,4)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.03;
set(gca, 'Position',pos);
text(-2.2,1.3,'D', 'FontSize',12);

% JJA
clim = tmx(cru_yr>=syear & cru_yr<=eyear & cru_mo==8, :, :);

[r] = corr(epi, reshape(clim, [], nx*ny));
subplot(3,2,5)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.06;
set(gca, 'Position',pos);
text(-2.2,1.3,'E', 'FontSize',12);
text(-3.5,0,'Jun-Aug', 'FontSize',12, 'Rotation',90, 'HorizontalAlignment',...
    'center', 'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold');

[r] = corr(cpi, reshape(clim, [], nx*ny));
subplot(3,2,6)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.06;
set(gca, 'Position',pos);
text(-2.2,1.3,'F', 'FontSize',12);

cb = colorbar('southoutside');
cb.Position = [0.15    0.1    0.7    0.04];
cb.Ticks = -0.6:0.1:0.6;
cb.TickLength = 0.054;
set(cb, 'FontSize',8)
ylb = ylabel(cb, 'Pearson''s \itR', 'FontSize',11);
cbarrow;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/epi-cpi-tmx.tif')
close all;

%% Precipitation ~ EPI/CPI
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 5.75];
set(h, 'defaultAxesColorOrder',[0 0 0; 0,0,0]/255)

% DJF
clim = pre(cru_yr>=syear & cru_yr<=eyear & cru_mo==2, :, :);

[r] = corr(epi, reshape(clim, [], nx*ny));
subplot(3,2,1)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
text(-2.2,1.3,'A', 'FontSize',12);
text(-3.5,0,'Dec-Feb', 'FontSize',12, 'Rotation',90, 'HorizontalAlignment',...
    'center', 'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold');
ttl = title('Eastern Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 1.75;

[r] = corr(cpi, reshape(clim, [], nx*ny));
subplot(3,2,2)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
text(-2.2,1.3,'B', 'FontSize',12);
ttl = title('Central Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 1.75;

% MAM
clim = pre(cru_yr>=syear & cru_yr<=eyear & cru_mo==5, :, :);

[r] = corr(epi, reshape(clim, [], nx*ny));
subplot(3,2,3)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.03;
set(gca, 'Position',pos);
text(-2.2,1.3,'C', 'FontSize',12);
text(-3.5,0,'Mar-May', 'FontSize',12, 'Rotation',90, 'HorizontalAlignment',...
    'center', 'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold');

[r] = corr(cpi, reshape(clim, [], nx*ny));
subplot(3,2,4)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.03;
set(gca, 'Position',pos);
text(-2.2,1.3,'D', 'FontSize',12);

% JJA
clim = pre(cru_yr>=syear & cru_yr<=eyear & cru_mo==8, :, :);

[r] = corr(epi, reshape(clim, [], nx*ny));
subplot(3,2,5)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.06;
set(gca, 'Position',pos);
text(-2.2,1.3,'E', 'FontSize',12);
text(-3.5,0,'Jun-Aug', 'FontSize',12, 'Rotation',90, 'HorizontalAlignment',...
    'center', 'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold');

[r] = corr(cpi, reshape(clim, [], nx*ny));
subplot(3,2,6)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.06;
set(gca, 'Position',pos);
text(-2.2,1.3,'F', 'FontSize',12);

cb = colorbar('southoutside');
cb.Position = [0.15    0.1    0.7    0.04];
cb.Ticks = -0.6:0.1:0.6;
cb.TickLength = 0.054;
set(cb, 'FontSize',8)
ylb = ylabel(cb, 'Pearson''s \itR', 'FontSize',11);
cbarrow;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/epi-cpi-pre.tif')
close all;

%% VPD ~ EPI/CPI
h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 6.5 5.75];
set(h, 'defaultAxesColorOrder',[0 0 0; 0,0,0]/255)

% DJF
clim = vpd(cru_yr>=syear & cru_yr<=eyear & cru_mo==2, :, :);

[r] = corr(epi, reshape(clim, [], nx*ny));
subplot(3,2,1)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
text(-2.2,1.3,'A', 'FontSize',12);
text(-3.5,0,'Dec-Feb', 'FontSize',12, 'Rotation',90, 'HorizontalAlignment',...
    'center', 'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold');
ttl = title('Eastern Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 1.75;

[r] = corr(cpi, reshape(clim, [], nx*ny));
subplot(3,2,2)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
text(-2.2,1.3,'B', 'FontSize',12);
ttl = title('Central Pacific ENSO', 'FontSize',12);
ttl.Position(2) = 1.75;

% MAM
clim = vpd(cru_yr>=syear & cru_yr<=eyear & cru_mo==5, :, :);

[r] = corr(epi, reshape(clim, [], nx*ny));
subplot(3,2,3)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.03;
set(gca, 'Position',pos);
text(-2.2,1.3,'C', 'FontSize',12);
text(-3.5,0,'Mar-May', 'FontSize',12, 'Rotation',90, 'HorizontalAlignment',...
    'center', 'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold');

[r] = corr(cpi, reshape(clim, [], nx*ny));
subplot(3,2,4)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.03;
set(gca, 'Position',pos);
text(-2.2,1.3,'D', 'FontSize',12);

% JJA
clim = vpd(cru_yr>=syear & cru_yr<=eyear & cru_mo==8, :, :);

[r] = corr(epi, reshape(clim, [], nx*ny));
subplot(3,2,5)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(gca, flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.06;
set(gca, 'Position',pos);
text(-2.2,1.3,'E', 'FontSize',12);
text(-3.5,0,'Jun-Aug', 'FontSize',12, 'Rotation',90, 'HorizontalAlignment',...
    'center', 'VerticalAlignment','middle', 'FontSize',12, 'FontWeight','bold');

[r] = corr(cpi, reshape(clim, [], nx*ny));
subplot(3,2,6)
ax = axesm('winkel','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','off',...
        'ParallelLabel','off','GLineWidth',0.5,'Frame','on','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel',min(latlim)+0.11);
axis off;
axis image;
surfm(cru_lat, cru_lon, reshape(r, ny, nx))
caxis([-0.6 0.6])
colormap(flipud(clr));
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(2) = pos(2)+0.06;
set(gca, 'Position',pos);
text(-2.2,1.3,'F', 'FontSize',12);

cb = colorbar('southoutside');
cb.Position = [0.15    0.1    0.7    0.04];
cb.Ticks = -0.6:0.1:0.6;
cb.TickLength = 0.054;
set(cb, 'FontSize',8)
ylb = ylabel(cb, 'Pearson''s \itR', 'FontSize',11);
cbarrow;

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/epi-cpi-vpd.tif')
close all;


