% Read MODIS land cover and CRU temperature and convert to Ahlstrom 2015
% classes

%% MCD12C1, aggregated to three classes (forest, savanna/shrub, grass/crop)
mcd12c1 = double(hdfread('./data/MCD12C1.A2001001.006.2018053185512.hdf', 'Majority_Land_Cover_Type_3'));
mcd12c1_lat = (90-0.025):-0.05:(-90+0.025);
mcd12c1_lon = (-180+0.025):0.05:(180-0.025);

biome = zeros(size(mcd12c1));
biome(mcd12c1 >= 5 & mcd12c1 <= 8) = 1; % Forest
biome(mcd12c1 == 2 | mcd12c1 == 4) = 2; % Savanna & Shrub
biome(mcd12c1 == 1 | mcd12c1 == 3) = 3; % Grass & Crop

clear mcd12c1;

%% CRU temperature
tmp = ncread('D:\Data_Analysis\CRU\cru_ts4.01.1901.2016.tmp.dat.nc', 'tmp');
tmp = permute(tmp, [3 2 1]);
tmp(tmp==9.9692e+36) = NaN;

% 12-month running means of precip
windowSize = 12;
b = ones(1,windowSize) / windowSize;
a = 1;
tmp = filter(b, a, tmp, [], 1);
clear b a windowSize;

cru_yr = reshape(repmat(1901:2016, 12, 1), [], 1);
cru_mo = repmat(1:12, 1, length(1901:2016))';
cru_lat = flipud(double(ncread('D:\Data_Analysis\CRU\cru_ts4.01.1901.2016.tmp.dat.nc', 'lat')));
cru_lon = double(ncread('D:\Data_Analysis\CRU\cru_ts4.01.1901.2016.tmp.dat.nc', 'lon'));

tmp = flipud(squeeze(mean(tmp(cru_yr>=1981 & cru_yr<=2010 & cru_mo==12, :, :), 1)));

%% Loop through each point and get class, lat, and MAT to define biome
biome_half = NaN(size(tmp));

for i = 1:length(cru_lat)
    for j = 1:length(cru_lon)
        
        lat = cru_lat(i); lon = cru_lon(j);
        latidx = mcd12c1_lat >= (lat - 0.25) & mcd12c1_lat <= (lat + 0.25);
        lonidx = mcd12c1_lon >= (lon - 0.25) & mcd12c1_lon <= (lon + 0.25);
        
        biome_sub = reshape(biome(latidx, lonidx), 1, []);
        biome_maj = mode(biome_sub);
        
        mat = tmp(i, j);
        
        if biome_maj==1 && mat>=18
            biome_half(i, j) = 1; % Tropical forest
        elseif biome_maj==1 && mat<18
            biome_half(i, j) = 2; % Extratropical forest
        elseif biome_maj==2 && lat<45
            biome_half(i, j) = 5; % Semiarid
        elseif biome_maj==2 && lat>45
            biome_half(i, j) = 3; % Tundra & Arctic Shrub
        elseif biome_maj==3
            biome_half(i, j) = 4; % Grass and crop
        end
        
    end
end
lat = cru_lat; lon = cru_lon;
clear i j latidx lonidx biome_sub biome_maj mat cru_* tmp mcd12c1* biome;

save('./data/ahlstrom_regions.mat');

%% Make figure
clr = wesanderson('aquatic4');

latlim = [-65 75];
lonlim = [-180 180];
worldland = shaperead('landareas','UseGeoCoords', true);

h = figure('Color','w');
h.Units = 'inches';
h.Position = [1 1 7 3.5];

ax = axesm('miller','MapLatLimit',latlim,'MapLonLimit',lonlim,'grid',...
        'on','PLineLocation',30,'MLineLocation',60,'MeridianLabel','on',...
        'ParallelLabel','on','GLineWidth',0.5,'Frame','off','FFaceColor',...
        'none', 'FontName', 'Helvetica','GColor',[0.6 0.6 0.6],...
        'FLineWidth',1, 'FontColor',[0.5 0.5 0.5], 'MLabelParallel','north');
axis off;
axis image;
surfm(lat, lon, biome_half)
caxis([0.5 5.5])
colormap(gca, clr);
geoshow(worldland,'FaceColor','none','EdgeColor',[0.6 0.6 0.6])
pos = get(gca, 'Position');
pos(1) = pos(1)-0.04;
set(gca, 'Position',pos);
% text(-2.2,1.3,'C', 'FontSize',12);
cb = colorbar('eastoutside');
cb.TickLength = 0;
cb.TickLabels = {'Tropical Forest','Extratropical Forest','Tundra/Arctic Shrub','Grass/Crop','Semiarid'};

% Amazon
rlim = [-30 10; -80 -35];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)]), mean([rlim(2,1) rlim(2,2)]),'Amazonia',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Sahel
rlim = [5 15; -20 50];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)]), mean([rlim(2,1) rlim(2,2)]),'Sahel',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Tropical & subtropical Africa
rlim = [-30 5; 8 42];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)])+10, mean([rlim(2,1) rlim(2,2)]),'Tropical &',...
    'FontSize',7, 'HorizontalAlignment','center','VerticalAlignment','middle')
textm(mean([rlim(1,1) rlim(1,2)]), mean([rlim(2,1) rlim(2,2)]),'subtropical',...
    'FontSize',7, 'HorizontalAlignment','center','VerticalAlignment','middle')
textm(mean([rlim(1,1) rlim(1,2)])-10, mean([rlim(2,1) rlim(2,2)]),'Africa',...
    'FontSize',7, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Australia
rlim = [-40 -10; 110 155];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)]), mean([rlim(2,1) rlim(2,2)]),'Australia',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Western North America
rlim = [20 70; -165 -100];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)])+5, mean([rlim(2,1) rlim(2,2)]),'Western',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')
textm(mean([rlim(1,1) rlim(1,2)])-5, mean([rlim(2,1) rlim(2,2)]),'North America',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Eastern U.S.
rlim = [25 50; -100 -60];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)])+5, mean([rlim(2,1) rlim(2,2)]),'Eastern',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')
textm(mean([rlim(1,1) rlim(1,2)])-5, mean([rlim(2,1) rlim(2,2)]),'U.S.',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Europe
rlim = [35 60; -10 40];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)]), mean([rlim(2,1) rlim(2,2)]),'Europe',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

% Central Asia & SW Russia
rlim = [45 65; 50 100];
plotm([rlim(1,1) rlim(1,2)], [rlim(2,1) rlim(2,1)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,1)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,2) rlim(1,2)], [rlim(2,1) rlim(2,2)], 'k-','LineWidth',1.5)
plotm([rlim(1,1) rlim(1,2)], [rlim(2,2) rlim(2,2)], 'k-','LineWidth',1.5)
textm(mean([rlim(1,1) rlim(1,2)])+5, mean([rlim(2,1) rlim(2,2)]),'SW Russia &',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')
textm(mean([rlim(1,1) rlim(1,2)])-5, mean([rlim(2,1) rlim(2,2)]),'central Asia',...
    'FontSize',9, 'HorizontalAlignment','center','VerticalAlignment','middle')

set(gcf,'PaperPositionMode','auto')
print('-dtiff','-f1','-r300','./output/supplemental-region-map-v2.tif')
close all;



