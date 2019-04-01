% Get 1982-2015 GPP data from CCW and format for comparison to other data
% sources (MsTMIP, MOD17, and FluxCom)

%% Units for final GPP data (original: g C m-2 month-1)
% monthly gridded GPP: kgC m-2 day-1
% monthly global GPP: TgC day-1
% annual gridded GPP: kgC m-2 year-1
% annual global GPP: TgC year-1

%% Setup
syear = 1982; % First year of analysis
eyear = 2016; % Last year of analysis
scale = 10^-9; % kg --> Tg
e = referenceEllipsoid('World Geodetic System 1984');
biomes = ncread('D:\Data_Analysis\MsTMIP\mstmip_driver_global_hd_biome_v1.nc4', 'biome_type');

%% Load CCW GPP data - gC m-2 month-1
years = syear:eyear;
cd('D:\Data_Analysis\CCW\GeoTIFFs');
[~, R] = geotiffread('GPP_GV42_1982.tif');
lat = (90 - R.CellExtentInLatitude/2):-R.CellExtentInLatitude:(-90 + R.CellExtentInLatitude/2);
lon = (-180 + R.CellExtentInLatitude/2):R.CellExtentInLatitude:(180 - R.CellExtentInLatitude/2);
ny = length(lat); nx = length(lon);

yr = reshape(repmat(syear:eyear, 12, 1), [], 1);
nt = length(yr);
mo = repmat(1:12, 1, length(syear:eyear))';
ndys = repmat([31,28,31,30,31,30,31,31,30,31,30,31]', 1, ny, nx); % Number of days, excluding leap years
ndys = permute(ndys, [2 3 1]);

GPP = NaN(ny, nx, nt);

for i = 1:length(years)
    
    fn = ['GPP_GV42_',num2str(years(i)),'.tif'];
    gpp = 0.001 * geotiffread(fn) ./ ndys; % from gC m-2 month-1 --> kgC m-2 day-1
    
    GPP(:, :, yr==years(i)) = gpp;
    
end
clear i gpp fn;
cd('D:\Publications\Dannenberg_et_al_CPElNinoGlobalGPP');

%% Convert CCW from 1/12 deg to 1/2 deg
% Initialize new 1/2 deg grid
newlat = 89.75:-0.5:-89.75;
newlon = -179.75:0.5:179.75;
GPP_halfdeg = NaN(length(newlat), length(newlon), nt);

% Area (in m^2) of the 1/12 deg grid
[LON, LAT] = meshgrid(lon, lat);
area = areaquad(reshape(LAT-(1/24),[],1),reshape(LON-(1/24),[],1),reshape(LAT+(1/24),[],1),reshape(LON+(1/24),[],1),e);
area12 = reshape(area, ny, nx); 

% Area (in m^2) of the 1/2 deg grid
[LON, LAT] = meshgrid(newlon, newlat);
area = areaquad(reshape(LAT-(1/4),[],1),reshape(LON-(1/4),[],1),reshape(LAT+(1/4),[],1),reshape(LON+(1/4),[],1),e);
area = reshape(area, length(newlat), length(newlon)); 
clear LON LAT;

% Loop through new grid and calculate GPP at 1/2 deg scale (kgC m-2 day-1)
for i = 1:length(newlat)
    for j = 1:length(newlon)
        
        latidx = lat>=(newlat(i)-0.25) & lat<=(newlat(i)+0.25);
        lonidx = lon>=(newlon(j)-0.25) & lon<=(newlon(j)+0.25);
        gpp = reshape(GPP(latidx, lonidx, :) .* repmat(area12(latidx, lonidx),1,1,nt), [], nt);
        
        temp = sum(gpp) / area(i,j);
        GPP_halfdeg(i,j,:) = temp;
        
    end
end

GPP_halfdeg(repmat(biomes',1,1,nt)==0) = NaN;

lat = newlat; ny = length(lat);
lon = newlon; nx = length(lon);
GPP = GPP_halfdeg;
ndys = repmat([31,28,31,30,31,30,31,31,30,31,30,31]', length(years), ny, nx); % Number of days, excluding leap years
ndys = permute(ndys, [2 3 1]);

clear newlat newlon GPP_halfdeg i j area12 gpp latidx lonidx biomes;

%% Aggregate to monthly and annual scales
% Monthly gridded mean
GPP_monthly = GPP .* ndys; % daily average --> monthly total GPP (kgC m-2 month-1)

% Multimodel annual gridded mean
windowSize = 12;
b = ones(1,windowSize);
a = 1;
GPP_annual = filter(b, a, GPP_monthly, [], 3); % 12-month running sums (kgC m-2 yr-1)
GPP_annual = GPP_annual(:, :, mo==12); % Get calendar year sum
clear GPP_monthly a b windowSize ndys;

%% Calculate global GPP at monthly and annual scale
yrs = years;
GPP_global_monthly = NaN(size(GPP_annual, 3), 12);
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_global_monthly(i,j) = nansum(nansum( gpp.*area )) * scale; % TgC day-1
    end
end

GPP_global_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_global_annual(i) = nansum(nansum( gpp.*area )) * scale; % TgC yr-1
end

clear i j k gpp yrs area e eyear syear R nt nx ny scale;

GPP_monthly = GPP; clear GPP; % Monthly GPP in kgC m-2 day-1

save('./data/gpp_ccw.mat');


%% Calculate regional GPP at monthly and annual scale

yrs = years;
GPP = GPP_monthly;
scale = 10^-9; % kg --> Tg
e = referenceEllipsoid('World Geodetic System 1984');
[LON, LAT] = meshgrid(lon, lat);
area = areaquad(reshape(LAT-(1/4),[],1),reshape(LON-(1/4),[],1),reshape(LAT+(1/4),[],1),reshape(LON+(1/4),[],1),e);
area = reshape(area, length(lat), length(lon)); 
clear LON LAT;

% Amazon
rlim = [-30 10; -80 -35];
GPP_amazon_monthly = NaN(size(GPP_annual, 3), 12);
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_amazon_monthly(i,j) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC day-1
    end
end
GPP_amazon_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_amazon_annual(i) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC yr-1
end

% Sahel
rlim = [5 15; -20 50];
GPP_sahel_monthly = NaN(size(GPP_annual, 3), 12);
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_sahel_monthly(i,j) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC day-1
    end
end
GPP_sahel_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_sahel_annual(i) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC yr-1
end

% Tropical and subtropical Africa
rlim = [-30 5; 10 40];
GPP_africa_monthly = NaN(size(GPP_annual, 3), 12);
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_africa_monthly(i,j) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC day-1
    end
end
GPP_africa_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_africa_annual(i) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC yr-1
end

% Australia
rlim = [-40 -10; 110 155];
GPP_austr_monthly = NaN(size(GPP_annual, 3), 12);
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_austr_monthly(i,j) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC day-1
    end
end
GPP_austr_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_austr_annual(i) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC yr-1
end

% Western North America
rlim = [20 70; -165 -100];
GPP_westna_monthly = NaN(size(GPP_annual, 3), 12);
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_westna_monthly(i,j) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC day-1
    end
end
GPP_westna_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_westna_annual(i) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC yr-1
end

% Eastern U.S.
rlim = [25 50; -100 -60];
GPP_eastus_monthly = NaN(size(GPP_annual, 3), 12);
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_eastus_monthly(i,j) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC day-1
    end
end
GPP_eastus_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_eastus_annual(i) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC yr-1
end

% Europe
rlim = [35 60; -10 40];
GPP_europe_monthly = NaN(size(GPP_annual, 3), 12);
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_europe_monthly(i,j) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC day-1
    end
end
GPP_europe_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_europe_annual(i) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC yr-1
end

% Central Asia and southern Russia
rlim = [45 65; 50 100];
GPP_casia_monthly = NaN(size(GPP_annual, 3), 12);
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_casia_monthly(i,j) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC day-1
    end
end
GPP_casia_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_casia_annual(i) = nansum(nansum( gpp(latidx,lonidx).*area(latidx,lonidx) )) * scale; % TgC yr-1
end

clear i j k gpp yrs area e eyear syear R nt nx ny scale GPP GPP_annual GPP_global* latidx lonidx lat lon yr mo rlim;

save('./data/gpp_ccw_regional.mat');


%% Calculate regional GPP (Ahlstrom et al. 2015 version) at monthly and annual scale

yrs = years;
GPP = GPP_monthly;
scale = 10^-9; % kg --> Tg
e = referenceEllipsoid('World Geodetic System 1984');
[LON, LAT] = meshgrid(lon, lat);
area = areaquad(reshape(LAT-(1/4),[],1),reshape(LON-(1/4),[],1),reshape(LAT+(1/4),[],1),reshape(LON+(1/4),[],1),e);
area = reshape(area, length(lat), length(lon)); 
clear LON LAT;

load ./data/ahlstrom_regions.mat;

% Tropical forest
GPP_tropical_monthly = NaN(size(GPP_annual, 3), 12);
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_tropical_monthly(i,j) = nansum(nansum( gpp(biome_half==1).*area(biome_half==1) )) * scale; % TgC day-1
    end
end
GPP_tropical_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_tropical_annual(i) = nansum(nansum( gpp(biome_half==1).*area(biome_half==1) )) * scale; % TgC yr-1
end

% Extratropical forest
GPP_extratropical_monthly = NaN(size(GPP_annual, 3), 12);
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_extratropical_monthly(i,j) = nansum(nansum( gpp(biome_half==2).*area(biome_half==2) )) * scale; % TgC day-1
    end
end
GPP_extratropical_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_extratropical_annual(i) = nansum(nansum( gpp(biome_half==2).*area(biome_half==2) )) * scale; % TgC yr-1
end

% Tundra & Arctic shrubland
GPP_tundra_monthly = NaN(size(GPP_annual, 3), 12);
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_tundra_monthly(i,j) = nansum(nansum( gpp(biome_half==3).*area(biome_half==3) )) * scale; % TgC day-1
    end
end
GPP_tundra_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_tundra_annual(i) = nansum(nansum( gpp(biome_half==3).*area(biome_half==3) )) * scale; % TgC yr-1
end

% Grass/crop
GPP_grass_monthly = NaN(size(GPP_annual, 3), 12);
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_grass_monthly(i,j) = nansum(nansum( gpp(biome_half==4).*area(biome_half==4) )) * scale; % TgC day-1
    end
end
GPP_grass_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_grass_annual(i) = nansum(nansum( gpp(biome_half==4).*area(biome_half==4) )) * scale; % TgC yr-1
end

% Semiarid
GPP_semiarid_monthly = NaN(size(GPP_annual, 3), 12);
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        
        gpp = GPP(:, :, yr==yrs(i) & mo==j);
        GPP_semiarid_monthly(i,j) = nansum(nansum( gpp(biome_half==5).*area(biome_half==5) )) * scale; % TgC day-1
    end
end
GPP_semiarid_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    gpp = GPP_annual(:,:,i);
    GPP_semiarid_annual(i) = nansum(nansum( gpp(biome_half==5).*area(biome_half==5) )) * scale; % TgC yr-1
end


clear i j k gpp yrs area e eyear syear R nt nx ny scale GPP GPP_monthly GPP_annual GPP_global* latidx lonidx lat lon yr mo rlim biome*;

save('./data/gpp_ccw_ahlstrom.mat');


