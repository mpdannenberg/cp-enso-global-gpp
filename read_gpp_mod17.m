% Get 1982-2015 GPP data from MOD17 and format for comparison to other data
% sources (MsTMIP, CCW, and FluxCom)

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
cd('D:\Data_Analysis\GIMMS3g_GPP');
info = ncinfo('gpp_CRUNCEP_V4P1_Standard_1982_Monthly_GEO_30min.nc');
lat = ncread('gpp_CRUNCEP_V4P1_Standard_1982_Monthly_GEO_30min.nc', 'lat');
lon = ncread('gpp_CRUNCEP_V4P1_Standard_1982_Monthly_GEO_30min.nc', 'lon');
ny = length(lat); nx = length(lon);

yr = reshape(repmat(syear:eyear, 12, 1), [], 1);
nt = length(yr);
mo = repmat(1:12, 1, length(syear:eyear))';
ndys = repmat([31,28,31,30,31,30,31,31,30,31,30,31]', 1, ny, nx); % Number of days, excluding leap years
ndys = permute(ndys, [2 3 1]);

GPP = NaN(ny, nx, nt);

for i = 1:length(years)
    
    fn = ['gpp_CRUNCEP_V4P1_Standard_',num2str(years(i)),'_Monthly_GEO_30min.nc'];
    gpp = ncread(fn, 'GPP');
    gpp = 0.001 * permute(gpp, [2 1 3]) ./ ndys; % from gC m-2 month-1 --> kgC m-2 day-1
    
    GPP(:, :, yr==years(i)) = gpp;
    
end
GPP(repmat(biomes',1,1,nt)==0) = NaN;
clear i gpp fn;
cd('D:\Publications\Dannenberg_et_al_CPElNinoGlobalGPP');

%% Area (in m^2) of the 1/2 deg grid
[LON, LAT] = meshgrid(lon, lat);
area = areaquad(reshape(LAT-(1/4),[],1),reshape(LON-(1/4),[],1),reshape(LAT+(1/4),[],1),reshape(LON+(1/4),[],1),e);
area = reshape(area, length(lat), length(lon)); 
clear LON LAT;

ndys = repmat([31,28,31,30,31,30,31,31,30,31,30,31]', length(years), ny, nx); % Number of days, excluding leap years
ndys = permute(ndys, [2 3 1]);

%% Aggregate to monthly and annual scales
% Monthly gridded mean
GPP_monthly = GPP .* ndys; % daily average --> monthly total GPP (kgC m-2 month-1)
GPP_monthly(isnan(GPP_monthly)) = 0; % For purposes of annual sums (i.e., ignore NaN in filter below)

% Multimodel annual gridded mean
windowSize = 12;
b = ones(1,windowSize);
a = 1;
GPP_annual = filter(b, a, GPP_monthly, [], 3); % 12-month running sums (kgC m-2 yr-1)
GPP_annual = GPP_annual(:, :, mo==12); % Get calendar year sum
GPP_annual(repmat(biomes',1,1,length(years))==0) = NaN; % Return ocean values to NaN

clear GPP_monthly a b windowSize ndys biomes;

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

clear i j k gpp yrs area e eyear syear R nt nx ny scale info;

GPP_monthly = GPP; clear GPP; % Monthly GPP in kgC m-2 day-1

save('./data/gpp_mod17.mat');

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
rlim = [-25 10; -80 -35];
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
rlim = [-35 5; 8 42];
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

% Tropical Asia
rlim = [-10 25; 80 150];
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

save('./data/gpp_mod17_regional.mat');



