% Get 1950-2010 Rh data from MsTMIP and format for comparison to other data
% sources (MOD17, CCW, and FluxCom)

%% Units for final Rh data (original: kgC m-2 s-1)
% monthly gridded Rh: kgC m-2 day-1
% monthly global Rh: TgC day-1
% annual gridded Rh: kgC m-2 year-1
% annual global Rh: TgC year-1

%% Setup
syear = 1951; % First year of analysis
eyear = 2010; % End year of analysis
scale = 10^-9; % kg --> Tg
models = {'BIOME-BGC','CLM4','CLM4VIC','DLEM','GTEC',...
    'ISAM','LPJ-wsl','ORCHIDEE-LSCE','SiB3','SiBCASA','TEM6','VEGAS2.1',...
    'VISIT'};

%% Load MsTMIP Rh data
cd('D:\Data_Analysis\MsTMIP');
lat = ncread('BIOME-BGC_SG1_Monthly_HeteroResp.nc4','lat');
lon = ncread('BIOME-BGC_SG1_Monthly_HeteroResp.nc4','lon');
ny = length(lat); nx = length(lon);

yr = reshape(repmat(1901:2010, 12, 1), [], 1);
idx = yr >= syear;

yr = yr(idx); nt = length(yr);
mo = repmat(1:12, 1, length(syear:eyear))';
ndys = repmat([31,28,31,30,31,30,31,31,30,31,30,31], 1, length(syear:eyear))'; % Number of days, excluding leap years

Rh = NaN(ny, nx, nt, length(models));
for i = 1:length(models)
    
    rh = ncread([models{i},'_SG1_Monthly_HeteroResp.nc4'],'HeteroResp') * 60 * 60 * 24; % from kgC m-2 s-1 --> kgC m-2 day-1
    rh(rh==-9999) = NaN;
    
    Rh(:, :, :, i) = permute(rh(:, :, idx), [2 1 3]);
    
end
clear i rh idx;
cd('D:\Publications\Dannenberg_et_al_CPElNinoGlobalGPP');

%% Aggregate to monthly and annual scales
% Multimodel monthly gridded mean
Rh_monthly = bsxfun(@times,Rh,reshape(ndys,1,1,[])); % daily average --> monthly total Rh
save('./data/rh_mstmip_full.mat', 'Rh','Rh_monthly', '-v7.3');
clear Rh Rh_monthly;
Rh = matfile('./data/rh_mstmip_full.mat');

% Multimodel annual gridded mean
windowSize = 12;
b = ones(1,windowSize);
a = 1;
Rh_annual = filter(b, a, Rh.Rh_monthly, [], 3); % 12-month running sums (kgC m-2 yr-1)
Rh_annual = Rh_annual(:, :, mo==12, :); % Get calendar year sum
Rh_annual_mean = nanmean(Rh_annual, 4);
clear a b windowSize ndys;

%% Calculate global Rh at monthly and annual scale
[LON, LAT] = meshgrid(lon, lat);
e = referenceEllipsoid('World Geodetic System 1984');
area = areaquad(reshape(LAT-0.25,[],1),reshape(LON-0.25,[],1),reshape(LAT+0.25,[],1),reshape(LON+0.25,[],1),e);
area = reshape(area, ny, nx); 
clear LON LAT e;

yrs = syear:eyear;
Rh_global_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_global_monthly(i,j,k) = nansum(nansum( rh.*area )) * scale; % TgC day-1
        end
    end
end
Rh_global_monthly_mean = nanmean(Rh_global_monthly, 3);

Rh_global_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_global_annual(i, k) = nansum(nansum( rh.*area )) * scale; % TgC yr-1
    end
end
Rh_global_annual_mean = nanmean(Rh_global_annual, 2);
years = yrs;

clear i j k rh yrs area nt nx ny scale syear eyear;

save('./data/rh_mstmip.mat', 'Rh_annual_mean','Rh_global_annual',...
    'Rh_global_annual_mean','Rh_global_monthly','Rh_global_monthly_mean',...
    'lat','lon','models','years');

%% Calculate regional Rh at monthly and annual scale

yrs = years;
scale = 10^-9; % kg --> Tg
e = referenceEllipsoid('World Geodetic System 1984');
[LON, LAT] = meshgrid(lon, lat);
area = areaquad(reshape(LAT-(1/4),[],1),reshape(LON-(1/4),[],1),reshape(LAT+(1/4),[],1),reshape(LON+(1/4),[],1),e);
area = reshape(area, length(lat), length(lon)); 
clear LON LAT;

% Amazon
rlim = [-25 10; -80 -35];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
Rh_amazon_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_amazon_monthly(i,j,k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
Rh_amazon_monthly_mean = nanmean(Rh_amazon_monthly, 3);
Rh_amazon_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_amazon_annual(i, k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
Rh_amazon_annual_mean = nanmean(Rh_amazon_annual, 2);

% Sahel
rlim = [5 15; -20 50];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
Rh_sahel_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_sahel_monthly(i,j,k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
Rh_sahel_monthly_mean = nanmean(Rh_sahel_monthly, 3);
Rh_sahel_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_sahel_annual(i, k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
Rh_sahel_annual_mean = nanmean(Rh_sahel_annual, 2);

% Tropical and subtropical Africa
rlim = [-35 5; 8 42];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
Rh_africa_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_africa_monthly(i,j,k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
Rh_africa_monthly_mean = nanmean(Rh_africa_monthly, 3);
Rh_africa_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_africa_annual(i, k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
Rh_africa_annual_mean = nanmean(Rh_africa_annual, 2);

% Australia
rlim = [-40 -10; 110 155];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
Rh_austr_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_austr_monthly(i,j,k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
Rh_austr_monthly_mean = nanmean(Rh_austr_monthly, 3);
Rh_austr_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_austr_annual(i, k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
Rh_austr_annual_mean = nanmean(Rh_austr_annual, 2);

% Western North America
rlim = [20 70; -165 -100];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
Rh_westna_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_westna_monthly(i,j,k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
Rh_westna_monthly_mean = nanmean(Rh_westna_monthly, 3);
Rh_westna_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_westna_annual(i, k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
Rh_westna_annual_mean = nanmean(Rh_westna_annual, 2);

% Eastern U.S.
rlim = [25 50; -100 -60];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
Rh_eastus_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_eastus_monthly(i,j,k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
Rh_eastus_monthly_mean = nanmean(Rh_eastus_monthly, 3);
Rh_eastus_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_eastus_annual(i, k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
Rh_eastus_annual_mean = nanmean(Rh_eastus_annual, 2);

% Europe
rlim = [35 60; -10 40];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
Rh_europe_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_europe_monthly(i,j,k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
Rh_europe_monthly_mean = nanmean(Rh_europe_monthly, 3);
Rh_europe_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_europe_annual(i, k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
Rh_europe_annual_mean = nanmean(Rh_europe_annual, 2);

% Tropical Asia
rlim = [-10 25; 80 150];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
Rh_casia_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_casia_monthly(i,j,k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
Rh_casia_monthly_mean = nanmean(Rh_casia_monthly, 3);
Rh_casia_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_casia_annual(i, k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
Rh_casia_annual_mean = nanmean(Rh_casia_annual, 2);

% Tropical BAND
rlim = [-23.5 23.5; -180 180];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
Rh_tropics_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_tropics_monthly(i,j,k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
Rh_tropics_monthly_mean = nanmean(Rh_tropics_monthly, 3);
Rh_tropics_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_tropics_annual(i, k) = nansum(nansum( rh(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
Rh_tropics_annual_mean = nanmean(Rh_tropics_annual, 2);

%% Calculate regional Rh (Ahlstrom et al. 2015 version) at monthly and annual scale

yrs = years;
scale = 10^-9; % kg --> Tg
e = referenceEllipsoid('World Geodetic System 1984');
[LON, LAT] = meshgrid(lon, lat);
area = areaquad(reshape(LAT-(1/4),[],1),reshape(LON-(1/4),[],1),reshape(LAT+(1/4),[],1),reshape(LON+(1/4),[],1),e);
area = reshape(area, length(lat), length(lon)); 
clear LON LAT;

load ./data/ahlstrom_regions.mat;
biome_half = flipud(biome_half);

% Tropical forest
Rh_tropical_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_tropical_monthly(i,j,k) = nansum(nansum( rh(biome_half==1).*area(biome_half==1) )) * scale; % TgC day-1
        end
    end
end
Rh_tropical_monthly_mean = nanmean(Rh_tropical_monthly, 3);
Rh_tropical_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_tropical_annual(i, k) = nansum(nansum( rh(biome_half==1).*area(biome_half==1) )) * scale; % TgC yr-1
    end
end
Rh_tropical_annual_mean = nanmean(Rh_tropical_annual, 2);

% Extratropical forest
Rh_extratropical_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_extratropical_monthly(i,j,k) = nansum(nansum( rh(biome_half==2).*area(biome_half==2) )) * scale; % TgC day-1
        end
    end
end
Rh_extratropical_monthly_mean = nanmean(Rh_extratropical_monthly, 3);
Rh_extratropical_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_extratropical_annual(i, k) = nansum(nansum( rh(biome_half==2).*area(biome_half==2) )) * scale; % TgC yr-1
    end
end
Rh_extratropical_annual_mean = nanmean(Rh_extratropical_annual, 2);

% Tundra/Arctic shrubland
Rh_tundra_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_tundra_monthly(i,j,k) = nansum(nansum( rh(biome_half==3).*area(biome_half==3) )) * scale; % TgC day-1
        end
    end
end
Rh_tundra_monthly_mean = nanmean(Rh_tundra_monthly, 3);
Rh_tundra_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_tundra_annual(i, k) = nansum(nansum( rh(biome_half==3).*area(biome_half==3) )) * scale; % TgC yr-1
    end
end
Rh_tundra_annual_mean = nanmean(Rh_tundra_annual, 2);

% Grass/crop
Rh_grass_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_grass_monthly(i,j,k) = nansum(nansum( rh(biome_half==4).*area(biome_half==4) )) * scale; % TgC day-1
        end
    end
end
Rh_grass_monthly_mean = nanmean(Rh_grass_monthly, 3);
Rh_grass_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_grass_annual(i, k) = nansum(nansum( rh(biome_half==4).*area(biome_half==4) )) * scale; % TgC yr-1
    end
end
Rh_grass_annual_mean = nanmean(Rh_grass_annual, 2);

% Semiarid
Rh_semiarid_monthly = NaN(size(Rh_annual, 3), 12, length(models));
for i = 1:size(Rh_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            rh = Rh.Rh(:, :, find(yr==yrs(i) & mo==j), k);
            Rh_semiarid_monthly(i,j,k) = nansum(nansum( rh(biome_half==5).*area(biome_half==5) )) * scale; % TgC day-1
        end
    end
end
Rh_semiarid_monthly_mean = nanmean(Rh_semiarid_monthly, 3);
Rh_semiarid_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        rh = Rh_annual(:,:,i,k);
        Rh_semiarid_annual(i, k) = nansum(nansum( rh(biome_half==5).*area(biome_half==5) )) * scale; % TgC yr-1
    end
end
Rh_semiarid_annual_mean = nanmean(Rh_semiarid_annual, 2);


clear i j k rh yrs area e eyear syear R nt nx ny scale Rh Rh_monthly Rh_annual Rh_global* latidx lonidx lat lon yr mo rlim biome*;

save('./data/rh_mstmip_regional.mat');


