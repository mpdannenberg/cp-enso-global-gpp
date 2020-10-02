% Get 1950-2010 GPP data from MsTMIP and format for comparison to other data
% sources (MOD17, CCW, and FluxCom)

%% Units for final GPP data (original: kgC m-2 s-1)
% monthly gridded GPP: kgC m-2 day-1
% monthly global GPP: TgC day-1
% annual gridded GPP: kgC m-2 year-1
% annual global GPP: TgC year-1

%% Setup
syear = 1951; % First year of analysis
eyear = 2010; % End year of analysis
scale = 10^-9; % kg --> Tg
models = {'BIOME-BGC','CLM4','CLM4VIC','DLEM','GTEC',...
    'ISAM','LPJ-wsl','ORCHIDEE-LSCE','SiB3','SiBCASA','TEM6','VEGAS2.1',...
    'VISIT'};

%% Load MsTMIP GPP data
cd('D:\Data_Analysis\MsTMIP');
lat = ncread('BIOME-BGC_SG1_Monthly_GPP.nc4','lat');
lon = ncread('BIOME-BGC_SG1_Monthly_GPP.nc4','lon');
ny = length(lat); nx = length(lon);

yr = reshape(repmat(1901:2010, 12, 1), [], 1);
idx = yr >= syear;

yr = yr(idx); nt = length(yr);
mo = repmat(1:12, 1, length(syear:eyear))';
ndys = repmat([31,28,31,30,31,30,31,31,30,31,30,31], 1, length(syear:eyear))'; % Number of days, excluding leap years

GPP = NaN(ny, nx, nt, length(models));
for i = 1:length(models)
    
    gpp = ncread([models{i},'_SG1_Monthly_GPP.nc4'],'GPP') * 60 * 60 * 24; % from kgC m-2 s-1 --> kgC m-2 day-1
    gpp(gpp==-9999) = NaN;
    
    GPP(:, :, :, i) = permute(gpp(:, :, idx), [2 1 3]);
    
end
clear i gpp idx;
cd('D:\Publications\Dannenberg_et_al_CPElNinoGlobalGPP');

%% Aggregate to monthly and annual scales
% Multimodel monthly gridded mean
GPP_monthly = bsxfun(@times,GPP,reshape(ndys,1,1,[])); % daily average --> monthly total GPP
save('./data/gpp_mstmip_full.mat', 'GPP','GPP_monthly', '-v7.3');
clear GPP GPP_monthly;
GPP = matfile('./data/gpp_mstmip_full.mat');

% Multimodel annual gridded mean
windowSize = 12;
b = ones(1,windowSize);
a = 1;
GPP_annual = filter(b, a, GPP.GPP_monthly, [], 3); % 12-month running sums (kgC m-2 yr-1)
GPP_annual(:,:,1:(windowSize-1),:) = NaN;
GPP_shyear = GPP_annual(:, :, mo==6, :); % Get July-June sum
GPP_annual = GPP_annual(:, :, mo==12, :); % Get calendar year sum
GPP_annual_mean = nanmean(GPP_annual, 4);
GPP_shyear_mean = nanmean(GPP_shyear, 4);
clear a b windowSize ndys;

%% Calculate global GPP at monthly and annual scale
[LON, LAT] = meshgrid(lon, lat);
e = referenceEllipsoid('World Geodetic System 1984');
area = areaquad(reshape(LAT-0.25,[],1),reshape(LON-0.25,[],1),reshape(LAT+0.25,[],1),reshape(LON+0.25,[],1),e);
area = reshape(area, ny, nx); 
clear LON LAT e;

yrs = syear:eyear;
GPP_global_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
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

GPP_global_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_global_shyear(i, k) = nansum(nansum( gpp.*area )) * scale; % TgC yr-1
    end
end
GPP_global_shyear_mean = nanmean(GPP_global_shyear, 2);
GPP_global_shyear_mean(1) = NaN;
GPP_global_shyear(1,:) = NaN;

years = yrs;

clear i j k gpp yrs area nt nx ny scale syear eyear;

save('./data/gpp_mstmip.mat', 'GPP_annual_mean','GPP_global_annual',...
    'GPP_global_annual_mean','GPP_global_monthly','GPP_global_monthly_mean',...
    'GPP_global_shyear','GPP_global_shyear_mean','lat','lon','models','years');

%% Calculate regional GPP at monthly and annual scale

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
GPP_amazon_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_amazon_monthly(i,j,k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
GPP_amazon_monthly_mean = nanmean(GPP_amazon_monthly, 3);
GPP_amazon_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_amazon_annual(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_amazon_annual_mean = nanmean(GPP_amazon_annual, 2);
GPP_amazon_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_amazon_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_amazon_shyear_mean = nanmean(GPP_amazon_shyear, 2);
GPP_amazon_shyear_mean(1) = NaN;
GPP_amazon_shyear(1,:) = NaN;

% Sahel
rlim = [5 15; -20 50];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
GPP_sahel_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_sahel_monthly(i,j,k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
GPP_sahel_monthly_mean = nanmean(GPP_sahel_monthly, 3);
GPP_sahel_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_sahel_annual(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_sahel_annual_mean = nanmean(GPP_sahel_annual, 2);
GPP_sahel_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_sahel_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_sahel_shyear_mean = nanmean(GPP_sahel_shyear, 2);
GPP_sahel_shyear_mean(1) = NaN;
GPP_sahel_shyear(1,:) = NaN;

% Tropical and subtropical Africa
rlim = [-35 5; 8 42];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
GPP_africa_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_africa_monthly(i,j,k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
GPP_africa_monthly_mean = nanmean(GPP_africa_monthly, 3);
GPP_africa_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_africa_annual(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_africa_annual_mean = nanmean(GPP_africa_annual, 2);
GPP_africa_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_africa_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_africa_shyear_mean = nanmean(GPP_africa_shyear, 2);
GPP_africa_shyear_mean(1) = NaN;
GPP_africa_shyear(1,:) = NaN;

% Australia
rlim = [-40 -10; 110 155];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
GPP_austr_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_austr_monthly(i,j,k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
GPP_austr_monthly_mean = nanmean(GPP_austr_monthly, 3);
GPP_austr_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_austr_annual(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_austr_annual_mean = nanmean(GPP_austr_annual, 2);
GPP_austr_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_austr_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_austr_shyear_mean = nanmean(GPP_austr_shyear, 2);
GPP_austr_shyear_mean(1) = NaN;
GPP_austr_shyear(1,:) = NaN;

% Western North America
rlim = [20 70; -165 -100];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
GPP_westna_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_westna_monthly(i,j,k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
GPP_westna_monthly_mean = nanmean(GPP_westna_monthly, 3);
GPP_westna_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_westna_annual(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_westna_annual_mean = nanmean(GPP_westna_annual, 2);
GPP_westna_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_westna_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_westna_shyear_mean = nanmean(GPP_westna_shyear, 2);
GPP_westna_shyear_mean(1) = NaN;
GPP_westna_shyear(1,:) = NaN;

% Eastern U.S.
rlim = [25 50; -100 -60];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
GPP_eastus_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_eastus_monthly(i,j,k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
GPP_eastus_monthly_mean = nanmean(GPP_eastus_monthly, 3);
GPP_eastus_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_eastus_annual(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_eastus_annual_mean = nanmean(GPP_eastus_annual, 2);
GPP_eastus_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_eastus_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_eastus_shyear_mean = nanmean(GPP_eastus_shyear, 2);
GPP_eastus_shyear_mean(1) = NaN;
GPP_eastus_shyear(1,:) = NaN;

% Europe
rlim = [35 60; -10 40];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
GPP_europe_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_europe_monthly(i,j,k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
GPP_europe_monthly_mean = nanmean(GPP_europe_monthly, 3);
GPP_europe_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_europe_annual(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_europe_annual_mean = nanmean(GPP_europe_annual, 2);
GPP_europe_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_europe_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_europe_shyear_mean = nanmean(GPP_europe_shyear, 2);
GPP_europe_shyear_mean(1) = NaN;
GPP_europe_shyear(1,:) = NaN;

% Tropical Asia
rlim = [-10 25; 80 150];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
GPP_casia_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_casia_monthly(i,j,k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
GPP_casia_monthly_mean = nanmean(GPP_casia_monthly, 3);
GPP_casia_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_casia_annual(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_casia_annual_mean = nanmean(GPP_casia_annual, 2);
GPP_casia_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_casia_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_casia_shyear_mean = nanmean(GPP_casia_shyear, 2);
GPP_casia_shyear_mean(1) = NaN;
GPP_casia_shyear(1,:) = NaN;

% Tropical BAND
rlim = [-23.5 23.5; -180 180];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
GPP_tropics_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_tropics_monthly(i,j,k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
GPP_tropics_monthly_mean = nanmean(GPP_tropics_monthly, 3);
GPP_tropics_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_tropics_annual(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_tropics_annual_mean = nanmean(GPP_tropics_annual, 2);
GPP_tropics_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_tropics_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_tropics_shyear_mean = nanmean(GPP_tropics_shyear, 2);
GPP_tropics_shyear_mean(1) = NaN;
GPP_tropics_shyear(1,:) = NaN;

%% Calculate regional GPP (Ahlstrom et al. 2015 version) at monthly and annual scale

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
GPP_tropical_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_tropical_monthly(i,j,k) = nansum(nansum( gpp(biome_half==1).*area(biome_half==1) )) * scale; % TgC day-1
        end
    end
end
GPP_tropical_monthly_mean = nanmean(GPP_tropical_monthly, 3);
GPP_tropical_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_tropical_annual(i, k) = nansum(nansum( gpp(biome_half==1).*area(biome_half==1) )) * scale; % TgC yr-1
    end
end
GPP_tropical_annual_mean = nanmean(GPP_tropical_annual, 2);
GPP_tropical_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_tropical_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_tropical_shyear_mean = nanmean(GPP_tropical_shyear, 2);
GPP_tropical_shyear_mean(1) = NaN;
GPP_tropical_shyear(1,:) = NaN;

% Extratropical forest
GPP_extratropical_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_extratropical_monthly(i,j,k) = nansum(nansum( gpp(biome_half==2).*area(biome_half==2) )) * scale; % TgC day-1
        end
    end
end
GPP_extratropical_monthly_mean = nanmean(GPP_extratropical_monthly, 3);
GPP_extratropical_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_extratropical_annual(i, k) = nansum(nansum( gpp(biome_half==2).*area(biome_half==2) )) * scale; % TgC yr-1
    end
end
GPP_extratropical_annual_mean = nanmean(GPP_extratropical_annual, 2);
GPP_extratropical_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_extratropical_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_extratropical_shyear_mean = nanmean(GPP_extratropical_shyear, 2);
GPP_extratropical_shyear_mean(1) = NaN;
GPP_extratropical_shyear(1,:) = NaN;

% Tundra/Arctic shrubland
GPP_tundra_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_tundra_monthly(i,j,k) = nansum(nansum( gpp(biome_half==3).*area(biome_half==3) )) * scale; % TgC day-1
        end
    end
end
GPP_tundra_monthly_mean = nanmean(GPP_tundra_monthly, 3);
GPP_tundra_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_tundra_annual(i, k) = nansum(nansum( gpp(biome_half==3).*area(biome_half==3) )) * scale; % TgC yr-1
    end
end
GPP_tundra_annual_mean = nanmean(GPP_tundra_annual, 2);
GPP_tundra_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_tundra_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_tundra_shyear_mean = nanmean(GPP_tundra_shyear, 2);
GPP_tundra_shyear_mean(1) = NaN;
GPP_tundra_shyear(1,:) = NaN;

% Grass/crop
GPP_grass_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_grass_monthly(i,j,k) = nansum(nansum( gpp(biome_half==4).*area(biome_half==4) )) * scale; % TgC day-1
        end
    end
end
GPP_grass_monthly_mean = nanmean(GPP_grass_monthly, 3);
GPP_grass_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_grass_annual(i, k) = nansum(nansum( gpp(biome_half==4).*area(biome_half==4) )) * scale; % TgC yr-1
    end
end
GPP_grass_annual_mean = nanmean(GPP_grass_annual, 2);
GPP_grass_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_grass_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_grass_shyear_mean = nanmean(GPP_grass_shyear, 2);
GPP_grass_shyear_mean(1) = NaN;
GPP_grass_shyear(1,:) = NaN;

% Semiarid
GPP_semiarid_monthly = NaN(size(GPP_annual, 3), 12, length(models));
for i = 1:size(GPP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            gpp = GPP.GPP(:, :, find(yr==yrs(i) & mo==j), k);
            GPP_semiarid_monthly(i,j,k) = nansum(nansum( gpp(biome_half==5).*area(biome_half==5) )) * scale; % TgC day-1
        end
    end
end
GPP_semiarid_monthly_mean = nanmean(GPP_semiarid_monthly, 3);
GPP_semiarid_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_semiarid_annual(i, k) = nansum(nansum( gpp(biome_half==5).*area(biome_half==5) )) * scale; % TgC yr-1
    end
end
GPP_semiarid_annual_mean = nanmean(GPP_semiarid_annual, 2);
GPP_semiarid_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_shyear(:,:,i,k);
        GPP_semiarid_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
GPP_semiarid_shyear_mean = nanmean(GPP_semiarid_shyear, 2);
GPP_semiarid_shyear_mean(1) = NaN;
GPP_semiarid_shyear(1,:) = NaN;


clear idx i j k gpp yrs area e eyear syear R nt nx ny scale GPP GPP_monthly GPP_annual GPP_shyear GPP_global* latidx lonidx lat lon yr mo rlim biome*;

%% 6 month lags for monthly analyses
GPPlag = lagmatrix(GPP_africa_monthly_mean, 1); GPP_africa_monthly_mean = [GPPlag(:,7:12) GPP_africa_monthly_mean];
GPPlag = NaN(size(GPP_africa_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_africa_monthly(:,:,i), 1); end
GPP_africa_monthly = cat(2, GPPlag(:,7:12,:), GPP_africa_monthly);

GPPlag = lagmatrix(GPP_amazon_monthly_mean, 1); GPP_amazon_monthly_mean = [GPPlag(:,7:12) GPP_amazon_monthly_mean];
GPPlag = NaN(size(GPP_amazon_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_amazon_monthly(:,:,i), 1); end
GPP_amazon_monthly = cat(2, GPPlag(:,7:12,:), GPP_amazon_monthly);

GPPlag = lagmatrix(GPP_austr_monthly_mean, 1); GPP_austr_monthly_mean = [GPPlag(:,7:12) GPP_austr_monthly_mean];
GPPlag = NaN(size(GPP_austr_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_austr_monthly(:,:,i), 1); end
GPP_austr_monthly = cat(2, GPPlag(:,7:12,:), GPP_austr_monthly);

GPPlag = lagmatrix(GPP_casia_monthly_mean, 1); GPP_casia_monthly_mean = [GPPlag(:,7:12) GPP_casia_monthly_mean];
GPPlag = NaN(size(GPP_casia_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_casia_monthly(:,:,i), 1); end
GPP_casia_monthly = cat(2, GPPlag(:,7:12,:), GPP_casia_monthly);

GPPlag = lagmatrix(GPP_eastus_monthly_mean, 1); GPP_eastus_monthly_mean = [GPPlag(:,7:12) GPP_eastus_monthly_mean];
GPPlag = NaN(size(GPP_eastus_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_eastus_monthly(:,:,i), 1); end
GPP_eastus_monthly = cat(2, GPPlag(:,7:12,:), GPP_eastus_monthly);

GPPlag = lagmatrix(GPP_europe_monthly_mean, 1); GPP_europe_monthly_mean = [GPPlag(:,7:12) GPP_europe_monthly_mean];
GPPlag = NaN(size(GPP_europe_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_europe_monthly(:,:,i), 1); end
GPP_europe_monthly = cat(2, GPPlag(:,7:12,:), GPP_europe_monthly);

GPPlag = lagmatrix(GPP_extratropical_monthly_mean, 1); GPP_extratropical_monthly_mean = [GPPlag(:,7:12) GPP_extratropical_monthly_mean];
GPPlag = NaN(size(GPP_extratropical_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_extratropical_monthly(:,:,i), 1); end
GPP_extratropical_monthly = cat(2, GPPlag(:,7:12,:), GPP_extratropical_monthly);

GPPlag = lagmatrix(GPP_grass_monthly_mean, 1); GPP_grass_monthly_mean = [GPPlag(:,7:12) GPP_grass_monthly_mean];
GPPlag = NaN(size(GPP_grass_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_grass_monthly(:,:,i), 1); end
GPP_grass_monthly = cat(2, GPPlag(:,7:12,:), GPP_grass_monthly);

GPPlag = lagmatrix(GPP_sahel_monthly_mean, 1); GPP_sahel_monthly_mean = [GPPlag(:,7:12) GPP_sahel_monthly_mean];
GPPlag = NaN(size(GPP_sahel_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_sahel_monthly(:,:,i), 1); end
GPP_sahel_monthly = cat(2, GPPlag(:,7:12,:), GPP_sahel_monthly);

GPPlag = lagmatrix(GPP_semiarid_monthly_mean, 1); GPP_semiarid_monthly_mean = [GPPlag(:,7:12) GPP_semiarid_monthly_mean];
GPPlag = NaN(size(GPP_semiarid_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_semiarid_monthly(:,:,i), 1); end
GPP_semiarid_monthly = cat(2, GPPlag(:,7:12,:), GPP_semiarid_monthly);

GPPlag = lagmatrix(GPP_tropical_monthly_mean, 1); GPP_tropical_monthly_mean = [GPPlag(:,7:12) GPP_tropical_monthly_mean];
GPPlag = NaN(size(GPP_tropical_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_tropical_monthly(:,:,i), 1); end
GPP_tropical_monthly = cat(2, GPPlag(:,7:12,:), GPP_tropical_monthly);

GPPlag = lagmatrix(GPP_tropics_monthly_mean, 1); GPP_tropics_monthly_mean = [GPPlag(:,7:12) GPP_tropics_monthly_mean];
GPPlag = NaN(size(GPP_tropics_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_tropics_monthly(:,:,i), 1); end
GPP_tropics_monthly = cat(2, GPPlag(:,7:12,:), GPP_tropics_monthly);

GPPlag = lagmatrix(GPP_tundra_monthly_mean, 1); GPP_tundra_monthly_mean = [GPPlag(:,7:12) GPP_tundra_monthly_mean];
GPPlag = NaN(size(GPP_tundra_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_tundra_monthly(:,:,i), 1); end
GPP_tundra_monthly = cat(2, GPPlag(:,7:12,:), GPP_tundra_monthly);

GPPlag = lagmatrix(GPP_westna_monthly_mean, 1); GPP_westna_monthly_mean = [GPPlag(:,7:12) GPP_westna_monthly_mean];
GPPlag = NaN(size(GPP_westna_monthly));
for i=1:length(models); GPPlag(:,:,i) = lagmatrix(GPP_westna_monthly(:,:,i), 1); end
GPP_westna_monthly = cat(2, GPPlag(:,7:12,:), GPP_westna_monthly);

clear i GPPlag;

save('./data/gpp_mstmip_regional.mat');


