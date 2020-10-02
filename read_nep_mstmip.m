% Get 1950-2010 NEP data from MsTMIP and format for comparison to other data
% sources (MOD17, CCW, and FluxCom)

%% Units for final NEP data (original: kgC m-2 s-1)
% monthly gridded NEP: kgC m-2 day-1
% monthly global NEP: TgC day-1
% annual gridded NEP: kgC m-2 year-1
% annual global NEP: TgC year-1

%% Setup
syear = 1951; % First year of analysis
eyear = 2010; % End year of analysis
scale = 10^-9; % kg --> Tg
models = {'BIOME-BGC','CLM4','CLM4VIC','DLEM','GTEC',...
    'ISAM','LPJ-wsl','ORCHIDEE-LSCE','SiB3','SiBCASA','TEM6','VEGAS2.1',...
    'VISIT'};

%% Load MsTMIP NEP data
cd('D:\Data_Analysis\MsTMIP');
lat = ncread('BIOME-BGC_SG1_Monthly_NEE.nc4','lat');
lon = ncread('BIOME-BGC_SG1_Monthly_NEE.nc4','lon');
ny = length(lat); nx = length(lon);

yr = reshape(repmat(1901:2010, 12, 1), [], 1);
idx = yr >= syear;

yr = yr(idx); nt = length(yr);
mo = repmat(1:12, 1, length(syear:eyear))';
ndys = repmat([31,28,31,30,31,30,31,31,30,31,30,31], 1, length(syear:eyear))'; % Number of days, excluding leap years

NEP = NaN(ny, nx, nt, length(models));
for i = 1:length(models)
    
    nep = ncread([models{i},'_SG1_Monthly_NEE.nc4'],'NEE') * 60 * 60 * 24; % from kgC m-2 s-1 --> kgC m-2 day-1
    nep(nep==-9999) = NaN;
    
    NEP(:, :, :, i) = -1 * permute(nep(:, :, idx), [2 1 3]); % NEE-->NEP and arrange dimensions
    
end
clear i nep idx;
cd('D:\Publications\Dannenberg_et_al_CPElNinoGlobalGPP');

%% Aggregate to monthly and annual scales
% Multimodel monthly gridded mean
NEP_monthly = bsxfun(@times,NEP,reshape(ndys,1,1,[])); % daily average --> monthly total NEP
save('./data/nep_mstmip_full.mat', 'NEP','NEP_monthly', '-v7.3');
clear NEP NEP_monthly;
NEP = matfile('./data/nep_mstmip_full.mat');

% Multimodel annual gridded mean
windowSize = 12;
b = ones(1,windowSize);
a = 1;
NEP_annual = filter(b, a, NEP.NEP_monthly, [], 3); % 12-month running sums (kgC m-2 yr-1)
NEP_annual(:,:,1:(windowSize-1),:) = NaN;
NEP_shyear = NEP_annual(:, :, mo==6, :); % Get July-June sum
NEP_annual = NEP_annual(:, :, mo==12, :); % Get calendar year sum
NEP_annual_mean = nanmean(NEP_annual, 4);
NEP_shyear_mean = nanmean(NEP_shyear, 4);
clear NEP_monthly a b windowSize ndys;

%% Calculate global NEP at monthly and annual scale
[LON, LAT] = meshgrid(lon, lat);
e = referenceEllipsoid('World Geodetic System 1984');
area = areaquad(reshape(LAT-0.25,[],1),reshape(LON-0.25,[],1),reshape(LAT+0.25,[],1),reshape(LON+0.25,[],1),e);
area = reshape(area, ny, nx); 
clear LON LAT e;

yrs = syear:eyear;
NEP_global_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_global_monthly(i,j,k) = nansum(nansum( nep.*area )) * scale; % TgC day-1
        end
    end
end
NEP_global_monthly_mean = nanmean(NEP_global_monthly, 3);

NEP_global_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_global_annual(i, k) = nansum(nansum( nep.*area )) * scale; % TgC yr-1
    end
end
NEP_global_annual_mean = nanmean(NEP_global_annual, 2);

NEP_global_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_shyear(:,:,i,k);
        NEP_global_shyear(i, k) = nansum(nansum( nep.*area )) * scale; % TgC yr-1
    end
end
NEP_global_shyear_mean = nanmean(NEP_global_shyear, 2);
NEP_global_shyear_mean(1) = NaN;
NEP_global_shyear(1,:) = NaN;

years = yrs;

clear i j k nep yrs area nt nx ny scale syear eyear info;

save('./data/nep_mstmip.mat', 'NEP_annual_mean','NEP_global_annual',...
    'NEP_global_annual_mean','NEP_global_monthly','NEP_global_monthly_mean',...
    'NEP_global_shyear_mean','NEP_global_shyear','lat','lon','models','years');

%% Calculate regional NEP at monthly and annual scale

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
NEP_amazon_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_amazon_monthly(i,j,k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
NEP_amazon_monthly_mean = nanmean(NEP_amazon_monthly, 3);
NEP_amazon_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_amazon_annual(i, k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_amazon_annual_mean = nanmean(NEP_amazon_annual, 2);
NEP_amazon_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = NEP_shyear(:,:,i,k);
        NEP_amazon_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_amazon_shyear_mean = nanmean(NEP_amazon_shyear, 2);
NEP_amazon_shyear_mean(1) = NaN;
NEP_amazon_shyear(1,:) = NaN;

% Sahel
rlim = [5 15; -20 50];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
NEP_sahel_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_sahel_monthly(i,j,k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
NEP_sahel_monthly_mean = nanmean(NEP_sahel_monthly, 3);
NEP_sahel_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_sahel_annual(i, k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_sahel_annual_mean = nanmean(NEP_sahel_annual, 2);
NEP_sahel_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = NEP_shyear(:,:,i,k);
        NEP_sahel_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_sahel_shyear_mean = nanmean(NEP_sahel_shyear, 2);
NEP_sahel_shyear_mean(1) = NaN;
NEP_sahel_shyear(1,:) = NaN;

% Tropical & southern Africa
rlim = [-35 5; 8 42];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
NEP_africa_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_africa_monthly(i,j,k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
NEP_africa_monthly_mean = nanmean(NEP_africa_monthly, 3);
NEP_africa_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_africa_annual(i, k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_africa_annual_mean = nanmean(NEP_africa_annual, 2);
NEP_africa_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = NEP_shyear(:,:,i,k);
        NEP_africa_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_africa_shyear_mean = nanmean(NEP_africa_shyear, 2);
NEP_africa_shyear_mean(1) = NaN;
NEP_africa_shyear(1,:) = NaN;

% Australia
rlim = [-40 -10; 110 155];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
NEP_austr_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_austr_monthly(i,j,k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
NEP_austr_monthly_mean = nanmean(NEP_austr_monthly, 3);
NEP_austr_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_austr_annual(i, k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_austr_annual_mean = nanmean(NEP_austr_annual, 2);
NEP_austr_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = NEP_shyear(:,:,i,k);
        NEP_austr_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_austr_shyear_mean = nanmean(NEP_austr_shyear, 2);
NEP_austr_shyear_mean(1) = NaN;
NEP_austr_shyear(1,:) = NaN;

% Western North America
rlim = [20 70; -165 -100];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
NEP_westna_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_westna_monthly(i,j,k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
NEP_westna_monthly_mean = nanmean(NEP_westna_monthly, 3);
NEP_westna_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_westna_annual(i, k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_westna_annual_mean = nanmean(NEP_westna_annual, 2);

% Eastern U.S.
rlim = [25 50; -100 -60];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
NEP_eastus_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_eastus_monthly(i,j,k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
NEP_eastus_monthly_mean = nanmean(NEP_eastus_monthly, 3);
NEP_eastus_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_eastus_annual(i, k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_eastus_annual_mean = nanmean(NEP_eastus_annual, 2);
NEP_westna_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = NEP_shyear(:,:,i,k);
        NEP_westna_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_westna_shyear_mean = nanmean(NEP_westna_shyear, 2);
NEP_westna_shyear_mean(1) = NaN;
NEP_westna_shyear(1,:) = NaN;

% Europe
rlim = [35 60; -10 40];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
NEP_europe_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_europe_monthly(i,j,k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
NEP_europe_monthly_mean = nanmean(NEP_europe_monthly, 3);
NEP_europe_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_europe_annual(i, k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_europe_annual_mean = nanmean(NEP_europe_annual, 2);
NEP_europe_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = NEP_shyear(:,:,i,k);
        NEP_europe_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_europe_shyear_mean = nanmean(NEP_europe_shyear, 2);
NEP_europe_shyear_mean(1) = NaN;
NEP_europe_shyear(1,:) = NaN;

% Tropical Asia
rlim = [-10 25; 80 150];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
NEP_casia_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_casia_monthly(i,j,k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
NEP_casia_monthly_mean = nanmean(NEP_casia_monthly, 3);
NEP_casia_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_casia_annual(i, k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_casia_annual_mean = nanmean(NEP_casia_annual, 2);
NEP_casia_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = NEP_shyear(:,:,i,k);
        NEP_casia_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_casia_shyear_mean = nanmean(NEP_casia_shyear, 2);
NEP_casia_shyear_mean(1) = NaN;
NEP_casia_shyear(1,:) = NaN;

% Tropical BAND
rlim = [-23.5 23.5; -180 180];
latidx = lat>=min(rlim(1,:)) & lat<=max(rlim(1,:));
lonidx = lon>=min(rlim(2,:)) & lon<=max(rlim(2,:));
NEP_tropics_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_tropics_monthly(i,j,k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC day-1
        end
    end
end
NEP_tropics_monthly_mean = nanmean(NEP_tropics_monthly, 3);
NEP_tropics_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_tropics_annual(i, k) = nansum(nansum( nep(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_tropics_annual_mean = nanmean(NEP_tropics_annual, 2);
NEP_tropics_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = NEP_shyear(:,:,i,k);
        NEP_tropics_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_tropics_shyear_mean = nanmean(NEP_tropics_shyear, 2);
NEP_tropics_shyear_mean(1) = NaN;
NEP_tropics_shyear(1,:) = NaN;

%% Calculate regional NEP (Ahlstrom et al. 2015 version) at monthly and annual scale

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
NEP_tropical_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_tropical_monthly(i,j,k) = nansum(nansum( nep(biome_half==1).*area(biome_half==1) )) * scale; % TgC day-1
        end
    end
end
NEP_tropical_monthly_mean = nanmean(NEP_tropical_monthly, 3);
NEP_tropical_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_tropical_annual(i, k) = nansum(nansum( nep(biome_half==1).*area(biome_half==1) )) * scale; % TgC yr-1
    end
end
NEP_tropical_annual_mean = nanmean(NEP_tropical_annual, 2);
NEP_tropical_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = NEP_shyear(:,:,i,k);
        NEP_tropical_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_tropical_shyear_mean = nanmean(NEP_tropical_shyear, 2);
NEP_tropical_shyear_mean(1) = NaN;
NEP_tropical_shyear(1,:) = NaN;

% Extratropical forest
NEP_extratropical_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_extratropical_monthly(i,j,k) = nansum(nansum( nep(biome_half==2).*area(biome_half==2) )) * scale; % TgC day-1
        end
    end
end
NEP_extratropical_monthly_mean = nanmean(NEP_extratropical_monthly, 3);
NEP_extratropical_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_extratropical_annual(i, k) = nansum(nansum( nep(biome_half==2).*area(biome_half==2) )) * scale; % TgC yr-1
    end
end
NEP_extratropical_annual_mean = nanmean(NEP_extratropical_annual, 2);
NEP_extratropical_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = NEP_shyear(:,:,i,k);
        NEP_extratropical_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_extratropical_shyear_mean = nanmean(NEP_extratropical_shyear, 2);
NEP_extratropical_shyear_mean(1) = NaN;
NEP_extratropical_shyear(1,:) = NaN;

% Tundra/Arctic shrubland
NEP_tundra_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_tundra_monthly(i,j,k) = nansum(nansum( nep(biome_half==3).*area(biome_half==3) )) * scale; % TgC day-1
        end
    end
end
NEP_tundra_monthly_mean = nanmean(NEP_tundra_monthly, 3);
NEP_tundra_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_tundra_annual(i, k) = nansum(nansum( nep(biome_half==3).*area(biome_half==3) )) * scale; % TgC yr-1
    end
end
NEP_tundra_annual_mean = nanmean(NEP_tundra_annual, 2);
NEP_tundra_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = NEP_shyear(:,:,i,k);
        NEP_tundra_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_tundra_shyear_mean = nanmean(NEP_tundra_shyear, 2);
NEP_tundra_shyear_mean(1) = NaN;
NEP_tundra_shyear(1,:) = NaN;

% Grass/crop
NEP_grass_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_grass_monthly(i,j,k) = nansum(nansum( nep(biome_half==4).*area(biome_half==4) )) * scale; % TgC day-1
        end
    end
end
NEP_grass_monthly_mean = nanmean(NEP_grass_monthly, 3);
NEP_grass_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_grass_annual(i, k) = nansum(nansum( nep(biome_half==4).*area(biome_half==4) )) * scale; % TgC yr-1
    end
end
NEP_grass_annual_mean = nanmean(NEP_grass_annual, 2);
NEP_grass_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = NEP_shyear(:,:,i,k);
        NEP_grass_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_grass_shyear_mean = nanmean(NEP_grass_shyear, 2);
NEP_grass_shyear_mean(1) = NaN;
NEP_grass_shyear(1,:) = NaN;

% Semiarid
NEP_semiarid_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP.NEP(:, :, find(yr==yrs(i) & mo==j), k);
            NEP_semiarid_monthly(i,j,k) = nansum(nansum( nep(biome_half==5).*area(biome_half==5) )) * scale; % TgC day-1
        end
    end
end
NEP_semiarid_monthly_mean = nanmean(NEP_semiarid_monthly, 3);
NEP_semiarid_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        nep = NEP_annual(:,:,i,k);
        NEP_semiarid_annual(i, k) = nansum(nansum( nep(biome_half==5).*area(biome_half==5) )) * scale; % TgC yr-1
    end
end
NEP_semiarid_annual_mean = nanmean(NEP_semiarid_annual, 2);
NEP_semiarid_shyear = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = NEP_shyear(:,:,i,k);
        NEP_semiarid_shyear(i, k) = nansum(nansum( gpp(latidx, lonidx).*area(latidx, lonidx) )) * scale; % TgC yr-1
    end
end
NEP_semiarid_shyear_mean = nanmean(NEP_semiarid_shyear, 2);
NEP_semiarid_shyear_mean(1) = NaN;
NEP_semiarid_shyear(1,:) = NaN;


clear i j k nep yrs area e eyear syear R nt nx ny scale biome* NEP NEP_monthly NEP_annual NEP_global* latidx lonidx lat lon yr mo rlim idx;

%% 6 month lags for monthly analyses
NEPlag = lagmatrix(NEP_africa_monthly_mean, 1); NEP_africa_monthly_mean = [NEPlag(:,7:12) NEP_africa_monthly_mean];
NEPlag = NaN(size(NEP_africa_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_africa_monthly(:,:,i), 1); end
NEP_africa_monthly = cat(2, NEPlag(:,7:12,:), NEP_africa_monthly);

NEPlag = lagmatrix(NEP_amazon_monthly_mean, 1); NEP_amazon_monthly_mean = [NEPlag(:,7:12) NEP_amazon_monthly_mean];
NEPlag = NaN(size(NEP_amazon_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_amazon_monthly(:,:,i), 1); end
NEP_amazon_monthly = cat(2, NEPlag(:,7:12,:), NEP_amazon_monthly);

NEPlag = lagmatrix(NEP_austr_monthly_mean, 1); NEP_austr_monthly_mean = [NEPlag(:,7:12) NEP_austr_monthly_mean];
NEPlag = NaN(size(NEP_austr_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_austr_monthly(:,:,i), 1); end
NEP_austr_monthly = cat(2, NEPlag(:,7:12,:), NEP_austr_monthly);

NEPlag = lagmatrix(NEP_casia_monthly_mean, 1); NEP_casia_monthly_mean = [NEPlag(:,7:12) NEP_casia_monthly_mean];
NEPlag = NaN(size(NEP_casia_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_casia_monthly(:,:,i), 1); end
NEP_casia_monthly = cat(2, NEPlag(:,7:12,:), NEP_casia_monthly);

NEPlag = lagmatrix(NEP_eastus_monthly_mean, 1); NEP_eastus_monthly_mean = [NEPlag(:,7:12) NEP_eastus_monthly_mean];
NEPlag = NaN(size(NEP_eastus_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_eastus_monthly(:,:,i), 1); end
NEP_eastus_monthly = cat(2, NEPlag(:,7:12,:), NEP_eastus_monthly);

NEPlag = lagmatrix(NEP_europe_monthly_mean, 1); NEP_europe_monthly_mean = [NEPlag(:,7:12) NEP_europe_monthly_mean];
NEPlag = NaN(size(NEP_europe_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_europe_monthly(:,:,i), 1); end
NEP_europe_monthly = cat(2, NEPlag(:,7:12,:), NEP_europe_monthly);

NEPlag = lagmatrix(NEP_extratropical_monthly_mean, 1); NEP_extratropical_monthly_mean = [NEPlag(:,7:12) NEP_extratropical_monthly_mean];
NEPlag = NaN(size(NEP_extratropical_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_extratropical_monthly(:,:,i), 1); end
NEP_extratropical_monthly = cat(2, NEPlag(:,7:12,:), NEP_extratropical_monthly);

NEPlag = lagmatrix(NEP_grass_monthly_mean, 1); NEP_grass_monthly_mean = [NEPlag(:,7:12) NEP_grass_monthly_mean];
NEPlag = NaN(size(NEP_grass_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_grass_monthly(:,:,i), 1); end
NEP_grass_monthly = cat(2, NEPlag(:,7:12,:), NEP_grass_monthly);

NEPlag = lagmatrix(NEP_sahel_monthly_mean, 1); NEP_sahel_monthly_mean = [NEPlag(:,7:12) NEP_sahel_monthly_mean];
NEPlag = NaN(size(NEP_sahel_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_sahel_monthly(:,:,i), 1); end
NEP_sahel_monthly = cat(2, NEPlag(:,7:12,:), NEP_sahel_monthly);

NEPlag = lagmatrix(NEP_semiarid_monthly_mean, 1); NEP_semiarid_monthly_mean = [NEPlag(:,7:12) NEP_semiarid_monthly_mean];
NEPlag = NaN(size(NEP_semiarid_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_semiarid_monthly(:,:,i), 1); end
NEP_semiarid_monthly = cat(2, NEPlag(:,7:12,:), NEP_semiarid_monthly);

NEPlag = lagmatrix(NEP_tropical_monthly_mean, 1); NEP_tropical_monthly_mean = [NEPlag(:,7:12) NEP_tropical_monthly_mean];
NEPlag = NaN(size(NEP_tropical_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_tropical_monthly(:,:,i), 1); end
NEP_tropical_monthly = cat(2, NEPlag(:,7:12,:), NEP_tropical_monthly);

NEPlag = lagmatrix(NEP_tropics_monthly_mean, 1); NEP_tropics_monthly_mean = [NEPlag(:,7:12) NEP_tropics_monthly_mean];
NEPlag = NaN(size(NEP_tropics_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_tropics_monthly(:,:,i), 1); end
NEP_tropics_monthly = cat(2, NEPlag(:,7:12,:), NEP_tropics_monthly);

NEPlag = lagmatrix(NEP_tundra_monthly_mean, 1); NEP_tundra_monthly_mean = [NEPlag(:,7:12) NEP_tundra_monthly_mean];
NEPlag = NaN(size(NEP_tundra_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_tundra_monthly(:,:,i), 1); end
NEP_tundra_monthly = cat(2, NEPlag(:,7:12,:), NEP_tundra_monthly);

NEPlag = lagmatrix(NEP_westna_monthly_mean, 1); NEP_westna_monthly_mean = [NEPlag(:,7:12) NEP_westna_monthly_mean];
NEPlag = NaN(size(NEP_westna_monthly));
for i=1:length(models); NEPlag(:,:,i) = lagmatrix(NEP_westna_monthly(:,:,i), 1); end
NEP_westna_monthly = cat(2, NEPlag(:,7:12,:), NEP_westna_monthly);

clear i NEPlag;

save('./data/nep_mstmip_regional.mat');

