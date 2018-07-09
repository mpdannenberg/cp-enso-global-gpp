% Examine relationships between modeled MsTMIP GPP and EPI/CPI over 
% 1951-2010 period

syear = 1951; % First year of analysis

models = {'BIOME-BGC','CLASS-CTEM-N','CLM4','CLM4VIC','DLEM','GTEC',...
    'ISAM','LPJ-wsl','ORCHIDEE-LSCE','SiB3','SiBCASA','TEM6','VEGAS2.1',...
    'VISIT'};
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

%% Aggregate to monthly and annual scales
% Multimodel monthly gridded mean
GPP_monthly = bsxfun(@times,GPP,reshape(ndys,1,1,[])); % daily average --> monthly total GPP
GPP_monthly_mean = nanmean(GPP_monthly, 4);

% Multimodel annual gridded mean
windowSize = 12;
b = ones(1,windowSize);
a = 1;
GPP_annual = filter(b, a, GPP_monthly, [], 3); % 12-month running sums (kgC m-2 yr-1)
GPP_annual = GPP_annual(:, :, mo==12, :); % Get calendar year sum
GPP_annual_mean = nanmean(GPP_annual, 4);

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
            GPP_global_monthly(i,j,k) = nansum(nansum( gpp.*area )); % kgC day-1
        end
    end
end
GPP_global_monthly_mean = nanmean(GPP_global_monthly, 3);

GPP_global_annual = NaN(length(yrs), length(models));
for i = 1:length(yrs)
    for k = 1:length(models)
        gpp = GPP_annual(:,:,i,k);
        GPP_global_annual(i, k) = nansum(nansum( gpp.*area )); % kgC yr-1
    end
end
GPP_global_annual_mean = nanmean(GPP_global_annual, 2);


%% Regress monthly & annual gpp against CPI and EPI indices
cd('C:\Users\dannenberg\Documents\Publications\Dannenberg_et_al_CPElNinoGlobalGPP\cp-enso-global-gpp');
load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr<=2010);
epi = epi(yr<=2010);
yr = yr(yr<=2010);



