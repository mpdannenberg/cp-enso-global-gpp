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
mo = repmat(1:12, 1, length(syear:eyear))';
ndys = repmat([31,28,31,30,31,30,31,31,30,31,30,31], 1, length(syear:eyear))'; % Number of days, excluding leap years

GPP = NaN(ny, nx, nt, length(models));
for i = 1:length(models)
    
    gpp = ncread([models{i},'_SG1_Monthly_GPP.nc4'],'GPP') * 60 * 60 * 24; % from kgC m-2 s-1 --> kgC m-2 day-1
    gpp(gpp==-9999) = NaN;
    
    GPP(:, :, :, i) = permute(gpp(:, :, idx), [2 1 3]);
    
end
clear i gpp idx;
cd('C:\Users\dannenberg\Documents\Publications\Dannenberg_et_al_CPElNinoGlobalGPP\cp-enso-global-gpp');

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

yrs = syear:eyear;
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
years = yrs;

clear i j k gpp yrs area GPP GPP_annual mo nt nx ny scale syear yr eyear;

save('./data/gpp_mstmip.mat');

