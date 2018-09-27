% Get 1982-2015 GPP data from MOD17 and format for comparison to other data
% sources (MsTMIP, CCW, and FluxCom)

%% Units for final GPP data (original: g C m-2 month-1)
% monthly gridded GPP: kgC m-2 day-1
% monthly global GPP: TgC day-1
% annual gridded GPP: kgC m-2 year-1
% annual global GPP: TgC year-1

%% Setup
syear = 1982; % First year of analysis
eyear = 2015; % Last year of analysis
scale = 10^-9; % kg --> Tg
e = referenceEllipsoid('World Geodetic System 1984');
biomes = ncread('C:\Users\dannenberg\Documents\Data_Analysis\MsTMIP\mstmip_driver_global_hd_biome_v1.nc4', 'biome_type');

%% Load CCW GPP data - gC m-2 month-1
years = syear:eyear;
cd('C:\Users\dannenberg\Documents\Data_Analysis\GIMMS3g_GPP');
info = ncinfo('gpp_V4_Standard_1982_Monthly_GEO_30min.nc');
lat = ncread('gpp_V4_Standard_1982_Monthly_GEO_30min.nc', 'lat');
lon = ncread('gpp_V4_Standard_1982_Monthly_GEO_30min.nc', 'lon');
ny = length(lat); nx = length(lon);

yr = reshape(repmat(syear:eyear, 12, 1), [], 1);
nt = length(yr);
mo = repmat(1:12, 1, length(syear:eyear))';
ndys = repmat([31,28,31,30,31,30,31,31,30,31,30,31]', 1, ny, nx); % Number of days, excluding leap years
ndys = permute(ndys, [2 3 1]);

GPP = NaN(ny, nx, nt);

for i = 1:length(years)
    
    fn = ['gpp_V4_Standard_',num2str(years(i)),'_Monthly_GEO_30min.nc'];
    gpp = ncread(fn, 'GPP');
    gpp = 0.001 * permute(gpp, [2 1 3]) ./ ndys; % from gC m-2 month-1 --> kgC m-2 day-1
    
    GPP(:, :, yr==years(i)) = gpp;
    
end
GPP(repmat(biomes',1,1,nt)==0) = NaN;
clear i gpp fn;
cd('C:\Users\dannenberg\Documents\Publications\Dannenberg_et_al_CPElNinoGlobalGPP\cp-enso-global-gpp');

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

