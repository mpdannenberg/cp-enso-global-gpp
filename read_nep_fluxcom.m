% Get 1982-2013 NEP data from FluxCom and format for comparison to other data
% sources (MsTMIP, CCW, and MOD17)

% Use only RF model to avoid overrepresenting FluxCom in data-driven model
% ensemble? Or take mean of all three?

% Hmm... global sums seem off by a factor of 10

%% Units for final NEP data (original: g C m-2 day-1)
% monthly gridded NEP: kgC m-2 day-1
% monthly global NEP: TgC day-1
% annual gridded NEP: kgC m-2 year-1
% annual global NEP: TgC year-1

%% Setup
syear = 1982; % First year of analysis
eyear = 2013; % Last year of analysis
scale = 10^-9; % kg --> Tg
e = referenceEllipsoid('World Geodetic System 1984');

%% Load CCW NEP data - gC m-2 month-1
years = syear:eyear;
cd('C:\Users\dannenberg\Documents\Data_Analysis\FluxCom');
info = ncinfo('NEE.RF.CRUNCEP_v6.monthly.1980.nc');
lat = ncread('NEE.RF.CRUNCEP_v6.monthly.1980.nc', 'lat');
lon = ncread('NEE.RF.CRUNCEP_v6.monthly.1980.nc', 'lon');
ny = length(lat); nx = length(lon);

yr = reshape(repmat(syear:eyear, 12, 1), [], 1);
nt = length(yr);
mo = repmat(1:12, 1, length(syear:eyear))';

NEP = NaN(ny, nx, nt);

for i = 1:length(years)
    
    fn = ['NEE.RF.CRUNCEP_v6.monthly.',num2str(years(i)),'.nc'];
    nep = ncread(fn, 'NEE');
    nep(nep == -9999) = NaN;
    nep = -1 * 0.001 * permute(nep, [2 1 3]); % from gC m-2 day-1 --> kgC m-2 day-1 _AND_ from NEE to NEP
    
    NEP(:, :, yr==years(i)) = nep;
    
end
clear i nep fn;
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
NEP_monthly = NEP .* ndys; % daily average --> monthly total NEP (kgC m-2 month-1)

% Multimodel annual gridded mean
windowSize = 12;
b = ones(1,windowSize);
a = 1;
NEP_annual = filter(b, a, NEP_monthly, [], 3); % 12-month running sums (kgC m-2 yr-1)
NEP_annual = NEP_annual(:, :, mo==12); % Get calendar year sum
clear NEP_monthly a b windowSize ndys;

%% Calculate global NEP at monthly and annual scale
yrs = years;
NEP_global_monthly = NaN(size(NEP_annual, 3), 12);
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        nep = NEP(:, :, yr==yrs(i) & mo==j);
        NEP_global_monthly(i,j) = nansum(nansum( nep.*area )) * scale; % TgC day-1
    end
end

NEP_global_annual = NaN(length(yrs), 1);
for i = 1:length(yrs)
    nep = NEP_annual(:,:,i);
    NEP_global_annual(i) = nansum(nansum( nep.*area )) * scale; % TgC yr-1
end

clear i j k nep yrs area e eyear syear R nt nx ny scale info;

NEP_monthly = NEP; clear NEP; % Monthly NEP in kgC m-2 day-1

save('./data/nep_fluxcom.mat');

