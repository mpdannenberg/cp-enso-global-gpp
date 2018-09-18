% Examine relationships between remotely sensed GPP and EPI/CPI over 
% 1982-2016 period

%% Units for final GPP data (original: kgC m-2 month-1)
% monthly gridded GPP: kgC m-2 day-1
% monthly global GPP: TgC day-1
% annual gridded GPP: kgC m-2 year-1
% annual global GPP: TgC year-1

%% Setup
syear = 1982; % First year of analysis
eyear = 2015; % Last year of analysis
scale = 10^-9; % kg --> Tg
e = referenceEllipsoid('World Geodetic System 1984');

%% Load CCW GPP data - gC m-2 month-1
years = syear:eyear;
cd('C:\Users\dannenberg\Documents\Data_Analysis\CCW\GeoTIFFs');
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
cd('C:\Users\dannenberg\Documents\Publications\Dannenberg_et_al_CPElNinoGlobalGPP\cp-enso-global-gpp');

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

lat = newlat; ny = length(lat);
lon = newlon; nx = length(lon);
GPP = GPP_halfdeg;
ndys = repmat([31,28,31,30,31,30,31,31,30,31,30,31]', length(years), ny, nx); % Number of days, excluding leap years
ndys = permute(ndys, [2 3 1]);

clear newlat newlon GPP_halfdeg i j area12 gpp latidx lonidx;

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

