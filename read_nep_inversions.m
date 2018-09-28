% Get 1979-2013 NEP data from CO2 inversions and format for comparison to other data
% sources (MOD17, CCW, and FluxCom)

%% Units for final NEP data (original: kgC m-2 s-1)
% monthly gridded NEP: kgC m-2 day-1
% monthly global NEP: TgC day-1
% annual gridded NEP: kgC m-2 year-1
% annual global NEP: TgC year-1

%% Setup
syear = 1979; % First year of analysis
eyear = 2013; % Last year of analysis
scale = 10^-9; % kg --> Tg
models = {'v15r2','v15r4','v13r1','v14r2'};

%% Load MsTMIP NEP data
years = syear:eyear;
cd('C:\Users\dannenberg\Documents\Data_Analysis\CAMS_Inversion');
info = ncinfo('z_cams_l_lsce_197901_v15r4_ra_sfc_mm_co2flux.nc');
lat = ncread('z_cams_l_lsce_197901_v15r4_ra_sfc_mm_co2flux.nc', 'latitude');
lon = ncread('z_cams_l_lsce_197901_v15r4_ra_sfc_mm_co2flux.nc', 'longitude');
area = ncread('z_cams_l_lsce_197901_v15r4_ra_sfc_mm_co2flux.nc', 'area');
lsf = ncread('z_cams_l_lsce_197901_v15r4_ra_sfc_mm_co2flux.nc', 'lsf');
ny = length(lat); nx = length(lon);

yr = reshape(repmat(syear:eyear, 12, 1), [], 1);
nt = length(yr);
mo = repmat(1:12, 1, length(syear:eyear))';
ndys = repmat([31,28,31,30,31,30,31,31,30,31,30,31], 1, length(syear:eyear))'; % Number of days, excluding leap years

NEP = NaN(ny, nx, nt, length(models));
for i = 1:nt
    
    % CAMS 
    cd('C:\Users\dannenberg\Documents\Data_Analysis\CAMS_Inversion');
    
    fn = ['z_cams_l_lsce_',sprintf('%04d%02d',yr(i),mo(i)),'_',models{1},'_ra_sfc_mm_co2flux.nc'];
    nep = ncread(fn, 'flux_apos_bio');
    nep = -1 * nep ./ ndys(i); % from kgC m-2 month-1 --> kgC m-2 day-1 _AND_ from NEE to NEP
    NEP(:, :, i, 1) = nep;
    
    fn = ['z_cams_l_lsce_',sprintf('%04d%02d',yr(i),mo(i)),'_',models{2},'_ra_sfc_mm_co2flux.nc'];
    nep = ncread(fn, 'flux_apos_bio');
    nep = -1 * nep ./ ndys(i); % from kgC m-2 month-1 --> kgC m-2 day-1 _AND_ from NEE to NEP
    NEP(:, :, i, 2) = nep;
    
    % MACC-III 
    cd('C:\Users\dannenberg\Documents\Data_Analysis\MACC-III_Inversion');
    
    fn = ['z_macc_l_lsce_',sprintf('%04d%02d',yr(i),mo(i)),'_',models{3},'_ra_sfc_mm_co2flux.nc'];
    nep = ncread(fn, 'flux_apos_bio');
    nep = -1 * nep ./ ndys(i); % from kgC m-2 month-1 --> kgC m-2 day-1 _AND_ from NEE to NEP
    NEP(:, :, i, 3) = nep;
    
    fn = ['z_macc_l_lsce_',sprintf('%04d%02d',yr(i),mo(i)),'_',models{4},'_ra_sfc_mm_co2flux.nc'];
    nep = ncread(fn, 'flux_apos_bio');
    nep = -1 * nep ./ ndys(i); % from kgC m-2 month-1 --> kgC m-2 day-1 _AND_ from NEE to NEP
    NEP(:, :, i, 4) = nep;
    
end
clear i nep idx fn;
cd('C:\Users\dannenberg\Documents\Publications\Dannenberg_et_al_CPElNinoGlobalGPP\cp-enso-global-gpp');

%% Aggregate to monthly and annual scales
% Multimodel monthly gridded mean
NEP_monthly = bsxfun(@times,NEP,reshape(ndys,1,1,[])); % daily average --> monthly total NEP

% Multimodel annual gridded mean
windowSize = 12;
b = ones(1,windowSize);
a = 1;
NEP_annual = filter(b, a, NEP_monthly, [], 3); % 12-month running sums (kgC m-2 yr-1)
NEP_annual = NEP_annual(:, :, mo==12, :); % Get calendar year sum
NEP_annual_mean = nanmean(NEP_annual, 4);
clear NEP_monthly a b windowSize ndys;

%% Calculate global NEP at monthly and annual scale
% [LON, LAT] = meshgrid(lon, lat);
% e = referenceEllipsoid('World Geodetic System 1984');
% area = areaquad(reshape(LAT-0.25,[],1),reshape(LON-0.25,[],1),reshape(LAT+0.25,[],1),reshape(LON+0.25,[],1),e);
% area = reshape(area, ny, nx); 
% clear LON LAT e;

yrs = syear:eyear;
NEP_global_monthly = NaN(size(NEP_annual, 3), 12, length(models));
for i = 1:size(NEP_annual, 3)
    for j = 1:12
        for k = 1:length(models)
            nep = NEP(:, :, yr==yrs(i) & mo==j, k);
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
years = yrs;

clear i j k nep yrs area NEP NEP_annual mo nt nx ny scale syear yr eyear info;

%% Detrend (strong long-term increases in NEP)

mdl = fitlm(years', NEP_global_annual_mean);
NEP_global_annual_mean = mdl.Residuals.Raw;

for k=1:4
    mdl = fitlm(years', NEP_global_annual(:,k));
    NEP_global_annual(:,k) = mdl.Residuals.Raw;
end

for i = 1:12
    mdl = fitlm(years', NEP_global_monthly_mean(:,i));
    NEP_global_monthly_mean(:,i) = mdl.Residuals.Raw;
    
    for k =1:4
        
        mdl = fitlm(years', NEP_global_monthly(:,i,k));
        NEP_global_monthly(:,i,k) = mdl.Residuals.Raw;
        
    end
end
clear i k mdl;

save('./data/nep_inversions.mat');

