% Create EP/CP Niño dataset from COBE-2 SST
t0 = datenum(1891, 1, 1);
latlim = [-20 20];
lonlim = [120 290];
syear = 1951;

warning('off');

windowSize = 6; % Number of months over which to average SSTs
endMonth = 2; % End month of integration period

%% Load COBE-2 SST data
info = ncinfo('./data/sst.mon.mean.nc');
lat = double(ncread('./data/sst.mon.mean.nc','lat'));
lon = double(ncread('./data/sst.mon.mean.nc','lon'));
t = ncread('./data/sst.mon.mean.nc','time');
sst = double(ncread('./data/sst.mon.mean.nc','sst'));
sst = permute(sst, [2 1 3]);
sst(sst>40) = NaN;
[yr, mo] = datevec(t0 + t); clear t0 t;
[nx, ny, ~] = size(sst);

%% Need to do monthly anomalies
ssta = NaN(size(sst));
for i = 1:12
    
    x = sst(:, :, mo==i);
    xclim = sst(:, :, mo==i & yr>=1981 & yr<=2010);
    xm = mean(xclim, 3);
    
    ssta(:,:,mo==i) = x - repmat(xm, 1,1,size(x, 3));
    
end
clear x xclim xm i;

%% Filter to 6-month mean SST
b = ones(1,windowSize)/windowSize;
a = 1;
ssta = filter(b, a, ssta, [], 3);
clear b a windowSize;

%% Include only time period of interest
idx = yr >= syear & mo == endMonth;
ssta = ssta(:, :, idx);
yr = yr(idx);
lat_weight = repmat(sqrt(cosd(lat)), 1, length(lon));
clear idx mo;

%% Regress SST data onto Nino-4 and Nino-1+2 indices to get independent CP- and EP- SST anomalies
nino4 = ssta(lat>=-5 & lat<=5, lon>=160 & lon<=210, :);
nino4 = reshape(permute(nino4, [3 1 2]), length(yr), []);
w = repmat(reshape(lat_weight(lat>=-5 & lat<=5, lon>=160 & lon<=210), 1, []), length(yr), 1);
% nino4 = mean(nino4, 2);
nino4 = sum(nino4.*w, 2) ./ sum(w, 2);

nino12 = ssta(lat>=-10 & lat<=0, lon>=270 & lon<=280, :);
nino12 = reshape(permute(nino12, [3 1 2]), length(yr), []);
idx = sum(isnan(nino12))==0;
w = repmat(reshape(lat_weight(lat>=-10 & lat<=0, lon>=270 & lon<=280), 1, []), length(yr), 1);
% nino12 = nanmean(nino12, 2);
nino12 = sum(nino12(:,idx).*w(:,idx), 2) ./ sum(w(:,idx), 2);

sst_tp = ssta(lat>=min(latlim) & lat<=max(latlim), lon>=min(lonlim) & lon<=max(lonlim), :);
[ny, nx, ~] = size(sst_tp);
sst_tp = reshape(permute(sst_tp, [3 1 2]), length(yr), []);
w_tp = reshape(lat_weight(lat>=min(latlim) & lat<=max(latlim), lon>=min(lonlim) & lon<=max(lonlim)), 1, []);
sst_idx = sum(isnan(sst_tp)) == 0;
sst_tp = sst_tp(:, sst_idx);
w_tp = w_tp(sst_idx);
sst_ep = NaN(size(sst_tp));
sst_cp = NaN(size(sst_tp));
for i = 1:size(sst_tp, 2)
    y = sst_tp(:, i);
    lm_ep = fitlm(nino4, y);
    sst_ep(:, i) = lm_ep.Residuals.Raw;
    lm_cp = fitlm(nino12, y);
    sst_cp(:, i) = lm_cp.Residuals.Raw;
end
[ep_coef, ep_pc, ep_latent] = pca(sst_ep, 'VariableWeights',w_tp);
[cp_coef, cp_pc, cp_latent] = pca(sst_cp, 'VariableWeights',w_tp);
ep_eof = NaN(size(sst_idx));
cp_eof = NaN(size(sst_idx));
ep_eof(sst_idx) = ep_coef(:, 1);
cp_eof(sst_idx) = cp_coef(:, 1);
ep_eof = reshape(ep_eof, ny, nx);
cp_eof = reshape(cp_eof, ny, nx);
ep_idx = ep_pc(:, 1);
cp_idx = cp_pc(:, 1);

