% Examine relationships between data-driven GPP and EPI/CPI over 
% 1982-2015 period

syear = 1982;
eyear = 2016;

models = {'CCW','MOD17'}; nm = length(models);
ny = 360;
nx = 720;
nt = length(syear:eyear);

% temporary arrays for all three models
temp_annual = NaN(ny, nx, nt, nm);
temp_shyear = NaN(ny, nx, nt, nm);
temp_global_annual = NaN(nt, nm);
temp_global_shyear = NaN(nt, nm);
temp_global_monthly = NaN(nt, 18, nm);

% CCW
load ./data/gpp_ccw.mat;
temp_annual(:,:,:,1) = GPP_annual(:,:,years>=syear & years<=eyear) - repmat(nanmean(GPP_annual(:,:,years>=syear & years<=eyear), 3),1,1,nt);
temp_shyear(:,:,:,1) = GPP_shyear(:,:,years>=syear & years<=eyear) - repmat(nanmean(GPP_shyear(:,:,years>=syear & years<=eyear), 3),1,1,nt);
temp_global_annual(:, 1) = GPP_global_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_global_annual(years>=syear & years<=eyear)),nt,1);
temp_global_shyear(:, 1) = GPP_global_shyear(years>=syear & years<=eyear) - repmat(nanmean(GPP_global_shyear(years>=syear & years<=eyear)),nt,1);
temp_global_monthly(:,:,1) = GPP_global_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_global_monthly(years>=syear & years<=eyear, :)),nt,1);

% MOD17
load ./data/gpp_mod17.mat;
temp_annual(:,:,:,2) = GPP_annual(:,:,years>=syear & years<=eyear) - repmat(nanmean(GPP_annual(:,:,years>=syear & years<=eyear), 3),1,1,nt);
temp_shyear(:,:,:,2) = GPP_shyear(:,:,years>=syear & years<=eyear) - repmat(nanmean(GPP_shyear(:,:,years>=syear & years<=eyear), 3),1,1,nt);
temp_global_annual(:, 2) = GPP_global_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_global_annual(years>=syear & years<=eyear)),nt,1);
temp_global_shyear(:, 2) = GPP_global_shyear(years>=syear & years<=eyear) - repmat(nanmean(GPP_global_shyear(years>=syear & years<=eyear)),nt,1);
temp_global_monthly(:,:,2) = GPP_global_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_global_monthly(years>=syear & years<=eyear, :)),nt,1);

GPP_annual_mean = squeeze(nanmean(temp_annual, 4));
GPP_global_annual = temp_global_annual;
GPP_global_annual_mean = nanmean(GPP_global_annual,2);
GPP_shyear_mean = squeeze(nanmean(temp_shyear, 4));
GPP_global_shyear = temp_global_shyear;
GPP_global_shyear_mean = nanmean(GPP_global_shyear,2);
GPP_global_monthly = temp_global_monthly;
GPP_global_monthly_mean = nanmean(GPP_global_monthly,3);
years = syear:eyear;

clear temp* GPP_annual GPP_shyear GPP_monthly mo yr;

%% Regress monthly & annual gpp against CPI and EPI indices
load ./data/cpi_epi_1951-2016.mat;
[ny,nx,~] = size(GPP_annual_mean);

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

EP_GPP_global_monthly_beta = NaN(18, length(models));
EP_GPP_global_monthly_beta_CI = NaN(18, length(models));
EP_GPP_global_monthly_mean_beta = NaN(18, 1);
EP_GPP_global_monthly_mean_beta_CI = NaN(18, 1);
EP_GPP_global_annual_beta = NaN(1, length(models));
EP_GPP_global_annual_beta_CI = NaN(1, length(models));
EP_GPP_global_shyear_beta = NaN(1, length(models));
EP_GPP_global_shyear_beta_CI = NaN(1, length(models));

mdl = fitlm(epi, GPP_global_annual_mean);
EP_GPP_global_annual_mean_beta = mdl.Coefficients.Estimate(2);
mdl = fitlm(epi, GPP_global_shyear_mean);
EP_GPP_global_shyear_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:18
    mdl = fitlm(epi, GPP_global_monthly_mean(:, i));
    EP_GPP_global_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    EP_GPP_global_monthly_mean_beta_CI(i) = 1.96*mdl.Coefficients.SE(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_global_monthly(:, i, j));
        EP_GPP_global_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        EP_GPP_global_monthly_beta_CI(i, j) = 1.96*mdl.Coefficients.SE(2);
        if i == 1
            mdl = fitlm(epi, GPP_global_annual(:, j));
            EP_GPP_global_annual_beta(j) = mdl.Coefficients.Estimate(2);
            EP_GPP_global_annual_beta_CI(j) = 1.96*mdl.Coefficients.SE(2);
            disp(['Global Jan-Dec ',models{j},' response to EPI (PgC yr-1): ',...
                num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
                num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
            mdl = fitlm(epi, GPP_global_shyear(:, j));
            EP_GPP_global_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            EP_GPP_global_shyear_beta_CI(j) = 1.96*mdl.Coefficients.SE(2);
            disp(['Global Jul-Jun ',models{j},' response to EPI (PgC yr-1): ',...
                num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
                num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);

        end
    end
    
end

CP_GPP_global_monthly_beta = NaN(18, length(models));
CP_GPP_global_monthly_beta_CI = NaN(18, length(models));
CP_GPP_global_monthly_mean_beta = NaN(18, 1);
CP_GPP_global_monthly_mean_beta_CI = NaN(18, 1);
CP_GPP_global_annual_beta = NaN(1, length(models));
CP_GPP_global_annual_beta_CI = NaN(1, length(models));
CP_GPP_global_shyear_beta = NaN(1, length(models));
CP_GPP_global_shyear_beta_CI = NaN(1, length(models));

mdl = fitlm(cpi, GPP_global_annual_mean);
CP_GPP_global_annual_mean_beta = mdl.Coefficients.Estimate(2);
mdl = fitlm(cpi, GPP_global_shyear_mean);
CP_GPP_global_shyear_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:18
    mdl = fitlm(cpi, GPP_global_monthly_mean(:, i));
    CP_GPP_global_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    CP_GPP_global_monthly_mean_beta_CI(i) = 1.96*mdl.Coefficients.SE(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_global_monthly(:, i, j));
        CP_GPP_global_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        CP_GPP_global_monthly_beta_CI(i, j) = 1.96*mdl.Coefficients.SE(2);
        if i == 1
            mdl = fitlm(cpi, GPP_global_annual(:, j));
            CP_GPP_global_annual_beta(j) = mdl.Coefficients.Estimate(2);
            CP_GPP_global_annual_beta_CI(j) = 1.96*mdl.Coefficients.SE(2);
            disp(['Global Jan-Dec ',models{j},' response to CPI (PgC yr-1): ',...
                num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
                num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
            mdl = fitlm(cpi, GPP_global_shyear(:, j));
            CP_GPP_global_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            CP_GPP_global_shyear_beta_CI(j) = 1.96*mdl.Coefficients.SE(2);
            disp(['Global Jul-Jun ',models{j},' response to CPI (PgC yr-1): ',...
                num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
                num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
        end
    end
    
end

clear i j mdl;

%% Regress EPI and CPI vs. gridded GPP
EP_GPP_annual_r = NaN(ny, nx);
EP_GPP_annual_p = NaN(ny, nx);
EP_GPP_annual_beta = NaN(ny, nx);
CP_GPP_annual_r = NaN(ny, nx);
CP_GPP_annual_p = NaN(ny, nx);
CP_GPP_annual_beta = NaN(ny, nx);
EP_GPP_shyear_r = NaN(ny, nx);
EP_GPP_shyear_p = NaN(ny, nx);
EP_GPP_shyear_beta = NaN(ny, nx);
CP_GPP_shyear_r = NaN(ny, nx);
CP_GPP_shyear_p = NaN(ny, nx);
CP_GPP_shyear_beta = NaN(ny, nx);

for i = 1:ny
    for j = 1:nx
        ts = squeeze(GPP_annual_mean(i,j,:));
        if sum(isnan(ts))==0
            [r,p] = corr(epi, ts, 'rows','pairwise');
            mdl = fitlm(epi, ts);
            EP_GPP_annual_r(i,j) = r;
            EP_GPP_annual_p(i,j) = p;
            EP_GPP_annual_beta(i,j) = mdl.Coefficients.Estimate(2);
            
            [r,p] = corr(cpi, ts, 'rows','pairwise');
            mdl = fitlm(cpi, ts);
            CP_GPP_annual_r(i,j) = r;
            CP_GPP_annual_p(i,j) = p;
            CP_GPP_annual_beta(i,j) = mdl.Coefficients.Estimate(2);
        end
        
        ts = squeeze(GPP_shyear_mean(i,j,:));
        if sum(isnan(ts))<=1
            [r,p] = corr(epi, ts, 'rows','pairwise');
            mdl = fitlm(epi, ts);
            EP_GPP_shyear_r(i,j) = r;
            EP_GPP_shyear_p(i,j) = p;
            EP_GPP_shyear_beta(i,j) = mdl.Coefficients.Estimate(2);
            
            [r,p] = corr(cpi, ts, 'rows','pairwise');
            mdl = fitlm(cpi, ts);
            CP_GPP_shyear_r(i,j) = r;
            CP_GPP_shyear_p(i,j) = p;
            CP_GPP_shyear_beta(i,j) = mdl.Coefficients.Estimate(2);
        end
    end
end
clear i j mdl ts r p;

%% Now do global tropics
% temporary arrays for all three models
temp_tropics_annual = NaN(nt, nm);
temp_tropics_shyear = NaN(nt, nm);
temp_tropics_monthly = NaN(nt, 18, nm);

% CCW
load ./data/gpp_ccw_regional.mat;
clear temp *africa* *amazon* *austr* *casia* *eastus* *europe* *extratropical* *grass* *sahel* *semiarid* *tropical* *tundra* *westna*;
temp_tropics_annual(:, 1) = GPP_tropics_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_tropics_annual(years>=syear & years<=eyear)),nt,1);
temp_tropics_shyear(:, 1) = GPP_tropics_shyear(years>=syear & years<=eyear) - repmat(nanmean(GPP_tropics_shyear(years>=syear & years<=eyear)),nt,1);
temp_tropics_monthly(:,:,1) = GPP_tropics_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_tropics_monthly(years>=syear & years<=eyear, :)),nt,1);

% MOD17
load ./data/gpp_mod17_regional.mat;
clear temp *africa* *amazon* *austr* *casia* *eastus* *europe* *extratropical* *grass* *sahel* *semiarid* *tropical* *tundra* *westna*;
temp_tropics_annual(:, 2) = GPP_tropics_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_tropics_annual(years>=syear & years<=eyear)),nt,1);
temp_tropics_shyear(:, 2) = GPP_tropics_shyear(years>=syear & years<=eyear) - repmat(nanmean(GPP_tropics_shyear(years>=syear & years<=eyear)),nt,1);
temp_tropics_monthly(:,:,2) = GPP_tropics_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_tropics_monthly(years>=syear & years<=eyear, :)),nt,1);

GPP_tropics_annual = temp_tropics_annual;
GPP_tropics_annual_mean = nanmean(GPP_tropics_annual,2);
GPP_tropics_shyear = temp_tropics_shyear;
GPP_tropics_shyear_mean = nanmean(GPP_tropics_shyear,2);
GPP_tropics_monthly = temp_tropics_monthly;
GPP_tropics_monthly_mean = nanmean(GPP_tropics_monthly,3);
years = syear:eyear;

clear temp* GPP_annual GPP_shyear GPP_monthly mo yr;

% Regress monthly & annual gpp against CPI and EPI indices
load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

EP_GPP_tropics_monthly_beta = NaN(18, length(models));
EP_GPP_tropics_monthly_beta_CI = NaN(18, length(models));
EP_GPP_tropics_monthly_mean_beta = NaN(18, 1);
EP_GPP_tropics_monthly_mean_beta_CI = NaN(18, 1);
EP_GPP_tropics_annual_beta = NaN(1, length(models));
EP_GPP_tropics_annual_beta_CI = NaN(1, length(models));
EP_GPP_tropics_shyear_beta = NaN(1, length(models));
EP_GPP_tropics_shyear_beta_CI = NaN(1, length(models));

mdl = fitlm(epi, GPP_tropics_annual_mean);
EP_GPP_tropics_annual_mean_beta = mdl.Coefficients.Estimate(2);
mdl = fitlm(epi, GPP_tropics_shyear_mean);
EP_GPP_tropics_shyear_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:18
    mdl = fitlm(epi, GPP_tropics_monthly_mean(:, i));
    EP_GPP_tropics_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    EP_GPP_tropics_monthly_mean_beta_CI(i) = 1.96*mdl.Coefficients.SE(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_tropics_monthly(:, i, j));
        EP_GPP_tropics_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        EP_GPP_tropics_monthly_beta_CI(i, j) = 1.96*mdl.Coefficients.SE(2);
        if i == 1
            mdl = fitlm(epi, GPP_tropics_annual(:, j));
            EP_GPP_tropics_annual_beta(j) = mdl.Coefficients.Estimate(2);
            EP_GPP_tropics_annual_beta_CI(j) = 1.96*mdl.Coefficients.SE(2);
            disp(['Global tropics Jan-Dec ',models{j},' response to EPI (PgC yr-1): ',...
                num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
                num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
            mdl = fitlm(epi, GPP_tropics_shyear(:, j));
            EP_GPP_tropics_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            EP_GPP_tropics_shyear_beta_CI(j) = 1.96*mdl.Coefficients.SE(2);
            disp(['Global tropics Jul-Jun ',models{j},' response to EPI (PgC yr-1): ',...
                num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
                num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
        end
    end
    
end

CP_GPP_tropics_monthly_beta = NaN(18, length(models));
CP_GPP_tropics_monthly_beta_CI = NaN(18, length(models));
CP_GPP_tropics_monthly_mean_beta = NaN(18, 1);
CP_GPP_tropics_monthly_mean_beta_CI = NaN(18, 1);
CP_GPP_tropics_annual_beta = NaN(1, length(models));
CP_GPP_tropics_annual_beta_CI = NaN(1, length(models));
CP_GPP_tropics_shyear_beta = NaN(1, length(models));
CP_GPP_tropics_shyear_beta_CI = NaN(1, length(models));

mdl = fitlm(cpi, GPP_tropics_annual_mean);
CP_GPP_tropics_annual_mean_beta = mdl.Coefficients.Estimate(2);
mdl = fitlm(cpi, GPP_tropics_shyear_mean);
CP_GPP_tropics_shyear_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:18
    mdl = fitlm(cpi, GPP_tropics_monthly_mean(:, i));
    CP_GPP_tropics_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    CP_GPP_tropics_monthly_mean_beta_CI(i) = 1.96*mdl.Coefficients.SE(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_tropics_monthly(:, i, j));
        CP_GPP_tropics_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        CP_GPP_tropics_monthly_beta_CI(i, j) = 1.96*mdl.Coefficients.SE(2);
        if i == 1
            mdl = fitlm(cpi, GPP_tropics_annual(:, j));
            CP_GPP_tropics_annual_beta(j) = mdl.Coefficients.Estimate(2);
            CP_GPP_tropics_annual_beta_CI(j) = 1.96*mdl.Coefficients.SE(2);
            disp(['Global tropics Jan-Dec ',models{j},' response to CPI (PgC yr-1): ',...
                num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
                num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
            mdl = fitlm(cpi, GPP_tropics_shyear(:, j));
            CP_GPP_tropics_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            CP_GPP_tropics_shyear_beta_CI(j) = 1.96*mdl.Coefficients.SE(2);
            disp(['Global tropics Jul-Jun ',models{j},' response to CPI (PgC yr-1): ',...
                num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
                num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
        end
    end
    
end

clear i j ts r p mdl nt nx nm ny scale syear GPP mo;

save('./data/cp_ep_gpp_lue.mat');

