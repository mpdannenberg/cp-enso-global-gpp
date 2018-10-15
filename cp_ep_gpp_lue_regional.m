% Examine relationships between data-driven GPP and EPI/CPI over 
% 1982-2015 period

syear = 1982;
eyear = 2016;

models = {'CCW','MOD17'}; nm = length(models);
nt = length(syear:eyear);


%% Regress monthly & annual gpp against CPI and EPI indices

%% Africa
% temporary arrays for all three models
temp_africa_annual = NaN(nt, nm);
temp_africa_monthly = NaN(nt, 12, nm);

% CCW
load ./data/gpp_ccw_regional.mat;
temp_africa_annual(:, 1) = GPP_africa_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_africa_annual(years>=syear & years<=eyear)),nt,1);
temp_africa_monthly(:,:,1) = GPP_africa_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_africa_monthly(years>=syear & years<=eyear, :)),nt,1);

% MOD17
load ./data/gpp_mod17_regional.mat;
temp_africa_annual(:, 2) = GPP_africa_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_africa_annual(years>=syear & years<=eyear)),nt,1);
temp_africa_monthly(:,:,2) = GPP_africa_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_africa_monthly(years>=syear & years<=eyear, :)),nt,1);

GPP_africa_annual = temp_africa_annual;
GPP_africa_annual_mean = nanmean(GPP_africa_annual,2);
GPP_africa_monthly = temp_africa_monthly;
GPP_africa_monthly_mean = nanmean(GPP_africa_monthly,3);
years = syear:eyear;

clear temp* mo yr;

load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

EP_GPP_africa_monthly_beta = NaN(12, length(models));
EP_GPP_africa_monthly_mean_beta = NaN(12, 1);
EP_GPP_africa_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_africa_annual_mean);
EP_GPP_africa_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(epi, GPP_africa_monthly_mean(:, i));
    EP_GPP_africa_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_africa_monthly(:, i, j));
        EP_GPP_africa_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_africa_annual(:, j));
            EP_GPP_africa_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_GPP_africa_monthly_beta = NaN(12, length(models));
CP_GPP_africa_monthly_mean_beta = NaN(12, 1);
CP_GPP_africa_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_africa_annual_mean);
CP_GPP_africa_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(cpi, GPP_africa_monthly_mean(:, i));
    CP_GPP_africa_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_africa_monthly(:, i, j));
        CP_GPP_africa_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_africa_annual(:, j));
            CP_GPP_africa_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

%% Africa
% temporary arrays for all three models
temp_amazon_annual = NaN(nt, nm);
temp_amazon_monthly = NaN(nt, 12, nm);

% CCW
load ./data/gpp_ccw_regional.mat;
temp_amazon_annual(:, 1) = GPP_amazon_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_amazon_annual(years>=syear & years<=eyear)),nt,1);
temp_amazon_monthly(:,:,1) = GPP_amazon_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_amazon_monthly(years>=syear & years<=eyear, :)),nt,1);

% MOD17
load ./data/gpp_mod17_regional.mat;
temp_amazon_annual(:, 2) = GPP_amazon_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_amazon_annual(years>=syear & years<=eyear)),nt,1);
temp_amazon_monthly(:,:,2) = GPP_amazon_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_amazon_monthly(years>=syear & years<=eyear, :)),nt,1);

GPP_amazon_annual = temp_amazon_annual;
GPP_amazon_annual_mean = nanmean(GPP_amazon_annual,2);
GPP_amazon_monthly = temp_amazon_monthly;
GPP_amazon_monthly_mean = nanmean(GPP_amazon_monthly,3);
years = syear:eyear;

clear temp* mo yr;

load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

EP_GPP_amazon_monthly_beta = NaN(12, length(models));
EP_GPP_amazon_monthly_mean_beta = NaN(12, 1);
EP_GPP_amazon_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_amazon_annual_mean);
EP_GPP_amazon_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(epi, GPP_amazon_monthly_mean(:, i));
    EP_GPP_amazon_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_amazon_monthly(:, i, j));
        EP_GPP_amazon_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_amazon_annual(:, j));
            EP_GPP_amazon_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_GPP_amazon_monthly_beta = NaN(12, length(models));
CP_GPP_amazon_monthly_mean_beta = NaN(12, 1);
CP_GPP_amazon_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_amazon_annual_mean);
CP_GPP_amazon_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(cpi, GPP_amazon_monthly_mean(:, i));
    CP_GPP_amazon_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_amazon_monthly(:, i, j));
        CP_GPP_amazon_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_amazon_annual(:, j));
            CP_GPP_amazon_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

%% Africa
% temporary arrays for all three models
temp_austr_annual = NaN(nt, nm);
temp_austr_monthly = NaN(nt, 12, nm);

% CCW
load ./data/gpp_ccw_regional.mat;
temp_austr_annual(:, 1) = GPP_austr_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_austr_annual(years>=syear & years<=eyear)),nt,1);
temp_austr_monthly(:,:,1) = GPP_austr_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_austr_monthly(years>=syear & years<=eyear, :)),nt,1);

% MOD17
load ./data/gpp_mod17_regional.mat;
temp_austr_annual(:, 2) = GPP_austr_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_austr_annual(years>=syear & years<=eyear)),nt,1);
temp_austr_monthly(:,:,2) = GPP_austr_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_austr_monthly(years>=syear & years<=eyear, :)),nt,1);

GPP_austr_annual = temp_austr_annual;
GPP_austr_annual_mean = nanmean(GPP_austr_annual,2);
GPP_austr_monthly = temp_austr_monthly;
GPP_austr_monthly_mean = nanmean(GPP_austr_monthly,3);
years = syear:eyear;

clear temp* mo yr;

load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

EP_GPP_austr_monthly_beta = NaN(12, length(models));
EP_GPP_austr_monthly_mean_beta = NaN(12, 1);
EP_GPP_austr_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_austr_annual_mean);
EP_GPP_austr_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(epi, GPP_austr_monthly_mean(:, i));
    EP_GPP_austr_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_austr_monthly(:, i, j));
        EP_GPP_austr_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_austr_annual(:, j));
            EP_GPP_austr_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_GPP_austr_monthly_beta = NaN(12, length(models));
CP_GPP_austr_monthly_mean_beta = NaN(12, 1);
CP_GPP_austr_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_austr_annual_mean);
CP_GPP_austr_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(cpi, GPP_austr_monthly_mean(:, i));
    CP_GPP_austr_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_austr_monthly(:, i, j));
        CP_GPP_austr_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_austr_annual(:, j));
            CP_GPP_austr_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

%% Africa
% temporary arrays for all three models
temp_casia_annual = NaN(nt, nm);
temp_casia_monthly = NaN(nt, 12, nm);

% CCW
load ./data/gpp_ccw_regional.mat;
temp_casia_annual(:, 1) = GPP_casia_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_casia_annual(years>=syear & years<=eyear)),nt,1);
temp_casia_monthly(:,:,1) = GPP_casia_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_casia_monthly(years>=syear & years<=eyear, :)),nt,1);

% MOD17
load ./data/gpp_mod17_regional.mat;
temp_casia_annual(:, 2) = GPP_casia_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_casia_annual(years>=syear & years<=eyear)),nt,1);
temp_casia_monthly(:,:,2) = GPP_casia_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_casia_monthly(years>=syear & years<=eyear, :)),nt,1);

GPP_casia_annual = temp_casia_annual;
GPP_casia_annual_mean = nanmean(GPP_casia_annual,2);
GPP_casia_monthly = temp_casia_monthly;
GPP_casia_monthly_mean = nanmean(GPP_casia_monthly,3);
years = syear:eyear;

clear temp* mo yr;

load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

EP_GPP_casia_monthly_beta = NaN(12, length(models));
EP_GPP_casia_monthly_mean_beta = NaN(12, 1);
EP_GPP_casia_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_casia_annual_mean);
EP_GPP_casia_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(epi, GPP_casia_monthly_mean(:, i));
    EP_GPP_casia_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_casia_monthly(:, i, j));
        EP_GPP_casia_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_casia_annual(:, j));
            EP_GPP_casia_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_GPP_casia_monthly_beta = NaN(12, length(models));
CP_GPP_casia_monthly_mean_beta = NaN(12, 1);
CP_GPP_casia_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_casia_annual_mean);
CP_GPP_casia_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(cpi, GPP_casia_monthly_mean(:, i));
    CP_GPP_casia_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_casia_monthly(:, i, j));
        CP_GPP_casia_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_casia_annual(:, j));
            CP_GPP_casia_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

%% Africa
% temporary arrays for all three models
temp_eastus_annual = NaN(nt, nm);
temp_eastus_monthly = NaN(nt, 12, nm);

% CCW
load ./data/gpp_ccw_regional.mat;
temp_eastus_annual(:, 1) = GPP_eastus_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_eastus_annual(years>=syear & years<=eyear)),nt,1);
temp_eastus_monthly(:,:,1) = GPP_eastus_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_eastus_monthly(years>=syear & years<=eyear, :)),nt,1);

% MOD17
load ./data/gpp_mod17_regional.mat;
temp_eastus_annual(:, 2) = GPP_eastus_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_eastus_annual(years>=syear & years<=eyear)),nt,1);
temp_eastus_monthly(:,:,2) = GPP_eastus_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_eastus_monthly(years>=syear & years<=eyear, :)),nt,1);

GPP_eastus_annual = temp_eastus_annual;
GPP_eastus_annual_mean = nanmean(GPP_eastus_annual,2);
GPP_eastus_monthly = temp_eastus_monthly;
GPP_eastus_monthly_mean = nanmean(GPP_eastus_monthly,3);
years = syear:eyear;

clear temp* mo yr;

load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

EP_GPP_eastus_monthly_beta = NaN(12, length(models));
EP_GPP_eastus_monthly_mean_beta = NaN(12, 1);
EP_GPP_eastus_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_eastus_annual_mean);
EP_GPP_eastus_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(epi, GPP_eastus_monthly_mean(:, i));
    EP_GPP_eastus_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_eastus_monthly(:, i, j));
        EP_GPP_eastus_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_eastus_annual(:, j));
            EP_GPP_eastus_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_GPP_eastus_monthly_beta = NaN(12, length(models));
CP_GPP_eastus_monthly_mean_beta = NaN(12, 1);
CP_GPP_eastus_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_eastus_annual_mean);
CP_GPP_eastus_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(cpi, GPP_eastus_monthly_mean(:, i));
    CP_GPP_eastus_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_eastus_monthly(:, i, j));
        CP_GPP_eastus_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_eastus_annual(:, j));
            CP_GPP_eastus_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

%% Africa
% temporary arrays for all three models
temp_europe_annual = NaN(nt, nm);
temp_europe_monthly = NaN(nt, 12, nm);

% CCW
load ./data/gpp_ccw_regional.mat;
temp_europe_annual(:, 1) = GPP_europe_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_europe_annual(years>=syear & years<=eyear)),nt,1);
temp_europe_monthly(:,:,1) = GPP_europe_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_europe_monthly(years>=syear & years<=eyear, :)),nt,1);

% MOD17
load ./data/gpp_mod17_regional.mat;
temp_europe_annual(:, 2) = GPP_europe_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_europe_annual(years>=syear & years<=eyear)),nt,1);
temp_europe_monthly(:,:,2) = GPP_europe_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_europe_monthly(years>=syear & years<=eyear, :)),nt,1);

GPP_europe_annual = temp_europe_annual;
GPP_europe_annual_mean = nanmean(GPP_europe_annual,2);
GPP_europe_monthly = temp_europe_monthly;
GPP_europe_monthly_mean = nanmean(GPP_europe_monthly,3);
years = syear:eyear;

clear temp* mo yr;

load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

EP_GPP_europe_monthly_beta = NaN(12, length(models));
EP_GPP_europe_monthly_mean_beta = NaN(12, 1);
EP_GPP_europe_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_europe_annual_mean);
EP_GPP_europe_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(epi, GPP_europe_monthly_mean(:, i));
    EP_GPP_europe_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_europe_monthly(:, i, j));
        EP_GPP_europe_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_europe_annual(:, j));
            EP_GPP_europe_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_GPP_europe_monthly_beta = NaN(12, length(models));
CP_GPP_europe_monthly_mean_beta = NaN(12, 1);
CP_GPP_europe_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_europe_annual_mean);
CP_GPP_europe_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(cpi, GPP_europe_monthly_mean(:, i));
    CP_GPP_europe_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_europe_monthly(:, i, j));
        CP_GPP_europe_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_europe_annual(:, j));
            CP_GPP_europe_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

%% Africa
% temporary arrays for all three models
temp_sahel_annual = NaN(nt, nm);
temp_sahel_monthly = NaN(nt, 12, nm);

% CCW
load ./data/gpp_ccw_regional.mat;
temp_sahel_annual(:, 1) = GPP_sahel_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_sahel_annual(years>=syear & years<=eyear)),nt,1);
temp_sahel_monthly(:,:,1) = GPP_sahel_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_sahel_monthly(years>=syear & years<=eyear, :)),nt,1);

% MOD17
load ./data/gpp_mod17_regional.mat;
temp_sahel_annual(:, 2) = GPP_sahel_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_sahel_annual(years>=syear & years<=eyear)),nt,1);
temp_sahel_monthly(:,:,2) = GPP_sahel_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_sahel_monthly(years>=syear & years<=eyear, :)),nt,1);

GPP_sahel_annual = temp_sahel_annual;
GPP_sahel_annual_mean = nanmean(GPP_sahel_annual,2);
GPP_sahel_monthly = temp_sahel_monthly;
GPP_sahel_monthly_mean = nanmean(GPP_sahel_monthly,3);
years = syear:eyear;

clear temp* mo yr;

load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

EP_GPP_sahel_monthly_beta = NaN(12, length(models));
EP_GPP_sahel_monthly_mean_beta = NaN(12, 1);
EP_GPP_sahel_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_sahel_annual_mean);
EP_GPP_sahel_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(epi, GPP_sahel_monthly_mean(:, i));
    EP_GPP_sahel_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_sahel_monthly(:, i, j));
        EP_GPP_sahel_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_sahel_annual(:, j));
            EP_GPP_sahel_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_GPP_sahel_monthly_beta = NaN(12, length(models));
CP_GPP_sahel_monthly_mean_beta = NaN(12, 1);
CP_GPP_sahel_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_sahel_annual_mean);
CP_GPP_sahel_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(cpi, GPP_sahel_monthly_mean(:, i));
    CP_GPP_sahel_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_sahel_monthly(:, i, j));
        CP_GPP_sahel_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_sahel_annual(:, j));
            CP_GPP_sahel_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

%% Africa
% temporary arrays for all three models
temp_westna_annual = NaN(nt, nm);
temp_westna_monthly = NaN(nt, 12, nm);

% CCW
load ./data/gpp_ccw_regional.mat;
temp_westna_annual(:, 1) = GPP_westna_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_westna_annual(years>=syear & years<=eyear)),nt,1);
temp_westna_monthly(:,:,1) = GPP_westna_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_westna_monthly(years>=syear & years<=eyear, :)),nt,1);

% MOD17
load ./data/gpp_mod17_regional.mat;
temp_westna_annual(:, 2) = GPP_westna_annual(years>=syear & years<=eyear) - repmat(nanmean(GPP_westna_annual(years>=syear & years<=eyear)),nt,1);
temp_westna_monthly(:,:,2) = GPP_westna_monthly(years>=syear & years<=eyear, :) - repmat(nanmean(GPP_westna_monthly(years>=syear & years<=eyear, :)),nt,1);

GPP_westna_annual = temp_westna_annual;
GPP_westna_annual_mean = nanmean(GPP_westna_annual,2);
GPP_westna_monthly = temp_westna_monthly;
GPP_westna_monthly_mean = nanmean(GPP_westna_monthly,3);
years = syear:eyear;

clear temp* mo yr;

load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

EP_GPP_westna_monthly_beta = NaN(12, length(models));
EP_GPP_westna_monthly_mean_beta = NaN(12, 1);
EP_GPP_westna_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_westna_annual_mean);
EP_GPP_westna_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(epi, GPP_westna_monthly_mean(:, i));
    EP_GPP_westna_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_westna_monthly(:, i, j));
        EP_GPP_westna_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_westna_annual(:, j));
            EP_GPP_westna_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_GPP_westna_monthly_beta = NaN(12, length(models));
CP_GPP_westna_monthly_mean_beta = NaN(12, 1);
CP_GPP_westna_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_westna_annual_mean);
CP_GPP_westna_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(cpi, GPP_westna_monthly_mean(:, i));
    CP_GPP_westna_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_westna_monthly(:, i, j));
        CP_GPP_westna_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_westna_annual(:, j));
            CP_GPP_westna_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl GPP* syear eyear nt yr nm;


%% Save
save('./data/cp_ep_gpp_lue_regional.mat');

