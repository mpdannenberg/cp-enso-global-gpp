% Examine relationships between modeled MsTMIP GPP and EPI/CPI over 
% 1951-2010 period

%% Regress monthly & annual gpp against CPI and EPI indices
load ./data/gpp_mstmip_regional.mat;
load ./data/cpi_epi_1951-2016.mat;

[ny,nx,~] = size(GPP_annual_mean);

cpi = cpi(yr<=2010);
epi = epi(yr<=2010);
yr = yr(yr<=2010);

% Africa
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

% Africa
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

% Africa
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

% Africa
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

% Africa
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

% Africa
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

% Africa
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

% Africa
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

clear i j mdl GPP*;


save('./data/cp_ep_gpp_mstmip_regional.mat');

