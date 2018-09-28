% Examine relationships between modeled MsTMIP GPP and EPI/CPI over 
% 1951-2010 period
syear = 1979;
eyear = 2013;

models = {'v15r2','v15r4','v13r1','v14r2'}; nm = length(models);

%% Regress monthly & annual gpp against CPI and EPI indices
load ./data/nep_inversions.mat;
load ./data/cpi_epi_1951-2016.mat;

[ny,nx,~] = size(NEP_annual_mean);

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

EP_NEP_global_monthly_beta = NaN(12, length(models));
EP_NEP_global_monthly_mean_beta = NaN(12, 1);
EP_NEP_global_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_global_annual_mean);
EP_NEP_global_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(epi, NEP_global_monthly_mean(:, i));
    EP_NEP_global_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_global_monthly(:, i, j));
        EP_NEP_global_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_global_annual(:, j));
            EP_NEP_global_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_NEP_global_monthly_beta = NaN(12, length(models));
CP_NEP_global_monthly_mean_beta = NaN(12, 1);
CP_NEP_global_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_global_annual_mean);
CP_NEP_global_annual_mean_beta = mdl.Coefficients.Estimate(2);
for i = 1:12
    mdl = fitlm(cpi, NEP_global_monthly_mean(:, i));
    CP_NEP_global_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_global_monthly(:, i, j));
        CP_NEP_global_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_global_annual(:, j));
            CP_NEP_global_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

%% Regress EPI and CPI vs. gridded NEP
EP_NEP_annual_r = NaN(ny, nx);
EP_NEP_annual_p = NaN(ny, nx);
EP_NEP_annual_beta = NaN(ny, nx);
CP_NEP_annual_r = NaN(ny, nx);
CP_NEP_annual_p = NaN(ny, nx);
CP_NEP_annual_beta = NaN(ny, nx);

for i = 1:ny
    for j = 1:nx
        ts = squeeze(NEP_annual_mean(i,j,:));
        if sum(isnan(ts))==0
            [r,p] = corr(epi, ts);
            mdl = fitlm(epi, ts);
            EP_NEP_annual_r(i,j) = r;
            EP_NEP_annual_p(i,j) = p;
            EP_NEP_annual_beta(i,j) = mdl.Coefficients.Estimate(2);
            
            [r,p] = corr(cpi, ts);
            mdl = fitlm(cpi, ts);
            CP_NEP_annual_r(i,j) = r;
            CP_NEP_annual_p(i,j) = p;
            CP_NEP_annual_beta(i,j) = mdl.Coefficients.Estimate(2);
        end
    end
end

clear i j ts r p mdl nt nx ny scale syear NEP mo;

save('./data/cp_ep_nep_inversions.mat');
