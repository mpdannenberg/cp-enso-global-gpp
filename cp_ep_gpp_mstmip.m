% Examine relationships between modeled MsTMIP GPP and EPI/CPI over 
% 1951-2010 period

%% Regress monthly & annual gpp against CPI and EPI indices
load ./data/gpp_mstmip.mat;
load ./data/cpi_epi_1951-2016.mat;

[ny,nx,~] = size(GPP_annual_mean);

cpi = cpi(yr<=2010);
epi = epi(yr<=2010);
yr = yr(yr<=2010);

EP_GPP_global_monthly_beta = NaN(12, length(models));
EP_GPP_global_monthly_mean_beta = NaN(12, 1);
EP_GPP_global_monthly_mean_beta_CI = NaN(12, 1);
EP_GPP_global_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_global_annual_mean);
EP_GPP_global_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_global_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Global MsTMIP mean response to EPI (PgC yr-1): ',...
    num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
    num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
for i = 1:12
    mdl = fitlm(epi, GPP_global_monthly_mean(:, i));
    EP_GPP_global_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    EP_GPP_global_monthly_mean_beta_CI(i) = 1.96*mdl.Coefficients.SE(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_global_monthly(:, i, j));
        EP_GPP_global_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_global_annual(:, j));
            EP_GPP_global_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_GPP_global_monthly_beta = NaN(12, length(models));
CP_GPP_global_monthly_mean_beta = NaN(12, 1);
CP_GPP_global_monthly_mean_beta_CI = NaN(12, 1);
CP_GPP_global_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_global_annual_mean);
CP_GPP_global_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_global_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Global MsTMIP mean response to CPI (PgC yr-1): ',...
    num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
    num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
for i = 1:12
    mdl = fitlm(cpi, GPP_global_monthly_mean(:, i));
    CP_GPP_global_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    CP_GPP_global_monthly_mean_beta_CI(i) = 1.96*mdl.Coefficients.SE(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_global_monthly(:, i, j));
        CP_GPP_global_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_global_annual(:, j));
            CP_GPP_global_annual_beta(j) = mdl.Coefficients.Estimate(2);
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

for i = 1:ny
    for j = 1:nx
        ts = squeeze(GPP_annual_mean(i,j,:));
        if sum(isnan(ts))==0
            [r,p] = corr(epi, ts);
            mdl = fitlm(epi, ts);
            EP_GPP_annual_r(i,j) = r;
            EP_GPP_annual_p(i,j) = p;
            EP_GPP_annual_beta(i,j) = mdl.Coefficients.Estimate(2);
            
            [r,p] = corr(cpi, ts);
            mdl = fitlm(cpi, ts);
            CP_GPP_annual_r(i,j) = r;
            CP_GPP_annual_p(i,j) = p;
            CP_GPP_annual_beta(i,j) = mdl.Coefficients.Estimate(2);
        end
    end
end

clear i j ts r p mdl nt nx ny scale syear GPP mo;

save('./data/cp_ep_gpp_mstmip.mat');

