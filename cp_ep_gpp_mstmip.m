% Examine relationships between modeled MsTMIP GPP and EPI/CPI over 
% 1951-2010 period

syear = 1951;
eyear = 2010;

%% Regress monthly & annual gpp against CPI and EPI indices
load ./data/gpp_mstmip.mat;
load ./data/cpi_epi_1951-2016.mat;

[ny,nx,~] = size(GPP_annual_mean);

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

GPP_annual_mean = GPP_annual_mean(:,:,years>=syear & years<=eyear);
GPP_global_annual = GPP_global_annual(years>=syear & years<=eyear, :);
GPP_global_annual_mean = GPP_global_annual_mean(years>=syear & years<=eyear, :);
GPP_shyear_mean = GPP_shyear_mean(:,:,years>=syear & years<=eyear);
GPP_global_shyear = GPP_global_shyear(years>=syear & years<=eyear, :);
GPP_global_shyear_mean = GPP_global_shyear_mean(years>=syear & years<=eyear, :);
GPP_global_monthly = GPP_global_monthly(years>=syear & years<=eyear, :, :);
GPP_global_monthly_mean = GPP_global_monthly_mean(years>=syear & years<=eyear, :);
years = years(years>=syear & years<=eyear)';

EP_GPP_global_monthly_beta = NaN(12, length(models));
EP_GPP_global_monthly_mean_beta = NaN(12, 1);
EP_GPP_global_monthly_mean_beta_CI = NaN(12, 1);
EP_GPP_global_annual_beta = NaN(1, length(models));
EP_GPP_global_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_global_annual_mean);
EP_GPP_global_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_global_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Global Jan-Dec MsTMIP mean response to EPI (PgC yr-1): ',...
    num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
    num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
mdl = fitlm(epi, GPP_global_shyear_mean);
EP_GPP_global_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_global_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Global Jul-Jun MsTMIP mean response to EPI (PgC yr-1): ',...
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
            mdl = fitlm(epi, GPP_global_shyear(:, j));
            EP_GPP_global_shyear_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_GPP_global_monthly_beta = NaN(12, length(models));
CP_GPP_global_monthly_mean_beta = NaN(12, 1);
CP_GPP_global_monthly_mean_beta_CI = NaN(12, 1);
CP_GPP_global_annual_beta = NaN(1, length(models));
CP_GPP_global_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_global_annual_mean);
CP_GPP_global_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_global_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Global Jan-Dec MsTMIP mean response to CPI (PgC yr-1): ',...
    num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
    num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
mdl = fitlm(cpi, GPP_global_shyear_mean);
CP_GPP_global_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_global_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Global Jul-Jun MsTMIP mean response to CPI (PgC yr-1): ',...
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
            mdl = fitlm(cpi, GPP_global_shyear(:, j));
            CP_GPP_global_shyear_beta(j) = mdl.Coefficients.Estimate(2);
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
        if sum(isnan(ts))==0
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


%% Regress monthly & annual gpp against CPI and EPI indices
load ./data/gpp_mstmip_regional.mat;
clear temp *africa* *amazon* *austr* *casia* *eastus* *europe* *extratropical* *grass* *sahel* *semiarid* *tropical* *tundra* *westna*;
load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

GPP_tropics_annual = GPP_tropics_annual(years>=syear & years<=eyear, :);
GPP_tropics_annual_mean = GPP_tropics_annual_mean(years>=syear & years<=eyear, :);
GPP_tropics_shyear = GPP_tropics_shyear(years>=syear & years<=eyear, :);
GPP_tropics_shyear_mean = GPP_tropics_shyear_mean(years>=syear & years<=eyear, :);
GPP_tropics_monthly = GPP_tropics_monthly(years>=syear & years<=eyear, :, :);
GPP_tropics_monthly_mean = GPP_tropics_monthly_mean(years>=syear & years<=eyear, :);
years = years(years>=syear & years<=eyear)';

EP_GPP_tropics_monthly_beta = NaN(12, length(models));
EP_GPP_tropics_monthly_mean_beta = NaN(12, 1);
EP_GPP_tropics_monthly_mean_beta_CI = NaN(12, 1);
EP_GPP_tropics_annual_beta = NaN(1, length(models));
EP_GPP_tropics_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_tropics_annual_mean);
EP_GPP_tropics_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_tropics_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Tropical Jan-Dec MsTMIP mean response to EPI (PgC yr-1): ',...
    num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
    num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
mdl = fitlm(epi, GPP_tropics_shyear_mean);
EP_GPP_tropics_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_tropics_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Tropical Jul-Jun MsTMIP mean response to EPI (PgC yr-1): ',...
    num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
    num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
for i = 1:12
    mdl = fitlm(epi, GPP_tropics_monthly_mean(:, i));
    EP_GPP_tropics_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    EP_GPP_tropics_monthly_mean_beta_CI(i) = 1.96*mdl.Coefficients.SE(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_tropics_monthly(:, i, j));
        EP_GPP_tropics_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_tropics_annual(:, j));
            EP_GPP_tropics_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_tropics_shyear(:, j));
            EP_GPP_tropics_shyear_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_GPP_tropics_monthly_beta = NaN(12, length(models));
CP_GPP_tropics_monthly_mean_beta = NaN(12, 1);
CP_GPP_tropics_monthly_mean_beta_CI = NaN(12, 1);
CP_GPP_tropics_annual_beta = NaN(1, length(models));
CP_GPP_tropics_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_tropics_annual_mean);
CP_GPP_tropics_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_tropics_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Tropical Jan-Dec MsTMIP mean response to CPI (PgC yr-1): ',...
    num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
    num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
mdl = fitlm(cpi, GPP_tropics_shyear_mean);
CP_GPP_tropics_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_tropics_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Tropical Jul-Jun MsTMIP mean response to CPI (PgC yr-1): ',...
    num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
    num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
for i = 1:12
    mdl = fitlm(cpi, GPP_tropics_monthly_mean(:, i));
    CP_GPP_tropics_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    CP_GPP_tropics_monthly_mean_beta_CI(i) = 1.96*mdl.Coefficients.SE(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_tropics_monthly(:, i, j));
        CP_GPP_tropics_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_tropics_annual(:, j));
            CP_GPP_tropics_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_tropics_shyear(:, j));
            CP_GPP_tropics_shyear_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear eyear i j ts r p mdl nt nx ny scale syear GPP mo;

save('./data/cp_ep_gpp_mstmip.mat');

