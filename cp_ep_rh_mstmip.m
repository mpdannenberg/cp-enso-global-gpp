% Examine relationships between modeled MsTMIP Rh and EPI/CPI over 
% 1951-2010 period

syear = 1951;
eyear = 2010;

%% Regress monthly & annual gpp against CPI and EPI indices
load ./data/rh_mstmip.mat;
load ./data/cpi_epi_1951-2016.mat;

[ny,nx,~] = size(Rh_annual_mean);

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

Rh_annual_mean = Rh_annual_mean(:,:,years>=syear & years<=eyear);
Rh_global_annual = Rh_global_annual(years>=syear & years<=eyear, :);
Rh_global_annual_mean = Rh_global_annual_mean(years>=syear & years<=eyear, :);
Rh_global_monthly = Rh_global_monthly(years>=syear & years<=eyear, :, :);
Rh_global_monthly_mean = Rh_global_monthly_mean(years>=syear & years<=eyear, :);
years = years(years>=syear & years<=eyear)';

EP_Rh_global_monthly_beta = NaN(12, length(models));
EP_Rh_global_monthly_mean_beta = NaN(12, 1);
EP_Rh_global_monthly_mean_beta_CI = NaN(12, 1);
EP_Rh_global_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_global_annual_mean);
EP_Rh_global_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_global_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Global MsTMIP mean response to EPI (PgC yr-1): ',...
    num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
    num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
for i = 1:12
    mdl = fitlm(epi, Rh_global_monthly_mean(:, i));
    EP_Rh_global_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    EP_Rh_global_monthly_mean_beta_CI(i) = 1.96*mdl.Coefficients.SE(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_global_monthly(:, i, j));
        EP_Rh_global_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_global_annual(:, j));
            EP_Rh_global_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Rh_global_monthly_beta = NaN(12, length(models));
CP_Rh_global_monthly_mean_beta = NaN(12, 1);
CP_Rh_global_monthly_mean_beta_CI = NaN(12, 1);
CP_Rh_global_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_global_annual_mean);
CP_Rh_global_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_global_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Global MsTMIP mean response to CPI (PgC yr-1): ',...
    num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
    num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
for i = 1:12
    mdl = fitlm(cpi, Rh_global_monthly_mean(:, i));
    CP_Rh_global_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    CP_Rh_global_monthly_mean_beta_CI(i) = 1.96*mdl.Coefficients.SE(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_global_monthly(:, i, j));
        CP_Rh_global_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_global_annual(:, j));
            CP_Rh_global_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

%% Regress EPI and CPI vs. gridded Rh
EP_Rh_annual_r = NaN(ny, nx);
EP_Rh_annual_p = NaN(ny, nx);
EP_Rh_annual_beta = NaN(ny, nx);
CP_Rh_annual_r = NaN(ny, nx);
CP_Rh_annual_p = NaN(ny, nx);
CP_Rh_annual_beta = NaN(ny, nx);

for i = 1:ny
    for j = 1:nx
        ts = squeeze(Rh_annual_mean(i,j,:));
        if sum(isnan(ts))==0
            [r,p] = corr(epi, ts);
            mdl = fitlm(epi, ts);
            EP_Rh_annual_r(i,j) = r;
            EP_Rh_annual_p(i,j) = p;
            EP_Rh_annual_beta(i,j) = mdl.Coefficients.Estimate(2);
            
            [r,p] = corr(cpi, ts);
            mdl = fitlm(cpi, ts);
            CP_Rh_annual_r(i,j) = r;
            CP_Rh_annual_p(i,j) = p;
            CP_Rh_annual_beta(i,j) = mdl.Coefficients.Estimate(2);
        end
    end
end
clear i j mdl ts r p;


%% Regress monthly & annual rh against CPI and EPI indices
load ./data/rh_mstmip_regional.mat;
clear temp *africa* *amazon* *austr* *casia* *eastus* *europe* *extratropical* *grass* *sahel* *semiarid* *tropical* *tundra* *westna*;
load ./data/cpi_epi_1951-2016.mat;

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

Rh_tropics_annual = Rh_tropics_annual(years>=syear & years<=eyear, :);
Rh_tropics_annual_mean = Rh_tropics_annual_mean(years>=syear & years<=eyear, :);
Rh_tropics_monthly = Rh_tropics_monthly(years>=syear & years<=eyear, :, :);
Rh_tropics_monthly_mean = Rh_tropics_monthly_mean(years>=syear & years<=eyear, :);
years = years(years>=syear & years<=eyear)';

EP_Rh_tropics_monthly_beta = NaN(12, length(models));
EP_Rh_tropics_monthly_mean_beta = NaN(12, 1);
EP_Rh_tropics_monthly_mean_beta_CI = NaN(12, 1);
EP_Rh_tropics_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_tropics_annual_mean);
EP_Rh_tropics_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_tropics_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Tropical MsTMIP mean response to EPI (PgC yr-1): ',...
    num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
    num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
for i = 1:12
    mdl = fitlm(epi, Rh_tropics_monthly_mean(:, i));
    EP_Rh_tropics_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    EP_Rh_tropics_monthly_mean_beta_CI(i) = 1.96*mdl.Coefficients.SE(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_tropics_monthly(:, i, j));
        EP_Rh_tropics_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_tropics_annual(:, j));
            EP_Rh_tropics_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Rh_tropics_monthly_beta = NaN(12, length(models));
CP_Rh_tropics_monthly_mean_beta = NaN(12, 1);
CP_Rh_tropics_monthly_mean_beta_CI = NaN(12, 1);
CP_Rh_tropics_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_tropics_annual_mean);
CP_Rh_tropics_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_tropics_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
disp(['Tropical MsTMIP mean response to CPI (PgC yr-1): ',...
    num2str(round(mdl.Coefficients.Estimate(2)/1000,2)),' +/- ',...
    num2str(round(1.96*mdl.Coefficients.SE(2)/1000,2))]);
for i = 1:12
    mdl = fitlm(cpi, Rh_tropics_monthly_mean(:, i));
    CP_Rh_tropics_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    CP_Rh_tropics_monthly_mean_beta_CI(i) = 1.96*mdl.Coefficients.SE(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_tropics_monthly(:, i, j));
        CP_Rh_tropics_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_tropics_annual(:, j));
            CP_Rh_tropics_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear eyear i j ts r p mdl nt nx ny scale syear Rh mo;

save('./data/cp_ep_rh_mstmip.mat');

