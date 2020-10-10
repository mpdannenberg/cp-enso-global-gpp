% Examine relationships between modeled MsTMIP GPP and EPI/CPI over 
% 1951-2010 period

syear = 1951;
eyear = 2010;

%% Regress monthly & annual gpp against CPI and EPI indices
load ./data/gpp_mstmip_regional.mat;
load ./data/cpi_epi_1951-2016.mat;

ep_gpp_95CI_annual = cell2table(cell(13,1), 'VariableNames',{'MsTMIP'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US',...
    'Europe','The Sahel','Western North America','Tropical Forests',...
    'Extratropical Forests','Tundra/Arctic Shrubland','Grass/Crops',...
    'Semiarid'});
cp_gpp_95CI_annual = cell2table(cell(13,1), 'VariableNames',{'MsTMIP'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US',...
    'Europe','The Sahel','Western North America','Tropical Forests',...
    'Extratropical Forests','Tundra/Arctic Shrubland','Grass/Crops',...
    'Semiarid'});
ep_gpp_95CI_shyear = cell2table(cell(13,1), 'VariableNames',{'MsTMIP'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US',...
    'Europe','The Sahel','Western North America','Tropical Forests',...
    'Extratropical Forests','Tundra/Arctic Shrubland','Grass/Crops',...
    'Semiarid'});
cp_gpp_95CI_shyear = cell2table(cell(13,1), 'VariableNames',{'MsTMIP'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US',...
    'Europe','The Sahel','Western North America','Tropical Forests',...
    'Extratropical Forests','Tundra/Arctic Shrubland','Grass/Crops',...
    'Semiarid'});

[ny,nx,~] = size(GPP_annual_mean);

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

idx = years>=syear & years<=eyear;

% Africa
EP_GPP_africa_monthly_beta = NaN(12, length(models));
EP_GPP_africa_monthly_mean_beta = NaN(12, 1);
EP_GPP_africa_annual_beta = NaN(1, length(models));
EP_GPP_africa_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_africa_annual_mean(idx));
EP_GPP_africa_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_africa_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_annual{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, GPP_africa_shyear_mean(idx));
EP_GPP_africa_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_africa_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_shyear{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, GPP_africa_monthly_mean(idx, i));
    EP_GPP_africa_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_africa_monthly(idx, i, j));
        EP_GPP_africa_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_africa_annual(idx, j));
            EP_GPP_africa_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_africa_shyear(idx, j));
            EP_GPP_africa_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_GPP_africa_monthly_beta = NaN(12, length(models));
CP_GPP_africa_monthly_mean_beta = NaN(12, 1);
CP_GPP_africa_annual_beta = NaN(1, length(models));
CP_GPP_africa_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_africa_annual_mean(idx));
CP_GPP_africa_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_africa_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_annual{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, GPP_africa_shyear_mean(idx));
CP_GPP_africa_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_africa_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_shyear{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, GPP_africa_monthly_mean(idx, i));
    CP_GPP_africa_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_africa_monthly(idx, i, j));
        CP_GPP_africa_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_africa_annual(idx, j));
            CP_GPP_africa_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_africa_shyear(idx, j));
            CP_GPP_africa_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Amazonia
EP_GPP_amazon_monthly_beta = NaN(12, length(models));
EP_GPP_amazon_monthly_mean_beta = NaN(12, 1);
EP_GPP_amazon_annual_beta = NaN(1, length(models));
EP_GPP_amazon_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_amazon_annual_mean(idx));
EP_GPP_amazon_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_amazon_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_annual{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, GPP_amazon_shyear_mean(idx));
EP_GPP_amazon_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_amazon_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_shyear{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, GPP_amazon_monthly_mean(idx, i));
    EP_GPP_amazon_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_amazon_monthly(idx, i, j));
        EP_GPP_amazon_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_amazon_annual(idx, j));
            EP_GPP_amazon_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_amazon_shyear(idx, j));
            EP_GPP_amazon_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_GPP_amazon_monthly_beta = NaN(12, length(models));
CP_GPP_amazon_monthly_mean_beta = NaN(12, 1);
CP_GPP_amazon_annual_beta = NaN(1, length(models));
CP_GPP_amazon_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_amazon_annual_mean(idx));
CP_GPP_amazon_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_amazon_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_annual{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, GPP_amazon_shyear_mean(idx));
CP_GPP_amazon_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_amazon_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_shyear{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, GPP_amazon_monthly_mean(idx, i));
    CP_GPP_amazon_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_amazon_monthly(idx, i, j));
        CP_GPP_amazon_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_amazon_annual(idx, j));
            CP_GPP_amazon_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_amazon_shyear(idx, j));
            CP_GPP_amazon_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Australia
EP_GPP_austr_monthly_beta = NaN(12, length(models));
EP_GPP_austr_monthly_mean_beta = NaN(12, 1);
EP_GPP_austr_annual_beta = NaN(1, length(models));
EP_GPP_austr_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_austr_annual_mean(idx));
EP_GPP_austr_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_austr_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_annual{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, GPP_austr_shyear_mean(idx));
EP_GPP_austr_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_austr_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_shyear{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, GPP_austr_monthly_mean(idx, i));
    EP_GPP_austr_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_austr_monthly(idx, i, j));
        EP_GPP_austr_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_austr_annual(idx, j));
            EP_GPP_austr_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_austr_shyear(idx, j));
            EP_GPP_austr_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_GPP_austr_monthly_beta = NaN(12, length(models));
CP_GPP_austr_monthly_mean_beta = NaN(12, 1);
CP_GPP_austr_annual_beta = NaN(1, length(models));
CP_GPP_austr_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_austr_annual_mean(idx));
CP_GPP_austr_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_austr_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_annual{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, GPP_austr_shyear_mean(idx));
CP_GPP_austr_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_austr_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_shyear{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, GPP_austr_monthly_mean(idx, i));
    CP_GPP_austr_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_austr_monthly(idx, i, j));
        CP_GPP_austr_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_austr_annual(idx, j));
            CP_GPP_austr_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_austr_shyear(idx, j));
            CP_GPP_austr_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Central Asia
EP_GPP_casia_monthly_beta = NaN(12, length(models));
EP_GPP_casia_monthly_mean_beta = NaN(12, 1);
EP_GPP_casia_annual_beta = NaN(1, length(models));
EP_GPP_casia_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_casia_annual_mean(idx));
EP_GPP_casia_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_casia_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_annual{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, GPP_casia_shyear_mean(idx));
EP_GPP_casia_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_casia_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_shyear{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, GPP_casia_monthly_mean(idx, i));
    EP_GPP_casia_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_casia_monthly(idx, i, j));
        EP_GPP_casia_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_casia_annual(idx, j));
            EP_GPP_casia_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_casia_shyear(idx, j));
            EP_GPP_casia_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_GPP_casia_monthly_beta = NaN(12, length(models));
CP_GPP_casia_monthly_mean_beta = NaN(12, 1);
CP_GPP_casia_annual_beta = NaN(1, length(models));
CP_GPP_casia_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_casia_annual_mean(idx));
CP_GPP_casia_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_casia_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_annual{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, GPP_casia_shyear_mean(idx));
CP_GPP_casia_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_casia_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_shyear{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, GPP_casia_monthly_mean(idx, i));
    CP_GPP_casia_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_casia_monthly(idx, i, j));
        CP_GPP_casia_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_casia_annual(idx, j));
            CP_GPP_casia_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_casia_shyear(idx, j));
            CP_GPP_casia_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Eastern US
EP_GPP_eastus_monthly_beta = NaN(12, length(models));
EP_GPP_eastus_monthly_mean_beta = NaN(12, 1);
EP_GPP_eastus_annual_beta = NaN(1, length(models));
EP_GPP_eastus_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_eastus_annual_mean(idx));
EP_GPP_eastus_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_eastus_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_annual{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, GPP_eastus_shyear_mean(idx));
EP_GPP_eastus_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_eastus_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_shyear{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, GPP_eastus_monthly_mean(idx, i));
    EP_GPP_eastus_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_eastus_monthly(idx, i, j));
        EP_GPP_eastus_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_eastus_annual(idx, j));
            EP_GPP_eastus_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_eastus_shyear(idx, j));
            EP_GPP_eastus_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_GPP_eastus_monthly_beta = NaN(12, length(models));
CP_GPP_eastus_monthly_mean_beta = NaN(12, 1);
CP_GPP_eastus_annual_beta = NaN(1, length(models));
CP_GPP_eastus_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_eastus_annual_mean(idx));
CP_GPP_eastus_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_eastus_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_annual{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, GPP_eastus_shyear_mean(idx));
CP_GPP_eastus_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_eastus_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_shyear{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, GPP_eastus_monthly_mean(idx, i));
    CP_GPP_eastus_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_eastus_monthly(idx, i, j));
        CP_GPP_eastus_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_eastus_annual(idx, j));
            CP_GPP_eastus_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_eastus_shyear(idx, j));
            CP_GPP_eastus_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Europe
EP_GPP_europe_monthly_beta = NaN(12, length(models));
EP_GPP_europe_monthly_mean_beta = NaN(12, 1);
EP_GPP_europe_annual_beta = NaN(1, length(models));
EP_GPP_europe_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_europe_annual_mean(idx));
EP_GPP_europe_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_europe_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_annual{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, GPP_europe_shyear_mean(idx));
EP_GPP_europe_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_europe_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_shyear{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, GPP_europe_monthly_mean(idx, i));
    EP_GPP_europe_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_europe_monthly(idx, i, j));
        EP_GPP_europe_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_europe_annual(idx, j));
            EP_GPP_europe_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_europe_shyear(idx, j));
            EP_GPP_europe_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_GPP_europe_monthly_beta = NaN(12, length(models));
CP_GPP_europe_monthly_mean_beta = NaN(12, 1);
CP_GPP_europe_annual_beta = NaN(1, length(models));
CP_GPP_europe_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_europe_annual_mean(idx));
CP_GPP_europe_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_europe_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_annual{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, GPP_europe_shyear_mean(idx));
CP_GPP_europe_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_europe_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_shyear{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, GPP_europe_monthly_mean(idx, i));
    CP_GPP_europe_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_europe_monthly(idx, i, j));
        CP_GPP_europe_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_europe_annual(idx, j));
            CP_GPP_europe_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_europe_shyear(idx, j));
            CP_GPP_europe_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% The Sahel
EP_GPP_sahel_monthly_beta = NaN(12, length(models));
EP_GPP_sahel_monthly_mean_beta = NaN(12, 1);
EP_GPP_sahel_annual_beta = NaN(1, length(models));
EP_GPP_sahel_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_sahel_annual_mean(idx));
EP_GPP_sahel_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_sahel_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_annual{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, GPP_sahel_shyear_mean(idx));
EP_GPP_sahel_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_sahel_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_shyear{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, GPP_sahel_monthly_mean(idx, i));
    EP_GPP_sahel_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_sahel_monthly(idx, i, j));
        EP_GPP_sahel_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_sahel_annual(idx, j));
            EP_GPP_sahel_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_sahel_shyear(idx, j));
            EP_GPP_sahel_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_GPP_sahel_monthly_beta = NaN(12, length(models));
CP_GPP_sahel_monthly_mean_beta = NaN(12, 1);
CP_GPP_sahel_annual_beta = NaN(1, length(models));
CP_GPP_sahel_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_sahel_annual_mean(idx));
CP_GPP_sahel_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_sahel_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_annual{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, GPP_sahel_shyear_mean(idx));
CP_GPP_sahel_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_sahel_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_shyear{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, GPP_sahel_monthly_mean(idx, i));
    CP_GPP_sahel_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_sahel_monthly(idx, i, j));
        CP_GPP_sahel_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_sahel_annual(idx, j));
            CP_GPP_sahel_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_sahel_shyear(idx, j));
            CP_GPP_sahel_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Western North America
EP_GPP_westna_monthly_beta = NaN(12, length(models));
EP_GPP_westna_monthly_mean_beta = NaN(12, 1);
EP_GPP_westna_annual_beta = NaN(1, length(models));
EP_GPP_westna_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_westna_annual_mean(idx));
EP_GPP_westna_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_westna_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_annual{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, GPP_westna_shyear_mean(idx));
EP_GPP_westna_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_westna_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_shyear{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, GPP_westna_monthly_mean(idx, i));
    EP_GPP_westna_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_westna_monthly(idx, i, j));
        EP_GPP_westna_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_westna_annual(idx, j));
            EP_GPP_westna_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_westna_shyear(idx, j));
            EP_GPP_westna_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_GPP_westna_monthly_beta = NaN(12, length(models));
CP_GPP_westna_monthly_mean_beta = NaN(12, 1);
CP_GPP_westna_annual_beta = NaN(1, length(models));
CP_GPP_westna_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_westna_annual_mean(idx));
CP_GPP_westna_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_westna_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_annual{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, GPP_westna_shyear_mean(idx));
CP_GPP_westna_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_westna_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_shyear{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, GPP_westna_monthly_mean(idx, i));
    CP_GPP_westna_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_westna_monthly(idx, i, j));
        CP_GPP_westna_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_westna_annual(idx, j));
            CP_GPP_westna_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_westna_shyear(idx, j));
            CP_GPP_westna_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Tropical forest
EP_GPP_tropical_monthly_beta = NaN(12, length(models));
EP_GPP_tropical_monthly_mean_beta = NaN(12, 1);
EP_GPP_tropical_annual_beta = NaN(1, length(models));
EP_GPP_tropical_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_tropical_annual_mean(idx));
EP_GPP_tropical_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_tropical_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_annual{9,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, GPP_tropical_shyear_mean(idx));
EP_GPP_tropical_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_tropical_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_shyear{9,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, GPP_tropical_monthly_mean(idx, i));
    EP_GPP_tropical_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_tropical_monthly(idx, i, j));
        EP_GPP_tropical_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_tropical_annual(idx, j));
            EP_GPP_tropical_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_tropical_shyear(idx, j));
            EP_GPP_tropical_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_GPP_tropical_monthly_beta = NaN(12, length(models));
CP_GPP_tropical_monthly_mean_beta = NaN(12, 1);
CP_GPP_tropical_annual_beta = NaN(1, length(models));
CP_GPP_tropical_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_tropical_annual_mean(idx));
CP_GPP_tropical_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_tropical_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_annual{9,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, GPP_tropical_shyear_mean(idx));
CP_GPP_tropical_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_tropical_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_shyear{9,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, GPP_tropical_monthly_mean(idx, i));
    CP_GPP_tropical_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_tropical_monthly(idx, i, j));
        CP_GPP_tropical_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_tropical_annual(idx, j));
            CP_GPP_tropical_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_tropical_shyear(idx, j));
            CP_GPP_tropical_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Extratropical forest
EP_GPP_extratropical_monthly_beta = NaN(12, length(models));
EP_GPP_extratropical_monthly_mean_beta = NaN(12, 1);
EP_GPP_extratropical_annual_beta = NaN(1, length(models));
EP_GPP_extratropical_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_extratropical_annual_mean(idx));
EP_GPP_extratropical_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_extratropical_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_annual{10,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, GPP_extratropical_shyear_mean(idx));
EP_GPP_extratropical_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_extratropical_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_shyear{10,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, GPP_extratropical_monthly_mean(idx, i));
    EP_GPP_extratropical_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_extratropical_monthly(idx, i, j));
        EP_GPP_extratropical_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_extratropical_annual(idx, j));
            EP_GPP_extratropical_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_extratropical_shyear(idx, j));
            EP_GPP_extratropical_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_GPP_extratropical_monthly_beta = NaN(12, length(models));
CP_GPP_extratropical_monthly_mean_beta = NaN(12, 1);
CP_GPP_extratropical_annual_beta = NaN(1, length(models));
CP_GPP_extratropical_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_extratropical_annual_mean(idx));
CP_GPP_extratropical_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_extratropical_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_annual{10,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, GPP_extratropical_shyear_mean(idx));
CP_GPP_extratropical_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_extratropical_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_shyear{10,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, GPP_extratropical_monthly_mean(idx, i));
    CP_GPP_extratropical_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_extratropical_monthly(idx, i, j));
        CP_GPP_extratropical_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_extratropical_annual(idx, j));
            CP_GPP_extratropical_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_extratropical_shyear(idx, j));
            CP_GPP_extratropical_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Tundra/Arctic shrub
EP_GPP_tundra_monthly_beta = NaN(12, length(models));
EP_GPP_tundra_monthly_mean_beta = NaN(12, 1);
EP_GPP_tundra_annual_beta = NaN(1, length(models));
EP_GPP_tundra_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_tundra_annual_mean(idx));
EP_GPP_tundra_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_tundra_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_annual{11,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, GPP_tundra_shyear_mean(idx));
EP_GPP_tundra_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_tundra_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_shyear{11,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, GPP_tundra_monthly_mean(idx, i));
    EP_GPP_tundra_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_tundra_monthly(idx, i, j));
        EP_GPP_tundra_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_tundra_annual(idx, j));
            EP_GPP_tundra_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_tundra_shyear(idx, j));
            EP_GPP_tundra_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_GPP_tundra_monthly_beta = NaN(12, length(models));
CP_GPP_tundra_monthly_mean_beta = NaN(12, 1);
CP_GPP_tundra_annual_beta = NaN(1, length(models));
CP_GPP_tundra_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_tundra_annual_mean(idx));
CP_GPP_tundra_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_tundra_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_annual{11,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, GPP_tundra_shyear_mean(idx));
CP_GPP_tundra_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_tundra_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_shyear{11,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, GPP_tundra_monthly_mean(idx, i));
    CP_GPP_tundra_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_tundra_monthly(idx, i, j));
        CP_GPP_tundra_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_tundra_annual(idx, j));
            CP_GPP_tundra_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_tundra_shyear(idx, j));
            CP_GPP_tundra_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Grass/crop
EP_GPP_grass_monthly_beta = NaN(12, length(models));
EP_GPP_grass_monthly_mean_beta = NaN(12, 1);
EP_GPP_grass_annual_beta = NaN(1, length(models));
EP_GPP_grass_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_grass_annual_mean(idx));
EP_GPP_grass_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_grass_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_annual{12,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, GPP_grass_shyear_mean(idx));
EP_GPP_grass_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_grass_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_shyear{12,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, GPP_grass_monthly_mean(idx, i));
    EP_GPP_grass_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_grass_monthly(idx, i, j));
        EP_GPP_grass_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_grass_annual(idx, j));
            EP_GPP_grass_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_grass_shyear(idx, j));
            EP_GPP_grass_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_GPP_grass_monthly_beta = NaN(12, length(models));
CP_GPP_grass_monthly_mean_beta = NaN(12, 1);
CP_GPP_grass_annual_beta = NaN(1, length(models));
CP_GPP_grass_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_grass_annual_mean(idx));
CP_GPP_grass_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_grass_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_annual{12,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, GPP_grass_shyear_mean(idx));
CP_GPP_grass_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_grass_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_shyear{12,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, GPP_grass_monthly_mean(idx, i));
    CP_GPP_grass_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_grass_monthly(idx, i, j));
        CP_GPP_grass_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_grass_annual(idx, j));
            CP_GPP_grass_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_grass_shyear(idx, j));
            CP_GPP_grass_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Semiarid
EP_GPP_semiarid_monthly_beta = NaN(12, length(models));
EP_GPP_semiarid_monthly_mean_beta = NaN(12, 1);
EP_GPP_semiarid_annual_beta = NaN(1, length(models));
EP_GPP_semiarid_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, GPP_semiarid_annual_mean(idx));
EP_GPP_semiarid_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_semiarid_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_annual{13,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, GPP_semiarid_shyear_mean(idx));
EP_GPP_semiarid_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_GPP_semiarid_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_gpp_95CI_shyear{13,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, GPP_semiarid_monthly_mean(idx, i));
    EP_GPP_semiarid_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, GPP_semiarid_monthly(idx, i, j));
        EP_GPP_semiarid_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, GPP_semiarid_annual(idx, j));
            EP_GPP_semiarid_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, GPP_semiarid_shyear(idx, j));
            EP_GPP_semiarid_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_GPP_semiarid_monthly_beta = NaN(12, length(models));
CP_GPP_semiarid_monthly_mean_beta = NaN(12, 1);
CP_GPP_semiarid_annual_beta = NaN(1, length(models));
CP_GPP_semiarid_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, GPP_semiarid_annual_mean(idx));
CP_GPP_semiarid_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_semiarid_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_annual{13,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, GPP_semiarid_shyear_mean(idx));
CP_GPP_semiarid_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_GPP_semiarid_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_gpp_95CI_shyear{13,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, GPP_semiarid_monthly_mean(idx, i));
    CP_GPP_semiarid_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, GPP_semiarid_monthly(idx, i, j));
        CP_GPP_semiarid_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, GPP_semiarid_annual(idx, j));
            CP_GPP_semiarid_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, GPP_semiarid_shyear(idx, j));
            CP_GPP_semiarid_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;


save('./data/cp_ep_gpp_mstmip_regional.mat');

