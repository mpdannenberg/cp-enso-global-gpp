% Examine relationships between modeled MsTMIP GPP and EPI/CPI over 
% 1951-2010 period

syear = 1951;
eyear = 2010;

%% Regress monthly & annual gpp against CPI and EPI indices
load ./data/rh_mstmip_regional.mat;
load ./data/cpi_epi_1951-2016.mat;

ep_rh_95CI = cell2table(cell(13,1), 'VariableNames',{'MsTMIP'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US',...
    'Europe','The Sahel','Western North America','Tropical Forests',...
    'Extratropical Forests','Tundra/Arctic Shrubland','Grass/Crops',...
    'Semiarid'});
cp_rh_95CI = cell2table(cell(13,1), 'VariableNames',{'MsTMIP'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US',...
    'Europe','The Sahel','Western North America','Tropical Forests',...
    'Extratropical Forests','Tundra/Arctic Shrubland','Grass/Crops',...
    'Semiarid'});

[ny,nx,~] = size(Rh_annual_mean);

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

idx = years>=syear & years<=eyear;

% Africa
EP_Rh_africa_monthly_beta = NaN(12, length(models));
EP_Rh_africa_monthly_mean_beta = NaN(12, 1);
EP_Rh_africa_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_africa_annual_mean(idx));
EP_Rh_africa_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_africa_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_rh_95CI{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Rh_africa_monthly_mean(idx, i));
    EP_Rh_africa_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_africa_monthly(idx, i, j));
        EP_Rh_africa_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_africa_annual(idx, j));
            EP_Rh_africa_annual_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_Rh_africa_monthly_beta = NaN(12, length(models));
CP_Rh_africa_monthly_mean_beta = NaN(12, 1);
CP_Rh_africa_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_africa_annual_mean(idx));
CP_Rh_africa_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_africa_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_rh_95CI{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Rh_africa_monthly_mean(idx, i));
    CP_Rh_africa_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_africa_monthly(idx, i, j));
        CP_Rh_africa_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_africa_annual(idx, j));
            CP_Rh_africa_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Amazonia
EP_Rh_amazon_monthly_beta = NaN(12, length(models));
EP_Rh_amazon_monthly_mean_beta = NaN(12, 1);
EP_Rh_amazon_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_amazon_annual_mean(idx));
EP_Rh_amazon_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_amazon_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_rh_95CI{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Rh_amazon_monthly_mean(idx, i));
    EP_Rh_amazon_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_amazon_monthly(idx, i, j));
        EP_Rh_amazon_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_amazon_annual(idx, j));
            EP_Rh_amazon_annual_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_Rh_amazon_monthly_beta = NaN(12, length(models));
CP_Rh_amazon_monthly_mean_beta = NaN(12, 1);
CP_Rh_amazon_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_amazon_annual_mean(idx));
CP_Rh_amazon_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_amazon_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_rh_95CI{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Rh_amazon_monthly_mean(idx, i));
    CP_Rh_amazon_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_amazon_monthly(idx, i, j));
        CP_Rh_amazon_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_amazon_annual(idx, j));
            CP_Rh_amazon_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Australia
EP_Rh_austr_monthly_beta = NaN(12, length(models));
EP_Rh_austr_monthly_mean_beta = NaN(12, 1);
EP_Rh_austr_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_austr_annual_mean(idx));
EP_Rh_austr_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_austr_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_rh_95CI{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Rh_austr_monthly_mean(idx, i));
    EP_Rh_austr_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_austr_monthly(idx, i, j));
        EP_Rh_austr_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_austr_annual(idx, j));
            EP_Rh_austr_annual_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_Rh_austr_monthly_beta = NaN(12, length(models));
CP_Rh_austr_monthly_mean_beta = NaN(12, 1);
CP_Rh_austr_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_austr_annual_mean(idx));
CP_Rh_austr_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_austr_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_rh_95CI{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Rh_austr_monthly_mean(idx, i));
    CP_Rh_austr_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_austr_monthly(idx, i, j));
        CP_Rh_austr_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_austr_annual(idx, j));
            CP_Rh_austr_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Central Asia
EP_Rh_casia_monthly_beta = NaN(12, length(models));
EP_Rh_casia_monthly_mean_beta = NaN(12, 1);
EP_Rh_casia_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_casia_annual_mean(idx));
EP_Rh_casia_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_casia_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_rh_95CI{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Rh_casia_monthly_mean(idx, i));
    EP_Rh_casia_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_casia_monthly(idx, i, j));
        EP_Rh_casia_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_casia_annual(idx, j));
            EP_Rh_casia_annual_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_Rh_casia_monthly_beta = NaN(12, length(models));
CP_Rh_casia_monthly_mean_beta = NaN(12, 1);
CP_Rh_casia_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_casia_annual_mean(idx));
CP_Rh_casia_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_casia_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_rh_95CI{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Rh_casia_monthly_mean(idx, i));
    CP_Rh_casia_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_casia_monthly(idx, i, j));
        CP_Rh_casia_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_casia_annual(idx, j));
            CP_Rh_casia_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Eastern US
EP_Rh_eastus_monthly_beta = NaN(12, length(models));
EP_Rh_eastus_monthly_mean_beta = NaN(12, 1);
EP_Rh_eastus_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_eastus_annual_mean(idx));
EP_Rh_eastus_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_eastus_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_rh_95CI{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Rh_eastus_monthly_mean(idx, i));
    EP_Rh_eastus_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_eastus_monthly(idx, i, j));
        EP_Rh_eastus_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_eastus_annual(idx, j));
            EP_Rh_eastus_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Rh_eastus_monthly_beta = NaN(12, length(models));
CP_Rh_eastus_monthly_mean_beta = NaN(12, 1);
CP_Rh_eastus_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_eastus_annual_mean(idx));
CP_Rh_eastus_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_eastus_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_rh_95CI{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Rh_eastus_monthly_mean(idx, i));
    CP_Rh_eastus_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_eastus_monthly(idx, i, j));
        CP_Rh_eastus_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_eastus_annual(idx, j));
            CP_Rh_eastus_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Europe
EP_Rh_europe_monthly_beta = NaN(12, length(models));
EP_Rh_europe_monthly_mean_beta = NaN(12, 1);
EP_Rh_europe_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_europe_annual_mean(idx));
EP_Rh_europe_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_europe_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_rh_95CI{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Rh_europe_monthly_mean(idx, i));
    EP_Rh_europe_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_europe_monthly(idx, i, j));
        EP_Rh_europe_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_europe_annual(idx, j));
            EP_Rh_europe_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Rh_europe_monthly_beta = NaN(12, length(models));
CP_Rh_europe_monthly_mean_beta = NaN(12, 1);
CP_Rh_europe_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_europe_annual_mean(idx));
CP_Rh_europe_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_europe_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_rh_95CI{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Rh_europe_monthly_mean(idx, i));
    CP_Rh_europe_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_europe_monthly(idx, i, j));
        CP_Rh_europe_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_europe_annual(idx, j));
            CP_Rh_europe_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% The Sahel
EP_Rh_sahel_monthly_beta = NaN(12, length(models));
EP_Rh_sahel_monthly_mean_beta = NaN(12, 1);
EP_Rh_sahel_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_sahel_annual_mean(idx));
EP_Rh_sahel_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_sahel_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_rh_95CI{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Rh_sahel_monthly_mean(idx, i));
    EP_Rh_sahel_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_sahel_monthly(idx, i, j));
        EP_Rh_sahel_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_sahel_annual(idx, j));
            EP_Rh_sahel_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Rh_sahel_monthly_beta = NaN(12, length(models));
CP_Rh_sahel_monthly_mean_beta = NaN(12, 1);
CP_Rh_sahel_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_sahel_annual_mean(idx));
CP_Rh_sahel_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_sahel_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_rh_95CI{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Rh_sahel_monthly_mean(idx, i));
    CP_Rh_sahel_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_sahel_monthly(idx, i, j));
        CP_Rh_sahel_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_sahel_annual(idx, j));
            CP_Rh_sahel_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Western North America
EP_Rh_westna_monthly_beta = NaN(12, length(models));
EP_Rh_westna_monthly_mean_beta = NaN(12, 1);
EP_Rh_westna_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_westna_annual_mean(idx));
EP_Rh_westna_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_westna_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_rh_95CI{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Rh_westna_monthly_mean(idx, i));
    EP_Rh_westna_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_westna_monthly(idx, i, j));
        EP_Rh_westna_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_westna_annual(idx, j));
            EP_Rh_westna_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Rh_westna_monthly_beta = NaN(12, length(models));
CP_Rh_westna_monthly_mean_beta = NaN(12, 1);
CP_Rh_westna_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_westna_annual_mean(idx));
CP_Rh_westna_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_westna_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_rh_95CI{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Rh_westna_monthly_mean(idx, i));
    CP_Rh_westna_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_westna_monthly(idx, i, j));
        CP_Rh_westna_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_westna_annual(idx, j));
            CP_Rh_westna_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Tropical forest
EP_Rh_tropical_monthly_beta = NaN(12, length(models));
EP_Rh_tropical_monthly_mean_beta = NaN(12, 1);
EP_Rh_tropical_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_tropical_annual_mean(idx));
EP_Rh_tropical_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_tropical_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_rh_95CI{9,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Rh_tropical_monthly_mean(idx, i));
    EP_Rh_tropical_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_tropical_monthly(idx, i, j));
        EP_Rh_tropical_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_tropical_annual(idx, j));
            EP_Rh_tropical_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Rh_tropical_monthly_beta = NaN(12, length(models));
CP_Rh_tropical_monthly_mean_beta = NaN(12, 1);
CP_Rh_tropical_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_tropical_annual_mean(idx));
CP_Rh_tropical_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_tropical_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_rh_95CI{9,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Rh_tropical_monthly_mean(idx, i));
    CP_Rh_tropical_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_tropical_monthly(idx, i, j));
        CP_Rh_tropical_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_tropical_annual(idx, j));
            CP_Rh_tropical_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Extratropical forest
EP_Rh_extratropical_monthly_beta = NaN(12, length(models));
EP_Rh_extratropical_monthly_mean_beta = NaN(12, 1);
EP_Rh_extratropical_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_extratropical_annual_mean(idx));
EP_Rh_extratropical_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_extratropical_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_rh_95CI{10,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Rh_extratropical_monthly_mean(idx, i));
    EP_Rh_extratropical_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_extratropical_monthly(idx, i, j));
        EP_Rh_extratropical_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_extratropical_annual(idx, j));
            EP_Rh_extratropical_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Rh_extratropical_monthly_beta = NaN(12, length(models));
CP_Rh_extratropical_monthly_mean_beta = NaN(12, 1);
CP_Rh_extratropical_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_extratropical_annual_mean(idx));
CP_Rh_extratropical_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_extratropical_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_rh_95CI{10,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Rh_extratropical_monthly_mean(idx, i));
    CP_Rh_extratropical_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_extratropical_monthly(idx, i, j));
        CP_Rh_extratropical_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_extratropical_annual(idx, j));
            CP_Rh_extratropical_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Tundra/Arctic shrub
EP_Rh_tundra_monthly_beta = NaN(12, length(models));
EP_Rh_tundra_monthly_mean_beta = NaN(12, 1);
EP_Rh_tundra_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_tundra_annual_mean(idx));
EP_Rh_tundra_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_tundra_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_rh_95CI{11,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Rh_tundra_monthly_mean(idx, i));
    EP_Rh_tundra_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_tundra_monthly(idx, i, j));
        EP_Rh_tundra_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_tundra_annual(idx, j));
            EP_Rh_tundra_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Rh_tundra_monthly_beta = NaN(12, length(models));
CP_Rh_tundra_monthly_mean_beta = NaN(12, 1);
CP_Rh_tundra_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_tundra_annual_mean(idx));
CP_Rh_tundra_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_tundra_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_rh_95CI{11,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Rh_tundra_monthly_mean(idx, i));
    CP_Rh_tundra_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_tundra_monthly(idx, i, j));
        CP_Rh_tundra_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_tundra_annual(idx, j));
            CP_Rh_tundra_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Grass/crop
EP_Rh_grass_monthly_beta = NaN(12, length(models));
EP_Rh_grass_monthly_mean_beta = NaN(12, 1);
EP_Rh_grass_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_grass_annual_mean(idx));
EP_Rh_grass_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_grass_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_rh_95CI{12,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Rh_grass_monthly_mean(idx, i));
    EP_Rh_grass_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_grass_monthly(idx, i, j));
        EP_Rh_grass_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_grass_annual(idx, j));
            EP_Rh_grass_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Rh_grass_monthly_beta = NaN(12, length(models));
CP_Rh_grass_monthly_mean_beta = NaN(12, 1);
CP_Rh_grass_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_grass_annual_mean(idx));
CP_Rh_grass_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_grass_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_rh_95CI{12,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Rh_grass_monthly_mean(idx, i));
    CP_Rh_grass_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_grass_monthly(idx, i, j));
        CP_Rh_grass_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_grass_annual(idx, j));
            CP_Rh_grass_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Semiarid
EP_Rh_semiarid_monthly_beta = NaN(12, length(models));
EP_Rh_semiarid_monthly_mean_beta = NaN(12, 1);
EP_Rh_semiarid_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Rh_semiarid_annual_mean(idx));
EP_Rh_semiarid_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Rh_semiarid_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_rh_95CI{13,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Rh_semiarid_monthly_mean(idx, i));
    EP_Rh_semiarid_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Rh_semiarid_monthly(idx, i, j));
        EP_Rh_semiarid_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Rh_semiarid_annual(idx, j));
            EP_Rh_semiarid_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Rh_semiarid_monthly_beta = NaN(12, length(models));
CP_Rh_semiarid_monthly_mean_beta = NaN(12, 1);
CP_Rh_semiarid_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Rh_semiarid_annual_mean(idx));
CP_Rh_semiarid_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Rh_semiarid_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_rh_95CI{13,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Rh_semiarid_monthly_mean(idx, i));
    CP_Rh_semiarid_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Rh_semiarid_monthly(idx, i, j));
        CP_Rh_semiarid_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Rh_semiarid_annual(idx, j));
            CP_Rh_semiarid_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl Rh*;


save('./data/cp_ep_rh_mstmip_regional.mat');

