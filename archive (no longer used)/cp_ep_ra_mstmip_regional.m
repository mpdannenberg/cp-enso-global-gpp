% Examine relationships between modeled MsTMIP GPP and EPI/CPI over 
% 1951-2010 period

syear = 1951;
eyear = 2010;

%% Regress monthly & annual gpp against CPI and EPI indices
load ./data/ra_mstmip_regional.mat;
load ./data/cpi_epi_1951-2016.mat;

ep_ra_95CI = cell2table(cell(13,1), 'VariableNames',{'MsTMIP'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US',...
    'Europe','The Sahel','Western North America','Tropical Forests',...
    'Extratropical Forests','Tundra/Arctic Shrubland','Grass/Crops',...
    'Semiarid'});
cp_ra_95CI = cell2table(cell(13,1), 'VariableNames',{'MsTMIP'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US',...
    'Europe','The Sahel','Western North America','Tropical Forests',...
    'Extratropical Forests','Tundra/Arctic Shrubland','Grass/Crops',...
    'Semiarid'});

[ny,nx,~] = size(Ra_annual_mean);

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

idx = years>=syear & years<=eyear;

% Africa
EP_Ra_africa_monthly_beta = NaN(12, length(models));
EP_Ra_africa_monthly_mean_beta = NaN(12, 1);
EP_Ra_africa_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Ra_africa_annual_mean(idx));
EP_Ra_africa_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Ra_africa_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_ra_95CI{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Ra_africa_monthly_mean(idx, i));
    EP_Ra_africa_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Ra_africa_monthly(idx, i, j));
        EP_Ra_africa_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Ra_africa_annual(idx, j));
            EP_Ra_africa_annual_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_Ra_africa_monthly_beta = NaN(12, length(models));
CP_Ra_africa_monthly_mean_beta = NaN(12, 1);
CP_Ra_africa_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Ra_africa_annual_mean(idx));
CP_Ra_africa_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Ra_africa_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_ra_95CI{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Ra_africa_monthly_mean(idx, i));
    CP_Ra_africa_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Ra_africa_monthly(idx, i, j));
        CP_Ra_africa_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Ra_africa_annual(idx, j));
            CP_Ra_africa_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Amazonia
EP_Ra_amazon_monthly_beta = NaN(12, length(models));
EP_Ra_amazon_monthly_mean_beta = NaN(12, 1);
EP_Ra_amazon_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Ra_amazon_annual_mean(idx));
EP_Ra_amazon_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Ra_amazon_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_ra_95CI{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Ra_amazon_monthly_mean(idx, i));
    EP_Ra_amazon_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Ra_amazon_monthly(idx, i, j));
        EP_Ra_amazon_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Ra_amazon_annual(idx, j));
            EP_Ra_amazon_annual_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_Ra_amazon_monthly_beta = NaN(12, length(models));
CP_Ra_amazon_monthly_mean_beta = NaN(12, 1);
CP_Ra_amazon_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Ra_amazon_annual_mean(idx));
CP_Ra_amazon_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Ra_amazon_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_ra_95CI{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Ra_amazon_monthly_mean(idx, i));
    CP_Ra_amazon_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Ra_amazon_monthly(idx, i, j));
        CP_Ra_amazon_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Ra_amazon_annual(idx, j));
            CP_Ra_amazon_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Australia
EP_Ra_austr_monthly_beta = NaN(12, length(models));
EP_Ra_austr_monthly_mean_beta = NaN(12, 1);
EP_Ra_austr_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Ra_austr_annual_mean(idx));
EP_Ra_austr_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Ra_austr_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_ra_95CI{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Ra_austr_monthly_mean(idx, i));
    EP_Ra_austr_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Ra_austr_monthly(idx, i, j));
        EP_Ra_austr_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Ra_austr_annual(idx, j));
            EP_Ra_austr_annual_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_Ra_austr_monthly_beta = NaN(12, length(models));
CP_Ra_austr_monthly_mean_beta = NaN(12, 1);
CP_Ra_austr_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Ra_austr_annual_mean(idx));
CP_Ra_austr_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Ra_austr_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_ra_95CI{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Ra_austr_monthly_mean(idx, i));
    CP_Ra_austr_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Ra_austr_monthly(idx, i, j));
        CP_Ra_austr_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Ra_austr_annual(idx, j));
            CP_Ra_austr_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Central Asia
EP_Ra_casia_monthly_beta = NaN(12, length(models));
EP_Ra_casia_monthly_mean_beta = NaN(12, 1);
EP_Ra_casia_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Ra_casia_annual_mean(idx));
EP_Ra_casia_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Ra_casia_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_ra_95CI{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Ra_casia_monthly_mean(idx, i));
    EP_Ra_casia_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Ra_casia_monthly(idx, i, j));
        EP_Ra_casia_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Ra_casia_annual(idx, j));
            EP_Ra_casia_annual_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_Ra_casia_monthly_beta = NaN(12, length(models));
CP_Ra_casia_monthly_mean_beta = NaN(12, 1);
CP_Ra_casia_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Ra_casia_annual_mean(idx));
CP_Ra_casia_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Ra_casia_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_ra_95CI{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Ra_casia_monthly_mean(idx, i));
    CP_Ra_casia_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Ra_casia_monthly(idx, i, j));
        CP_Ra_casia_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Ra_casia_annual(idx, j));
            CP_Ra_casia_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Eastern US
EP_Ra_eastus_monthly_beta = NaN(12, length(models));
EP_Ra_eastus_monthly_mean_beta = NaN(12, 1);
EP_Ra_eastus_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Ra_eastus_annual_mean(idx));
EP_Ra_eastus_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Ra_eastus_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_ra_95CI{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Ra_eastus_monthly_mean(idx, i));
    EP_Ra_eastus_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Ra_eastus_monthly(idx, i, j));
        EP_Ra_eastus_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Ra_eastus_annual(idx, j));
            EP_Ra_eastus_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Ra_eastus_monthly_beta = NaN(12, length(models));
CP_Ra_eastus_monthly_mean_beta = NaN(12, 1);
CP_Ra_eastus_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Ra_eastus_annual_mean(idx));
CP_Ra_eastus_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Ra_eastus_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_ra_95CI{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Ra_eastus_monthly_mean(idx, i));
    CP_Ra_eastus_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Ra_eastus_monthly(idx, i, j));
        CP_Ra_eastus_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Ra_eastus_annual(idx, j));
            CP_Ra_eastus_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Europe
EP_Ra_europe_monthly_beta = NaN(12, length(models));
EP_Ra_europe_monthly_mean_beta = NaN(12, 1);
EP_Ra_europe_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Ra_europe_annual_mean(idx));
EP_Ra_europe_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Ra_europe_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_ra_95CI{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Ra_europe_monthly_mean(idx, i));
    EP_Ra_europe_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Ra_europe_monthly(idx, i, j));
        EP_Ra_europe_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Ra_europe_annual(idx, j));
            EP_Ra_europe_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Ra_europe_monthly_beta = NaN(12, length(models));
CP_Ra_europe_monthly_mean_beta = NaN(12, 1);
CP_Ra_europe_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Ra_europe_annual_mean(idx));
CP_Ra_europe_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Ra_europe_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_ra_95CI{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Ra_europe_monthly_mean(idx, i));
    CP_Ra_europe_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Ra_europe_monthly(idx, i, j));
        CP_Ra_europe_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Ra_europe_annual(idx, j));
            CP_Ra_europe_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% The Sahel
EP_Ra_sahel_monthly_beta = NaN(12, length(models));
EP_Ra_sahel_monthly_mean_beta = NaN(12, 1);
EP_Ra_sahel_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Ra_sahel_annual_mean(idx));
EP_Ra_sahel_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Ra_sahel_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_ra_95CI{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Ra_sahel_monthly_mean(idx, i));
    EP_Ra_sahel_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Ra_sahel_monthly(idx, i, j));
        EP_Ra_sahel_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Ra_sahel_annual(idx, j));
            EP_Ra_sahel_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Ra_sahel_monthly_beta = NaN(12, length(models));
CP_Ra_sahel_monthly_mean_beta = NaN(12, 1);
CP_Ra_sahel_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Ra_sahel_annual_mean(idx));
CP_Ra_sahel_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Ra_sahel_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_ra_95CI{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Ra_sahel_monthly_mean(idx, i));
    CP_Ra_sahel_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Ra_sahel_monthly(idx, i, j));
        CP_Ra_sahel_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Ra_sahel_annual(idx, j));
            CP_Ra_sahel_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Western North America
EP_Ra_westna_monthly_beta = NaN(12, length(models));
EP_Ra_westna_monthly_mean_beta = NaN(12, 1);
EP_Ra_westna_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Ra_westna_annual_mean(idx));
EP_Ra_westna_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Ra_westna_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_ra_95CI{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Ra_westna_monthly_mean(idx, i));
    EP_Ra_westna_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Ra_westna_monthly(idx, i, j));
        EP_Ra_westna_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Ra_westna_annual(idx, j));
            EP_Ra_westna_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Ra_westna_monthly_beta = NaN(12, length(models));
CP_Ra_westna_monthly_mean_beta = NaN(12, 1);
CP_Ra_westna_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Ra_westna_annual_mean(idx));
CP_Ra_westna_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Ra_westna_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_ra_95CI{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Ra_westna_monthly_mean(idx, i));
    CP_Ra_westna_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Ra_westna_monthly(idx, i, j));
        CP_Ra_westna_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Ra_westna_annual(idx, j));
            CP_Ra_westna_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Tropical forest
EP_Ra_tropical_monthly_beta = NaN(12, length(models));
EP_Ra_tropical_monthly_mean_beta = NaN(12, 1);
EP_Ra_tropical_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Ra_tropical_annual_mean(idx));
EP_Ra_tropical_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Ra_tropical_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_ra_95CI{9,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Ra_tropical_monthly_mean(idx, i));
    EP_Ra_tropical_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Ra_tropical_monthly(idx, i, j));
        EP_Ra_tropical_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Ra_tropical_annual(idx, j));
            EP_Ra_tropical_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Ra_tropical_monthly_beta = NaN(12, length(models));
CP_Ra_tropical_monthly_mean_beta = NaN(12, 1);
CP_Ra_tropical_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Ra_tropical_annual_mean(idx));
CP_Ra_tropical_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Ra_tropical_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_ra_95CI{9,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Ra_tropical_monthly_mean(idx, i));
    CP_Ra_tropical_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Ra_tropical_monthly(idx, i, j));
        CP_Ra_tropical_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Ra_tropical_annual(idx, j));
            CP_Ra_tropical_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Extratropical forest
EP_Ra_extratropical_monthly_beta = NaN(12, length(models));
EP_Ra_extratropical_monthly_mean_beta = NaN(12, 1);
EP_Ra_extratropical_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Ra_extratropical_annual_mean(idx));
EP_Ra_extratropical_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Ra_extratropical_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_ra_95CI{10,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Ra_extratropical_monthly_mean(idx, i));
    EP_Ra_extratropical_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Ra_extratropical_monthly(idx, i, j));
        EP_Ra_extratropical_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Ra_extratropical_annual(idx, j));
            EP_Ra_extratropical_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Ra_extratropical_monthly_beta = NaN(12, length(models));
CP_Ra_extratropical_monthly_mean_beta = NaN(12, 1);
CP_Ra_extratropical_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Ra_extratropical_annual_mean(idx));
CP_Ra_extratropical_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Ra_extratropical_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_ra_95CI{10,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Ra_extratropical_monthly_mean(idx, i));
    CP_Ra_extratropical_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Ra_extratropical_monthly(idx, i, j));
        CP_Ra_extratropical_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Ra_extratropical_annual(idx, j));
            CP_Ra_extratropical_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Tundra/Arctic shrub
EP_Ra_tundra_monthly_beta = NaN(12, length(models));
EP_Ra_tundra_monthly_mean_beta = NaN(12, 1);
EP_Ra_tundra_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Ra_tundra_annual_mean(idx));
EP_Ra_tundra_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Ra_tundra_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_ra_95CI{11,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Ra_tundra_monthly_mean(idx, i));
    EP_Ra_tundra_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Ra_tundra_monthly(idx, i, j));
        EP_Ra_tundra_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Ra_tundra_annual(idx, j));
            EP_Ra_tundra_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Ra_tundra_monthly_beta = NaN(12, length(models));
CP_Ra_tundra_monthly_mean_beta = NaN(12, 1);
CP_Ra_tundra_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Ra_tundra_annual_mean(idx));
CP_Ra_tundra_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Ra_tundra_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_ra_95CI{11,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Ra_tundra_monthly_mean(idx, i));
    CP_Ra_tundra_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Ra_tundra_monthly(idx, i, j));
        CP_Ra_tundra_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Ra_tundra_annual(idx, j));
            CP_Ra_tundra_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Grass/crop
EP_Ra_grass_monthly_beta = NaN(12, length(models));
EP_Ra_grass_monthly_mean_beta = NaN(12, 1);
EP_Ra_grass_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Ra_grass_annual_mean(idx));
EP_Ra_grass_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Ra_grass_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_ra_95CI{12,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Ra_grass_monthly_mean(idx, i));
    EP_Ra_grass_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Ra_grass_monthly(idx, i, j));
        EP_Ra_grass_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Ra_grass_annual(idx, j));
            EP_Ra_grass_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Ra_grass_monthly_beta = NaN(12, length(models));
CP_Ra_grass_monthly_mean_beta = NaN(12, 1);
CP_Ra_grass_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Ra_grass_annual_mean(idx));
CP_Ra_grass_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Ra_grass_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_ra_95CI{12,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Ra_grass_monthly_mean(idx, i));
    CP_Ra_grass_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Ra_grass_monthly(idx, i, j));
        CP_Ra_grass_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Ra_grass_annual(idx, j));
            CP_Ra_grass_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Semiarid
EP_Ra_semiarid_monthly_beta = NaN(12, length(models));
EP_Ra_semiarid_monthly_mean_beta = NaN(12, 1);
EP_Ra_semiarid_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, Ra_semiarid_annual_mean(idx));
EP_Ra_semiarid_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_Ra_semiarid_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_ra_95CI{13,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, Ra_semiarid_monthly_mean(idx, i));
    EP_Ra_semiarid_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, Ra_semiarid_monthly(idx, i, j));
        EP_Ra_semiarid_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, Ra_semiarid_annual(idx, j));
            EP_Ra_semiarid_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_Ra_semiarid_monthly_beta = NaN(12, length(models));
CP_Ra_semiarid_monthly_mean_beta = NaN(12, 1);
CP_Ra_semiarid_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, Ra_semiarid_annual_mean(idx));
CP_Ra_semiarid_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_Ra_semiarid_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_ra_95CI{13,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, Ra_semiarid_monthly_mean(idx, i));
    CP_Ra_semiarid_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, Ra_semiarid_monthly(idx, i, j));
        CP_Ra_semiarid_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, Ra_semiarid_annual(idx, j));
            CP_Ra_semiarid_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl Ra*;


save('./data/cp_ep_ra_mstmip_regional.mat');

