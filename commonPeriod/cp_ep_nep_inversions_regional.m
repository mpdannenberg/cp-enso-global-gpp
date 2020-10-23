% Examine relationships between CAMS NEP and EPI/CPI over 
% 1982-2010 common period
syear = 1982;
eyear = 2010;

models = {'v18r2'}; nm = length(models);

%% Regress monthly & annual nep against CPI and EPI indices
load ../data/nep_inversions_regional.mat;
load ../data/cpi_epi_1951-2016.mat;

ep_nep_95CI_annual = cell2table(cell(8,1), 'VariableNames',{'Inversions'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US','Europe','The Sahel','Western North America'});
cp_nep_95CI_annual = cell2table(cell(8,1), 'VariableNames',{'Inversions'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US','Europe','The Sahel','Western North America'});
ep_nep_95CI_shyear = cell2table(cell(8,1), 'VariableNames',{'Inversions'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US','Europe','The Sahel','Western North America'});
cp_nep_95CI_shyear = cell2table(cell(8,1), 'VariableNames',{'Inversions'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US','Europe','The Sahel','Western North America'});

[ny,nx,~] = size(NEP_annual_mean);

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

idx = years>=syear & years<=eyear;

% Africa
EP_NEP_africa_monthly_beta = NaN(18, length(models));
EP_NEP_africa_monthly_mean_beta = NaN(18, 1);
EP_NEP_africa_annual_beta = NaN(1, length(models));
EP_NEP_africa_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_africa_annual_mean(idx));
EP_NEP_africa_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_africa_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_annual{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, NEP_africa_shyear_mean(idx));
EP_NEP_africa_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_africa_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_shyear{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(epi, NEP_africa_monthly_mean(idx, i));
    EP_NEP_africa_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_africa_monthly(idx, i, j));
        EP_NEP_africa_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_africa_annual(idx, j));
            EP_NEP_africa_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, NEP_africa_shyear(idx, j));
            EP_NEP_africa_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_NEP_africa_monthly_beta = NaN(18, length(models));
CP_NEP_africa_monthly_mean_beta = NaN(18, 1);
CP_NEP_africa_annual_beta = NaN(1, length(models));
CP_NEP_africa_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_africa_annual_mean(idx));
CP_NEP_africa_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_africa_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_annual{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, NEP_africa_shyear_mean(idx));
CP_NEP_africa_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_africa_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_shyear{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(cpi, NEP_africa_monthly_mean(idx, i));
    CP_NEP_africa_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_africa_monthly(idx, i, j));
        CP_NEP_africa_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_africa_annual(idx, j));
            CP_NEP_africa_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, NEP_africa_shyear(idx, j));
            CP_NEP_africa_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Amazonia
EP_NEP_amazon_monthly_beta = NaN(18, length(models));
EP_NEP_amazon_monthly_mean_beta = NaN(18, 1);
EP_NEP_amazon_annual_beta = NaN(1, length(models));
EP_NEP_amazon_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_amazon_annual_mean(idx));
EP_NEP_amazon_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_amazon_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_annual{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, NEP_amazon_shyear_mean(idx));
EP_NEP_amazon_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_amazon_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_shyear{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(epi, NEP_amazon_monthly_mean(idx, i));
    EP_NEP_amazon_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_amazon_monthly(idx, i, j));
        EP_NEP_amazon_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_amazon_annual(idx, j));
            EP_NEP_amazon_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, NEP_amazon_shyear(idx, j));
            EP_NEP_amazon_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_NEP_amazon_monthly_beta = NaN(18, length(models));
CP_NEP_amazon_monthly_mean_beta = NaN(18, 1);
CP_NEP_amazon_annual_beta = NaN(1, length(models));
CP_NEP_amazon_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_amazon_annual_mean(idx));
CP_NEP_amazon_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_amazon_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_annual{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, NEP_amazon_shyear_mean(idx));
CP_NEP_amazon_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_amazon_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_shyear{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(cpi, NEP_amazon_monthly_mean(idx, i));
    CP_NEP_amazon_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_amazon_monthly(idx, i, j));
        CP_NEP_amazon_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_amazon_annual(idx, j));
            CP_NEP_amazon_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, NEP_amazon_shyear(idx, j));
            CP_NEP_amazon_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Australia
EP_NEP_austr_monthly_beta = NaN(18, length(models));
EP_NEP_austr_monthly_mean_beta = NaN(18, 1);
EP_NEP_austr_annual_beta = NaN(1, length(models));
EP_NEP_austr_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_austr_annual_mean(idx));
EP_NEP_austr_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_austr_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_annual{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, NEP_austr_shyear_mean(idx));
EP_NEP_austr_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_austr_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_shyear{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(epi, NEP_austr_monthly_mean(idx, i));
    EP_NEP_austr_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_austr_monthly(idx, i, j));
        EP_NEP_austr_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_austr_annual(idx, j));
            EP_NEP_austr_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, NEP_austr_shyear(idx, j));
            EP_NEP_austr_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_NEP_austr_monthly_beta = NaN(18, length(models));
CP_NEP_austr_monthly_mean_beta = NaN(18, 1);
CP_NEP_austr_annual_beta = NaN(1, length(models));
CP_NEP_austr_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_austr_annual_mean(idx));
CP_NEP_austr_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_austr_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_annual{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, NEP_austr_shyear_mean(idx));
CP_NEP_austr_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_austr_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_shyear{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(cpi, NEP_austr_monthly_mean(idx, i));
    CP_NEP_austr_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_austr_monthly(idx, i, j));
        CP_NEP_austr_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_austr_annual(idx, j));
            CP_NEP_austr_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, NEP_austr_shyear(idx, j));
            CP_NEP_austr_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Central Asia
EP_NEP_casia_monthly_beta = NaN(18, length(models));
EP_NEP_casia_monthly_mean_beta = NaN(18, 1);
EP_NEP_casia_annual_beta = NaN(1, length(models));
EP_NEP_casia_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_casia_annual_mean(idx));
EP_NEP_casia_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_casia_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_annual{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, NEP_casia_shyear_mean(idx));
EP_NEP_casia_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_casia_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_shyear{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(epi, NEP_casia_monthly_mean(idx, i));
    EP_NEP_casia_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_casia_monthly(idx, i, j));
        EP_NEP_casia_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_casia_annual(idx, j));
            EP_NEP_casia_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, NEP_casia_shyear(idx, j));
            EP_NEP_casia_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_NEP_casia_monthly_beta = NaN(18, length(models));
CP_NEP_casia_monthly_mean_beta = NaN(18, 1);
CP_NEP_casia_annual_beta = NaN(1, length(models));
CP_NEP_casia_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_casia_annual_mean(idx));
CP_NEP_casia_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_casia_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_annual{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, NEP_casia_shyear_mean(idx));
CP_NEP_casia_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_casia_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_shyear{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(cpi, NEP_casia_monthly_mean(idx, i));
    CP_NEP_casia_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_casia_monthly(idx, i, j));
        CP_NEP_casia_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_casia_annual(idx, j));
            CP_NEP_casia_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, NEP_casia_shyear(idx, j));
            CP_NEP_casia_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Eastern US
EP_NEP_eastus_monthly_beta = NaN(18, length(models));
EP_NEP_eastus_monthly_mean_beta = NaN(18, 1);
EP_NEP_eastus_annual_beta = NaN(1, length(models));
EP_NEP_eastus_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_eastus_annual_mean(idx));
EP_NEP_eastus_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_eastus_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_annual{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, NEP_eastus_shyear_mean(idx));
EP_NEP_eastus_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_eastus_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_shyear{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(epi, NEP_eastus_monthly_mean(idx, i));
    EP_NEP_eastus_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_eastus_monthly(idx, i, j));
        EP_NEP_eastus_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_eastus_annual(idx, j));
            EP_NEP_eastus_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, NEP_eastus_shyear(idx, j));
            EP_NEP_eastus_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_NEP_eastus_monthly_beta = NaN(18, length(models));
CP_NEP_eastus_monthly_mean_beta = NaN(18, 1);
CP_NEP_eastus_annual_beta = NaN(1, length(models));
CP_NEP_eastus_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_eastus_annual_mean(idx));
CP_NEP_eastus_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_eastus_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_annual{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, NEP_eastus_shyear_mean(idx));
CP_NEP_eastus_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_eastus_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_shyear{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(cpi, NEP_eastus_monthly_mean(idx, i));
    CP_NEP_eastus_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_eastus_monthly(idx, i, j));
        CP_NEP_eastus_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_eastus_annual(idx, j));
            CP_NEP_eastus_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, NEP_eastus_shyear(idx, j));
            CP_NEP_eastus_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Europe
EP_NEP_europe_monthly_beta = NaN(18, length(models));
EP_NEP_europe_monthly_mean_beta = NaN(18, 1);
EP_NEP_europe_annual_beta = NaN(1, length(models));
EP_NEP_europe_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_europe_annual_mean(idx));
EP_NEP_europe_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_europe_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_annual{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, NEP_europe_shyear_mean(idx));
EP_NEP_europe_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_europe_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_shyear{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(epi, NEP_europe_monthly_mean(idx, i));
    EP_NEP_europe_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_europe_monthly(idx, i, j));
        EP_NEP_europe_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_europe_annual(idx, j));
            EP_NEP_europe_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, NEP_europe_shyear(idx, j));
            EP_NEP_europe_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_NEP_europe_monthly_beta = NaN(18, length(models));
CP_NEP_europe_monthly_mean_beta = NaN(18, 1);
CP_NEP_europe_annual_beta = NaN(1, length(models));
CP_NEP_europe_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_europe_annual_mean(idx));
CP_NEP_europe_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_europe_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_annual{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, NEP_europe_shyear_mean(idx));
CP_NEP_europe_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_europe_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_shyear{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(cpi, NEP_europe_monthly_mean(idx, i));
    CP_NEP_europe_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_europe_monthly(idx, i, j));
        CP_NEP_europe_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_europe_annual(idx, j));
            CP_NEP_europe_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, NEP_europe_shyear(idx, j));
            CP_NEP_europe_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% The Sahel
EP_NEP_sahel_monthly_beta = NaN(18, length(models));
EP_NEP_sahel_monthly_mean_beta = NaN(18, 1);
EP_NEP_sahel_annual_beta = NaN(1, length(models));
EP_NEP_sahel_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_sahel_annual_mean(idx));
EP_NEP_sahel_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_sahel_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_annual{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, NEP_sahel_shyear_mean(idx));
EP_NEP_sahel_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_sahel_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_shyear{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(epi, NEP_sahel_monthly_mean(idx, i));
    EP_NEP_sahel_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_sahel_monthly(idx, i, j));
        EP_NEP_sahel_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_sahel_annual(idx, j));
            EP_NEP_sahel_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, NEP_sahel_shyear(idx, j));
            EP_NEP_sahel_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_NEP_sahel_monthly_beta = NaN(18, length(models));
CP_NEP_sahel_monthly_mean_beta = NaN(18, 1);
CP_NEP_sahel_annual_beta = NaN(1, length(models));
CP_NEP_sahel_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_sahel_annual_mean(idx));
CP_NEP_sahel_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_sahel_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_annual{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, NEP_sahel_shyear_mean(idx));
CP_NEP_sahel_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_sahel_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_shyear{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(cpi, NEP_sahel_monthly_mean(idx, i));
    CP_NEP_sahel_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_sahel_monthly(idx, i, j));
        CP_NEP_sahel_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_sahel_annual(idx, j));
            CP_NEP_sahel_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, NEP_sahel_shyear(idx, j));
            CP_NEP_sahel_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;

% Western North America
EP_NEP_westna_monthly_beta = NaN(18, length(models));
EP_NEP_westna_monthly_mean_beta = NaN(18, 1);
EP_NEP_westna_annual_beta = NaN(1, length(models));
EP_NEP_westna_shyear_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_westna_annual_mean(idx));
EP_NEP_westna_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_westna_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_annual{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(epi, NEP_westna_shyear_mean(idx));
EP_NEP_westna_shyear_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_westna_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI_shyear{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(epi, NEP_westna_monthly_mean(idx, i));
    EP_NEP_westna_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_westna_monthly(idx, i, j));
        EP_NEP_westna_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_westna_annual(idx, j));
            EP_NEP_westna_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(epi, NEP_westna_shyear(idx, j));
            EP_NEP_westna_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_NEP_westna_monthly_beta = NaN(18, length(models));
CP_NEP_westna_monthly_mean_beta = NaN(18, 1);
CP_NEP_westna_annual_beta = NaN(1, length(models));
CP_NEP_westna_shyear_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_westna_annual_mean(idx));
CP_NEP_westna_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_westna_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_annual{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
mdl = fitlm(cpi, NEP_westna_shyear_mean(idx));
CP_NEP_westna_shyear_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_westna_shyear_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI_shyear{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:18
    mdl = fitlm(cpi, NEP_westna_monthly_mean(idx, i));
    CP_NEP_westna_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_westna_monthly(idx, i, j));
        CP_NEP_westna_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_westna_annual(idx, j));
            CP_NEP_westna_annual_beta(j) = mdl.Coefficients.Estimate(2);
            mdl = fitlm(cpi, NEP_westna_shyear(idx, j));
            CP_NEP_westna_shyear_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

clear i j mdl;


save('./data/cp_ep_nep_inversions_regional.mat');

