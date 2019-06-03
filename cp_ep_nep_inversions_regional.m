% Examine relationships between modeled MsTMIP GPP and EPI/CPI over 
% 1951-2010 period
syear = 1979;
eyear = 2013;

models = {'v15r2','v15r4','v13r1','v14r2'}; nm = length(models);

%% Regress monthly & annual gpp against CPI and EPI indices
load ./data/nep_inversions_regional.mat;
load ./data/cpi_epi_1951-2016.mat;

ep_nep_95CI = cell2table(cell(8,1), 'VariableNames',{'Inversions'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US','Europe','The Sahel','Western North America'});
cp_nep_95CI = cell2table(cell(8,1), 'VariableNames',{'Inversions'},...
    'RowNames',{'Africa','Amazonia','Australia','Central Asia','Eastern US','Europe','The Sahel','Western North America'});

[ny,nx,~] = size(NEP_annual_mean);

cpi = cpi(yr>=syear & yr<=eyear);
epi = epi(yr>=syear & yr<=eyear);
yr = yr(yr>=syear & yr<=eyear);

idx = years>=syear & years<=eyear;

% Africa
EP_NEP_africa_monthly_beta = NaN(12, length(models));
EP_NEP_africa_monthly_mean_beta = NaN(12, 1);
EP_NEP_africa_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_africa_annual_mean(idx));
EP_NEP_africa_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_africa_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, NEP_africa_monthly_mean(idx, i));
    EP_NEP_africa_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_africa_monthly(idx, i, j));
        EP_NEP_africa_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_africa_annual(idx, j));
            EP_NEP_africa_annual_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_NEP_africa_monthly_beta = NaN(12, length(models));
CP_NEP_africa_monthly_mean_beta = NaN(12, 1);
CP_NEP_africa_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_africa_annual_mean(idx));
CP_NEP_africa_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_africa_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI{1,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, NEP_africa_monthly_mean(idx, i));
    CP_NEP_africa_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_africa_monthly(idx, i, j));
        CP_NEP_africa_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_africa_annual(idx, j));
            CP_NEP_africa_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Amazonia
EP_NEP_amazon_monthly_beta = NaN(12, length(models));
EP_NEP_amazon_monthly_mean_beta = NaN(12, 1);
EP_NEP_amazon_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_amazon_annual_mean(idx));
EP_NEP_amazon_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_amazon_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, NEP_amazon_monthly_mean(idx, i));
    EP_NEP_amazon_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_amazon_monthly(idx, i, j));
        EP_NEP_amazon_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_amazon_annual(idx, j));
            EP_NEP_amazon_annual_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_NEP_amazon_monthly_beta = NaN(12, length(models));
CP_NEP_amazon_monthly_mean_beta = NaN(12, 1);
CP_NEP_amazon_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_amazon_annual_mean(idx));
CP_NEP_amazon_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_amazon_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI{2,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, NEP_amazon_monthly_mean(idx, i));
    CP_NEP_amazon_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_amazon_monthly(idx, i, j));
        CP_NEP_amazon_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_amazon_annual(idx, j));
            CP_NEP_amazon_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Australia
EP_NEP_austr_monthly_beta = NaN(12, length(models));
EP_NEP_austr_monthly_mean_beta = NaN(12, 1);
EP_NEP_austr_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_austr_annual_mean(idx));
EP_NEP_austr_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_austr_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, NEP_austr_monthly_mean(idx, i));
    EP_NEP_austr_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_austr_monthly(idx, i, j));
        EP_NEP_austr_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_austr_annual(idx, j));
            EP_NEP_austr_annual_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_NEP_austr_monthly_beta = NaN(12, length(models));
CP_NEP_austr_monthly_mean_beta = NaN(12, 1);
CP_NEP_austr_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_austr_annual_mean(idx));
CP_NEP_austr_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_austr_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI{3,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, NEP_austr_monthly_mean(idx, i));
    CP_NEP_austr_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_austr_monthly(idx, i, j));
        CP_NEP_austr_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_austr_annual(idx, j));
            CP_NEP_austr_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Central Asia
EP_NEP_casia_monthly_beta = NaN(12, length(models));
EP_NEP_casia_monthly_mean_beta = NaN(12, 1);
EP_NEP_casia_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_casia_annual_mean(idx));
EP_NEP_casia_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_casia_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, NEP_casia_monthly_mean(idx, i));
    EP_NEP_casia_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_casia_monthly(idx, i, j));
        EP_NEP_casia_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_casia_annual(idx, j));
            EP_NEP_casia_annual_beta(j) = mdl.Coefficients.Estimate(2);
            
        end
    end
    
end

CP_NEP_casia_monthly_beta = NaN(12, length(models));
CP_NEP_casia_monthly_mean_beta = NaN(12, 1);
CP_NEP_casia_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_casia_annual_mean(idx));
CP_NEP_casia_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_casia_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI{4,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, NEP_casia_monthly_mean(idx, i));
    CP_NEP_casia_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_casia_monthly(idx, i, j));
        CP_NEP_casia_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_casia_annual(idx, j));
            CP_NEP_casia_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Eastern US
EP_NEP_eastus_monthly_beta = NaN(12, length(models));
EP_NEP_eastus_monthly_mean_beta = NaN(12, 1);
EP_NEP_eastus_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_eastus_annual_mean(idx));
EP_NEP_eastus_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_eastus_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, NEP_eastus_monthly_mean(idx, i));
    EP_NEP_eastus_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_eastus_monthly(idx, i, j));
        EP_NEP_eastus_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_eastus_annual(idx, j));
            EP_NEP_eastus_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_NEP_eastus_monthly_beta = NaN(12, length(models));
CP_NEP_eastus_monthly_mean_beta = NaN(12, 1);
CP_NEP_eastus_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_eastus_annual_mean(idx));
CP_NEP_eastus_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_eastus_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI{5,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, NEP_eastus_monthly_mean(idx, i));
    CP_NEP_eastus_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_eastus_monthly(idx, i, j));
        CP_NEP_eastus_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_eastus_annual(idx, j));
            CP_NEP_eastus_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Europe
EP_NEP_europe_monthly_beta = NaN(12, length(models));
EP_NEP_europe_monthly_mean_beta = NaN(12, 1);
EP_NEP_europe_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_europe_annual_mean(idx));
EP_NEP_europe_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_europe_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, NEP_europe_monthly_mean(idx, i));
    EP_NEP_europe_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_europe_monthly(idx, i, j));
        EP_NEP_europe_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_europe_annual(idx, j));
            EP_NEP_europe_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_NEP_europe_monthly_beta = NaN(12, length(models));
CP_NEP_europe_monthly_mean_beta = NaN(12, 1);
CP_NEP_europe_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_europe_annual_mean(idx));
CP_NEP_europe_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_europe_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI{6,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, NEP_europe_monthly_mean(idx, i));
    CP_NEP_europe_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_europe_monthly(idx, i, j));
        CP_NEP_europe_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_europe_annual(idx, j));
            CP_NEP_europe_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% The Sahel
EP_NEP_sahel_monthly_beta = NaN(12, length(models));
EP_NEP_sahel_monthly_mean_beta = NaN(12, 1);
EP_NEP_sahel_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_sahel_annual_mean(idx));
EP_NEP_sahel_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_sahel_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, NEP_sahel_monthly_mean(idx, i));
    EP_NEP_sahel_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_sahel_monthly(idx, i, j));
        EP_NEP_sahel_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_sahel_annual(idx, j));
            EP_NEP_sahel_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_NEP_sahel_monthly_beta = NaN(12, length(models));
CP_NEP_sahel_monthly_mean_beta = NaN(12, 1);
CP_NEP_sahel_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_sahel_annual_mean(idx));
CP_NEP_sahel_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_sahel_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI{7,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, NEP_sahel_monthly_mean(idx, i));
    CP_NEP_sahel_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_sahel_monthly(idx, i, j));
        CP_NEP_sahel_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_sahel_annual(idx, j));
            CP_NEP_sahel_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl;

% Western North America
EP_NEP_westna_monthly_beta = NaN(12, length(models));
EP_NEP_westna_monthly_mean_beta = NaN(12, 1);
EP_NEP_westna_annual_beta = NaN(1, length(models));

mdl = fitlm(epi, NEP_westna_annual_mean(idx));
EP_NEP_westna_annual_mean_beta = mdl.Coefficients.Estimate(2);
EP_NEP_westna_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
ep_nep_95CI{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(epi, NEP_westna_monthly_mean(idx, i));
    EP_NEP_westna_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(epi, NEP_westna_monthly(idx, i, j));
        EP_NEP_westna_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(epi, NEP_westna_annual(idx, j));
            EP_NEP_westna_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

CP_NEP_westna_monthly_beta = NaN(12, length(models));
CP_NEP_westna_monthly_mean_beta = NaN(12, 1);
CP_NEP_westna_annual_beta = NaN(1, length(models));

mdl = fitlm(cpi, NEP_westna_annual_mean(idx));
CP_NEP_westna_annual_mean_beta = mdl.Coefficients.Estimate(2);
CP_NEP_westna_annual_mean_beta_CI = 1.96*mdl.Coefficients.SE(2);
cp_nep_95CI{8,1} = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
    mdl.Coefficients.Estimate(2)/1000, 1.96*mdl.Coefficients.SE(2)/1000));
for i = 1:12
    mdl = fitlm(cpi, NEP_westna_monthly_mean(idx, i));
    CP_NEP_westna_monthly_mean_beta(i) = mdl.Coefficients.Estimate(2);
    for j = 1:length(models)
        mdl = fitlm(cpi, NEP_westna_monthly(idx, i, j));
        CP_NEP_westna_monthly_beta(i, j) = mdl.Coefficients.Estimate(2);
        if i == 1
            mdl = fitlm(cpi, NEP_westna_annual(idx, j));
            CP_NEP_westna_annual_beta(j) = mdl.Coefficients.Estimate(2);
        end
    end
    
end

clear i j mdl NEP* lsf;


save('./data/cp_ep_nep_inversions_regional.mat');

