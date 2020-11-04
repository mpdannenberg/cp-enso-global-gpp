% Build GPP and NEP tables

Tannual = table('Size',[15 4],...
    'RowNames',{'Global','Tropics','Tropical and S. Africa','Tropical South America','Australia','Tropical Asia','Eastern U.S.','Europe','The Sahel','W. North America','Tropical Forests','Extratropical forests','Tundra & Arctic Shrubland','Grass/Crops','Semiarid'},...
    'VariableNames',{'EP_MsTMIP','EP_CAMS','CP_MsTMIP','CP_CAMS'},...
    'VariableTypes',{'string','string','string','string'});
Tshyear = table('Size',[15 4],...
    'RowNames',{'Global','Tropics','Tropical and S. Africa','Tropical South America','Australia','Tropical Asia','Eastern U.S.','Europe','The Sahel','W. North America','Tropical Forests','Extratropical forests','Tundra & Arctic Shrubland','Grass/Crops','Semiarid'},...
    'VariableNames',{'EP_MsTMIP','EP_CAMS','CP_MsTMIP','CP_CAMS'},...
    'VariableTypes',{'string','string','string','string'});

% NEP global MsTMIP
load ./data/cp_ep_nep_mstmip.mat;
Tannual.EP_MsTMIP(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_NEP_global_annual_mean_beta/1000, EP_NEP_global_annual_mean_beta_CI/1000));
Tannual.EP_MsTMIP(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_NEP_tropics_annual_mean_beta/1000, EP_NEP_tropics_annual_mean_beta_CI/1000));

Tshyear.EP_MsTMIP(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_NEP_global_shyear_mean_beta/1000, EP_NEP_global_shyear_mean_beta_CI/1000));
Tshyear.EP_MsTMIP(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_NEP_tropics_shyear_mean_beta/1000, EP_NEP_tropics_shyear_mean_beta_CI/1000));

Tannual.CP_MsTMIP(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_NEP_global_annual_mean_beta/1000, CP_NEP_global_annual_mean_beta_CI/1000));
Tannual.CP_MsTMIP(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_NEP_tropics_annual_mean_beta/1000, CP_NEP_tropics_annual_mean_beta_CI/1000));

Tshyear.CP_MsTMIP(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_NEP_global_shyear_mean_beta/1000, CP_NEP_global_shyear_mean_beta_CI/1000));
Tshyear.CP_MsTMIP(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_NEP_tropics_shyear_mean_beta/1000, CP_NEP_tropics_shyear_mean_beta_CI/1000));
clearvars -except Tannual Tshyear

% NEP regional LUE
load ./data/cp_ep_nep_mstmip_regional.mat;
Tannual.EP_MsTMIP(3:end) = ep_nep_95CI_annual.MsTMIP;

Tshyear.EP_MsTMIP(3:end) = ep_nep_95CI_shyear.MsTMIP;

Tannual.CP_MsTMIP(3:end) = cp_nep_95CI_annual.MsTMIP;

Tshyear.CP_MsTMIP(3:end) = cp_nep_95CI_shyear.MsTMIP;
clearvars -except Tannual Tshyear

% NEP global inversion
load ./data/cp_ep_nep_inversions.mat;
Tannual.EP_CAMS(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_NEP_global_annual_mean_beta/1000, EP_NEP_global_annual_mean_beta_CI/1000));
Tannual.EP_CAMS(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_NEP_tropics_annual_mean_beta/1000, EP_NEP_tropics_annual_mean_beta_CI/1000));

Tshyear.EP_CAMS(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_NEP_global_shyear_mean_beta/1000, EP_NEP_global_shyear_mean_beta_CI/1000));
Tshyear.EP_CAMS(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_NEP_tropics_shyear_mean_beta/1000, EP_NEP_tropics_shyear_mean_beta_CI/1000));

Tannual.CP_CAMS(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_NEP_global_annual_mean_beta/1000, CP_NEP_global_annual_mean_beta_CI/1000));
Tannual.CP_CAMS(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_NEP_tropics_annual_mean_beta/1000, CP_NEP_tropics_annual_mean_beta_CI/1000));

Tshyear.CP_CAMS(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_NEP_global_shyear_mean_beta/1000, CP_NEP_global_shyear_mean_beta_CI/1000));
Tshyear.CP_CAMS(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_NEP_tropics_shyear_mean_beta/1000, CP_NEP_tropics_shyear_mean_beta_CI/1000));
clearvars -except Tannual Tshyear

% NEP regional inversion
load ./data/cp_ep_nep_inversions_regional.mat;
Tannual.EP_CAMS(3:10) = ep_nep_95CI_annual.Inversions;

Tshyear.EP_CAMS(3:10) = ep_nep_95CI_shyear.Inversions;

Tannual.CP_CAMS(3:10) = cp_nep_95CI_annual.Inversions;

Tshyear.CP_CAMS(3:10) = cp_nep_95CI_shyear.Inversions;
clearvars -except Tannual Tshyear

writetable(Tannual, './output/nep-responses-annual.xlsx', 'WriteRowNames',1);
writetable(Tshyear, './output/nep-responses-shyear.xlsx', 'WriteRowNames',1);

