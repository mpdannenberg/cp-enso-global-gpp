% Build GPP and NEP tables

Tannual = table('Size',[15 6],...
    'RowNames',{'Global','Tropics','Tropical and S. Africa','Tropical South America','Australia','Tropical Asia','Eastern U.S.','Europe','The Sahel','W. North America','Tropical Forests','Extratropical forests','Tundra & Arctic Shrubland','Grass/Crops','Semiarid'},...
    'VariableNames',{'EP_MsTMIP','EP_MOD17','EP_CCW','CP_MsTMIP','CP_MOD17','CP_CCW'},...
    'VariableTypes',{'string','string','string','string','string','string'});
Tshyear = table('Size',[15 6],...
    'RowNames',{'Global','Tropics','Tropical and S. Africa','Tropical South America','Australia','Tropical Asia','Eastern U.S.','Europe','The Sahel','W. North America','Tropical Forests','Extratropical forests','Tundra & Arctic Shrubland','Grass/Crops','Semiarid'},...
    'VariableNames',{'EP_MsTMIP','EP_MOD17','EP_CCW','CP_MsTMIP','CP_MOD17','CP_CCW'},...
    'VariableTypes',{'string','string','string','string','string','string'});

% GPP global LUE
load ./data/cp_ep_gpp_lue.mat;
Tannual.EP_MOD17(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_GPP_global_annual_beta(2)/1000, EP_GPP_global_annual_beta_CI(2)/1000));
Tannual.EP_CCW(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_GPP_global_annual_beta(1)/1000, EP_GPP_global_annual_beta_CI(1)/1000));
Tannual.EP_MOD17(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_GPP_tropics_annual_beta(2)/1000, EP_GPP_tropics_annual_beta_CI(2)/1000));
Tannual.EP_CCW(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_GPP_tropics_annual_beta(1)/1000, EP_GPP_tropics_annual_beta_CI(1)/1000));

Tshyear.EP_MOD17(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_GPP_global_shyear_beta(2)/1000, EP_GPP_global_shyear_beta_CI(2)/1000));
Tshyear.EP_CCW(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_GPP_global_shyear_beta(1)/1000, EP_GPP_global_shyear_beta_CI(1)/1000));
Tshyear.EP_MOD17(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_GPP_tropics_shyear_beta(2)/1000, EP_GPP_tropics_shyear_beta_CI(2)/1000));
Tshyear.EP_CCW(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_GPP_tropics_shyear_beta(1)/1000, EP_GPP_tropics_shyear_beta_CI(1)/1000));

Tannual.CP_MOD17(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_GPP_global_annual_beta(2)/1000, CP_GPP_global_annual_beta_CI(2)/1000));
Tannual.CP_CCW(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_GPP_global_annual_beta(1)/1000, CP_GPP_global_annual_beta_CI(1)/1000));
Tannual.CP_MOD17(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_GPP_tropics_annual_beta(2)/1000, CP_GPP_tropics_annual_beta_CI(2)/1000));
Tannual.CP_CCW(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_GPP_tropics_annual_beta(1)/1000, CP_GPP_tropics_annual_beta_CI(1)/1000));

Tshyear.CP_MOD17(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_GPP_global_shyear_beta(2)/1000, CP_GPP_global_shyear_beta_CI(2)/1000));
Tshyear.CP_CCW(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_GPP_global_shyear_beta(1)/1000, CP_GPP_global_shyear_beta_CI(1)/1000));
Tshyear.CP_MOD17(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_GPP_tropics_shyear_beta(2)/1000, CP_GPP_tropics_shyear_beta_CI(2)/1000));
Tshyear.CP_CCW(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_GPP_tropics_shyear_beta(1)/1000, CP_GPP_tropics_shyear_beta_CI(1)/1000));
clearvars -except Tannual Tshyear

% GPP regional LUE
load ./data/cp_ep_gpp_lue_regional.mat;
Tannual.EP_MOD17(3:end) = ep_gpp_95CI_annual.MOD17;
Tannual.EP_CCW(3:end) = ep_gpp_95CI_annual.CCW;

Tshyear.EP_MOD17(3:end) = ep_gpp_95CI_shyear.MOD17;
Tshyear.EP_CCW(3:end) = ep_gpp_95CI_shyear.CCW;

Tannual.CP_MOD17(3:end) = cp_gpp_95CI_annual.MOD17;
Tannual.CP_CCW(3:end) = cp_gpp_95CI_annual.CCW;

Tshyear.CP_MOD17(3:end) = cp_gpp_95CI_shyear.MOD17;
Tshyear.CP_CCW(3:end) = cp_gpp_95CI_shyear.CCW;
clearvars -except Tannual Tshyear

% GPP global MsTMIP
load ./data/cp_ep_gpp_mstmip.mat;
Tannual.EP_MsTMIP(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_GPP_global_annual_mean_beta/1000, EP_GPP_global_annual_mean_beta_CI/1000));
Tannual.EP_MsTMIP(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_GPP_tropics_annual_mean_beta/1000, EP_GPP_tropics_annual_mean_beta_CI/1000));

Tshyear.EP_MsTMIP(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_GPP_global_shyear_mean_beta/1000, EP_GPP_global_shyear_mean_beta_CI/1000));
Tshyear.EP_MsTMIP(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                EP_GPP_tropics_shyear_mean_beta/1000, EP_GPP_tropics_shyear_mean_beta_CI/1000));

Tannual.CP_MsTMIP(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_GPP_global_annual_mean_beta/1000, CP_GPP_global_annual_mean_beta_CI/1000));
Tannual.CP_MsTMIP(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_GPP_tropics_annual_mean_beta/1000, CP_GPP_tropics_annual_mean_beta_CI/1000));

Tshyear.CP_MsTMIP(1) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_GPP_global_shyear_mean_beta/1000, CP_GPP_global_shyear_mean_beta_CI/1000));
Tshyear.CP_MsTMIP(2) = cellstr(sprintf(['%.02f ',char(177),' %.02f'],...
                CP_GPP_tropics_shyear_mean_beta/1000, CP_GPP_tropics_shyear_mean_beta_CI/1000));
clearvars -except Tannual Tshyear

% GPP regional LUE
load ./data/cp_ep_gpp_mstmip_regional.mat;
Tannual.EP_MsTMIP(3:end) = ep_gpp_95CI_annual.MsTMIP;

Tshyear.EP_MsTMIP(3:end) = ep_gpp_95CI_shyear.MsTMIP;

Tannual.CP_MsTMIP(3:end) = cp_gpp_95CI_annual.MsTMIP;

Tshyear.CP_MsTMIP(3:end) = cp_gpp_95CI_shyear.MsTMIP;
clearvars -except Tannual Tshyear

writetable(Tannual, './output/gpp-responses-annual.xlsx');
writetable(Tshyear, './output/gpp-responses-shyear.xlsx');

