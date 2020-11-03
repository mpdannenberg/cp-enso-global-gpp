% Build GPP and NEP tables

T = table('Size',[15 2],...
    'RowNames',{'Global','Tropics','Tropical and S. Africa','Tropical South America','Australia','Tropical Asia','Eastern U.S.','Europe','The Sahel','W. North America','Tropical Forests','Extratropical forests','Tundra & Arctic Shrubland','Grass/Crops','Semiarid'},...
    'VariableNames',{'Jul_Jun','Jan_Dec'},...
    'VariableTypes',{'single','single'});

% GPP global
load ./data/cp_ep_gpp_lue.mat;
global_diff_annual = abs(CP_GPP_global_annual_beta) - abs(EP_GPP_global_annual_beta);
global_diff_shyear = abs(CP_GPP_global_shyear_beta) - abs(EP_GPP_global_shyear_beta);
tropics_diff_annual = abs(CP_GPP_tropics_annual_beta) - abs(EP_GPP_tropics_annual_beta);
tropics_diff_shyear = abs(CP_GPP_tropics_shyear_beta) - abs(EP_GPP_tropics_shyear_beta);
clearvars -except T global_* tropics_*;

load ./data/cp_ep_gpp_mstmip.mat;
global_diff_annual = [global_diff_annual (abs(CP_GPP_global_annual_beta) - abs(EP_GPP_global_annual_beta))];
global_diff_shyear = [global_diff_shyear (abs(CP_GPP_global_shyear_beta) - abs(EP_GPP_global_shyear_beta))];
tropics_diff_annual = [tropics_diff_annual (abs(CP_GPP_tropics_annual_beta) - abs(EP_GPP_tropics_annual_beta))];
tropics_diff_shyear = [tropics_diff_shyear (abs(CP_GPP_tropics_shyear_beta) - abs(EP_GPP_tropics_shyear_beta))];
clearvars -except T global_* tropics_*;

T.Jan_Dec('Global') = sum(global_diff_annual > 0);
T.Jul_Jun('Global') = sum(global_diff_shyear > 0);
T.Jan_Dec('Tropics') = sum(tropics_diff_annual > 0);
T.Jul_Jun('Tropics') = sum(tropics_diff_shyear > 0);

clearvars -except T;

% GPP regional
load ./data/cp_ep_gpp_lue_regional.mat;
africa_diff_annual = abs(CP_GPP_africa_annual_beta) - abs(EP_GPP_africa_annual_beta);
africa_diff_shyear = abs(CP_GPP_africa_shyear_beta) - abs(EP_GPP_africa_shyear_beta);
amazon_diff_annual = abs(CP_GPP_amazon_annual_beta) - abs(EP_GPP_amazon_annual_beta);
amazon_diff_shyear = abs(CP_GPP_amazon_shyear_beta) - abs(EP_GPP_amazon_shyear_beta);
austr_diff_annual = abs(CP_GPP_austr_annual_beta) - abs(EP_GPP_austr_annual_beta);
austr_diff_shyear = abs(CP_GPP_austr_shyear_beta) - abs(EP_GPP_austr_shyear_beta);
casia_diff_annual = abs(CP_GPP_casia_annual_beta) - abs(EP_GPP_casia_annual_beta);
casia_diff_shyear = abs(CP_GPP_casia_shyear_beta) - abs(EP_GPP_casia_shyear_beta);
eastus_diff_annual = abs(CP_GPP_eastus_annual_beta) - abs(EP_GPP_eastus_annual_beta);
eastus_diff_shyear = abs(CP_GPP_eastus_shyear_beta) - abs(EP_GPP_eastus_shyear_beta);
europe_diff_annual = abs(CP_GPP_europe_annual_beta) - abs(EP_GPP_europe_annual_beta);
europe_diff_shyear = abs(CP_GPP_europe_shyear_beta) - abs(EP_GPP_europe_shyear_beta);
extratropical_diff_annual = abs(CP_GPP_extratropical_annual_beta) - abs(EP_GPP_extratropical_annual_beta);
extratropical_diff_shyear = abs(CP_GPP_extratropical_shyear_beta) - abs(EP_GPP_extratropical_shyear_beta);
grass_diff_annual = abs(CP_GPP_grass_annual_beta) - abs(EP_GPP_grass_annual_beta);
grass_diff_shyear = abs(CP_GPP_grass_shyear_beta) - abs(EP_GPP_grass_shyear_beta);
sahel_diff_annual = abs(CP_GPP_sahel_annual_beta) - abs(EP_GPP_sahel_annual_beta);
sahel_diff_shyear = abs(CP_GPP_sahel_shyear_beta) - abs(EP_GPP_sahel_shyear_beta);
semiarid_diff_annual = abs(CP_GPP_semiarid_annual_beta) - abs(EP_GPP_semiarid_annual_beta);
semiarid_diff_shyear = abs(CP_GPP_semiarid_shyear_beta) - abs(EP_GPP_semiarid_shyear_beta);
tropical_diff_annual = abs(CP_GPP_tropical_annual_beta) - abs(EP_GPP_tropical_annual_beta);
tropical_diff_shyear = abs(CP_GPP_tropical_shyear_beta) - abs(EP_GPP_tropical_shyear_beta);
tundra_diff_annual = abs(CP_GPP_tundra_annual_beta) - abs(EP_GPP_tundra_annual_beta);
tundra_diff_shyear = abs(CP_GPP_tundra_shyear_beta) - abs(EP_GPP_tundra_shyear_beta);
westna_diff_annual = abs(CP_GPP_westna_annual_beta) - abs(EP_GPP_westna_annual_beta);
westna_diff_shyear = abs(CP_GPP_westna_shyear_beta) - abs(EP_GPP_westna_shyear_beta);
clear CP_* EP_* GPP_* cp* ep* syear eyear nm nt years yr models;

load ./data/cp_ep_gpp_mstmip_regional.mat;
africa_diff_annual = [africa_diff_annual (abs(CP_GPP_africa_annual_beta) - abs(EP_GPP_africa_annual_beta))];
africa_diff_shyear = [africa_diff_shyear (abs(CP_GPP_africa_shyear_beta) - abs(EP_GPP_africa_shyear_beta))];
amazon_diff_annual = [amazon_diff_annual (abs(CP_GPP_amazon_annual_beta) - abs(EP_GPP_amazon_annual_beta))];
amazon_diff_shyear = [amazon_diff_shyear (abs(CP_GPP_amazon_shyear_beta) - abs(EP_GPP_amazon_shyear_beta))];
austr_diff_annual = [austr_diff_annual (abs(CP_GPP_austr_annual_beta) - abs(EP_GPP_austr_annual_beta))];
austr_diff_shyear = [austr_diff_shyear (abs(CP_GPP_austr_shyear_beta) - abs(EP_GPP_austr_shyear_beta))];
casia_diff_annual = [casia_diff_annual (abs(CP_GPP_casia_annual_beta) - abs(EP_GPP_casia_annual_beta))];
casia_diff_shyear = [casia_diff_shyear (abs(CP_GPP_casia_shyear_beta) - abs(EP_GPP_casia_shyear_beta))];
eastus_diff_annual = [eastus_diff_annual (abs(CP_GPP_eastus_annual_beta) - abs(EP_GPP_eastus_annual_beta))];
eastus_diff_shyear = [eastus_diff_shyear (abs(CP_GPP_eastus_shyear_beta) - abs(EP_GPP_eastus_shyear_beta))];
europe_diff_annual = [europe_diff_annual (abs(CP_GPP_europe_annual_beta) - abs(EP_GPP_europe_annual_beta))];
europe_diff_shyear = [europe_diff_shyear (abs(CP_GPP_europe_shyear_beta) - abs(EP_GPP_europe_shyear_beta))];
extratropical_diff_annual = [extratropical_diff_annual (abs(CP_GPP_extratropical_annual_beta) - abs(EP_GPP_extratropical_annual_beta))];
extratropical_diff_shyear = [extratropical_diff_shyear (abs(CP_GPP_extratropical_shyear_beta) - abs(EP_GPP_extratropical_shyear_beta))];
grass_diff_annual = [grass_diff_annual (abs(CP_GPP_grass_annual_beta) - abs(EP_GPP_grass_annual_beta))];
grass_diff_shyear = [grass_diff_shyear (abs(CP_GPP_grass_shyear_beta) - abs(EP_GPP_grass_shyear_beta))];
sahel_diff_annual = [sahel_diff_annual (abs(CP_GPP_sahel_annual_beta) - abs(EP_GPP_sahel_annual_beta))];
sahel_diff_shyear = [sahel_diff_shyear (abs(CP_GPP_sahel_shyear_beta) - abs(EP_GPP_sahel_shyear_beta))];
semiarid_diff_annual = [semiarid_diff_annual (abs(CP_GPP_semiarid_annual_beta) - abs(EP_GPP_semiarid_annual_beta))];
semiarid_diff_shyear = [semiarid_diff_shyear (abs(CP_GPP_semiarid_shyear_beta) - abs(EP_GPP_semiarid_shyear_beta))];
tropical_diff_annual = [tropical_diff_annual (abs(CP_GPP_tropical_annual_beta) - abs(EP_GPP_tropical_annual_beta))];
tropical_diff_shyear = [tropical_diff_shyear (abs(CP_GPP_tropical_shyear_beta) - abs(EP_GPP_tropical_shyear_beta))];
tundra_diff_annual = [tundra_diff_annual (abs(CP_GPP_tundra_annual_beta) - abs(EP_GPP_tundra_annual_beta))];
tundra_diff_shyear = [tundra_diff_shyear (abs(CP_GPP_tundra_shyear_beta) - abs(EP_GPP_tundra_shyear_beta))];
westna_diff_annual = [westna_diff_annual (abs(CP_GPP_westna_annual_beta) - abs(EP_GPP_westna_annual_beta))];
westna_diff_shyear = [westna_diff_shyear (abs(CP_GPP_westna_shyear_beta) - abs(EP_GPP_westna_shyear_beta))];
clear CP_* EP_* GPP_* cp* ep* syear eyear nm nt years yr models nx ny idx;

T.Jan_Dec('Tropical and S. Africa') = sum(africa_diff_annual > 0);
T.Jul_Jun('Tropical and S. Africa') = sum(africa_diff_shyear > 0);
T.Jan_Dec('Tropical South America') = sum(amazon_diff_annual > 0);
T.Jul_Jun('Tropical South America') = sum(amazon_diff_shyear > 0);
T.Jan_Dec('Australia') = sum(austr_diff_annual > 0);
T.Jul_Jun('Australia') = sum(austr_diff_shyear > 0);
T.Jan_Dec('Tropical Asia') = sum(casia_diff_annual > 0);
T.Jul_Jun('Tropical Asia') = sum(casia_diff_shyear > 0);
T.Jan_Dec('Eastern U.S.') = sum(eastus_diff_annual > 0);
T.Jul_Jun('Eastern U.S.') = sum(eastus_diff_shyear > 0);
T.Jan_Dec('Europe') = sum(europe_diff_annual > 0);
T.Jul_Jun('Europe') = sum(europe_diff_shyear > 0);
T.Jan_Dec('Extratropical forests') = sum(extratropical_diff_annual > 0);
T.Jul_Jun('Extratropical forests') = sum(extratropical_diff_shyear > 0);
T.Jan_Dec('Grass/Crops') = sum(grass_diff_annual > 0);
T.Jul_Jun('Grass/Crops') = sum(grass_diff_shyear > 0);
T.Jan_Dec('The Sahel') = sum(sahel_diff_annual > 0);
T.Jul_Jun('The Sahel') = sum(sahel_diff_shyear > 0);
T.Jan_Dec('Semiarid') = sum(semiarid_diff_annual > 0);
T.Jul_Jun('Semiarid') = sum(semiarid_diff_shyear > 0);
T.Jan_Dec('Tropical Forests') = sum(tropical_diff_annual > 0);
T.Jul_Jun('Tropical Forests') = sum(tropical_diff_shyear > 0);
T.Jan_Dec('Tundra & Arctic Shrubland') = sum(tundra_diff_annual > 0);
T.Jul_Jun('Tundra & Arctic Shrubland') = sum(tundra_diff_shyear > 0);
T.Jan_Dec('W. North America') = sum(westna_diff_annual > 0);
T.Jul_Jun('W. North America') = sum(westna_diff_shyear > 0);

clearvars -except T;
writetable(T, './output/gpp-model-agreement.xlsx');

