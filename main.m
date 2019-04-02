% Main file for carbon cycle impacts of CP vs. EP ENSO events

%% Calculate CP and EP indices from COBE-2 SST, and map/plot indices
cobe2_sst_cp_ep; % Figure 1

%% Read and process modeled GPP data
% Data-driven models
read_gpp_mod17;
read_gpp_ccw;

% Process-based models
read_gpp_mstmip; % Slow!

%% Read and process NEP data
% Data-driven models
read_nep_inversions;

% Process-based models
read_nep_mstmip; % Slow!

%% Estimate responses of GPP to CP/EP ENSO events - Table 1
% Global
cp_ep_gpp_lue; % slow!
cp_ep_gpp_mstmip; % slow!

% Regional
cp_ep_gpp_lue_regional;
cp_ep_gpp_mstmip_regional;

%% Estimate responses of NEP to CP/EP ENSO events - Table 2
% Global
cp_ep_nep_inversions; 
cp_ep_nep_mstmip; % slow!

% Regional
cp_ep_nep_inversions_regional;
cp_ep_nep_mstmip_regional;

%% Main figures
make_cp_ep_gpp_maps; % Figure 2
make_cp_ep_nep_maps; % Figure 3

%% Supplementary figures
make_cp_ep_gpp_maps_byYear; % Figure S1
supplemental_map_regions; % Figure S2
make_cp_ep_gpp_regional_plots_tropics; % Figure S3
make_cp_ep_gpp_regional_plots_temperate; % Figure S4
make_cp_ep_nep_regional_plots_tropics; % Figure S5
make_cp_ep_nep_regional_plots_temperate; % Figure S6
supplemental_mauna_loa; % Fig. S#

