% Main file for carbon cycle impacts of CP vs. EP ENSO events

%% Calculate CP and EP indices from COBE-2 SST, and map/plot indices
cobe2_sst_cp_ep; % Figure 1

%% Define regions based on MODIS land cover and CRU climate
make_mcd12c1_koppen_regions; % Fig. S1

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

%% Estimate responses of GPP to CP/EP ENSO events 
% Global
cp_ep_gpp_lue; % slow!
clear all;
cp_ep_gpp_mstmip; % slow!
clear all;

% Regional
cp_ep_gpp_lue_regional;
clear all;
cp_ep_gpp_mstmip_regional;
clear all;

%% Estimate responses of NEP to CP/EP ENSO events 
% Global
cp_ep_nep_inversions; 
clear all;
cp_ep_nep_mstmip; % slow!
clear all;

% Regional
cp_ep_nep_inversions_regional;
clear all;
cp_ep_nep_mstmip_regional;
clear all;

%% Main figures
make_cp_ep_wavelet_coherence; % Figure 2
make_cp_ep_gpp_maps; % Figure 3
make_cp_ep_ahlstrom_boxplots % Figure 4

%% Build supplementary tables
build_gpp_response_tables; % Tables S#-S#
build_nep_response_tables; % Tables S#-S#
build_gpp_model_agreement_table_figures; % Table S#
build_nep_model_agreement_table_figures; % Table S#

%% Make supplementary figures
supplemental_make_cp_ep_wavelet_coherence_with_phase_arrows; % Figure #
supplemental_make_cp_ep_regional_boxplots % Figure #
supplemental_make_cp_ep_nep_maps; % Figure #
supplemental_cp_ep_cru; % Figs. #-#
supplemental_mauna_loa; % Fig. #-#
supplemental_make_cp_ep_gpp_maps_byYear; % Fig. #

