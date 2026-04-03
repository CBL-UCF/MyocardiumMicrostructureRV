clc; clear;close all;

%% Generate microstructural connectivity table 

global RUC_w RUC_t RUC_l Mesoscale ModelName

ModelName = 'Micro_test_ORI';

% Define whether generating microscale RUC or Mesoscale RUC
Mesoscale = 1;

% Define the size of RUC 
RUC_w = 200;  
RUC_t = 80;
RUC_l = 500;

a_gen_contable

%% Generate mask for diff simulation and volume fraction calculation

% clear;
% clc;

ModelName = 'Meso_paper_figs';

%%Load required parameters from saved workspace
% workspace_address = append(ModelName,'_conData.mat');
% 
% load(workspace_address,'Final_branches');
% load(workspace_address,'non_shrink_vor_ruc');
% load(workspace_address,'RUC_con_table_branched');
% load(workspace_address,'RUC_con_table_sheetlet_wise');
% load(workspace_address,'RUC_con_v1_branched');
% load(workspace_address,'RUC_con_v2_branched');
% load(workspace_address,'RUC_v_coord');
% load(workspace_address,'Shrinked_after_inf');
% load(workspace_address,'Shrinked_after_inf_mid');
% load(workspace_address,'RUC_w','RUC_t','RUC_l');

% translate cells to origin and generate required datasets
d_preprocess_condata

save(append(ModelName,'_mechInput.mat'),"RUC_w","RUC_l","RUC_t","Modified_coordinates","Modified_sheetlet_regions","RUC_con_table_sheetlet_wise")

%%
% Generate mask
%e_generate_mask