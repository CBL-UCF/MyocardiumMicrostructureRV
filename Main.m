clc; clear;close all;

%% Generate microstructural connectivity table 

global RUC_w RUC_t RUC_l shtlt_w shtlt_t Mesoscale ModelName shtlt_shrink_ratios

ModelName = 'Micro_test_ORI';

% Define whether generating microscale RUC or Mesoscale RUC
Mesoscale = 1;

% Define the size of Representative Unit Cell (RUC)
RUC_w = 200;  
RUC_t = 80;
RUC_l = 500;

% Sheetlet width and thickness
shtlt_w = 200;
shtlt_t = 45;

% Shrink ratio for sheetlets in sheetlet and normal directions
% Ratio[1] = shrink ratio for cleavage planes normal to sheetlet direction
% Ratio[2] = shrink ratio for cleavage planes parallel to sheetlet direction
% Each ratio value range  - (0.3, 1)

shtlt_shrink_ratios = [0.955, 0.825]; 

% Cardiomyocyte dimensions
innitial_myocyte_diameter = 10;
myocyte_length_range = [100, 120]; % [minimum myocyte length, maximum myocyte length]
max_fibredir_deviation = innitial_myocyte_diameter / 6; % maximum deviation in sheetlet-normal plane per cardiomyocyte length


% Sheetlet innitial volume fraction
shtlt_ICV = 0.95;


% Technical values
randomization_breakpoint = 500000;

RUC_microstructure_Toolbox.generate_ruc_myoconnectivity(RUC_w, RUC_t, RUC_l, shtlt_w, shtlt_t, ...
    shtlt_shrink_ratios, innitial_myocyte_diameter, myocyte_length_range, shtlt_ICV, ...
    max_fibredir_deviation, randomization_breakpoint)

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