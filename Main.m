clc; clear; close all;

%% Generate microstructural connectivity table 

% global RUC_w RUC_t RUC_l shtlt_w shtlt_t Mesoscale ModelName shtlt_shrink_ratios

% This script generates or loads a representative unit cell (RUC)
% microstructure, applies coordinate transformations, prepares cell-level
% geometry information, and initializes mask generation.


ModelName = 'Meso_celldimensions';

% Flag to control whether connectivity is generated or loaded from file
Construct_connectivity = 1 ;

if Construct_connectivity

    % Define whether generating microscale RUC or Mesoscale RUC
    Mesoscale = 1;
    
    % Define the size of Representative Unit Cell (RUC)
    RUC_w = 700;  
    RUC_t = 700;
    RUC_l = 600;
    
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
    myocyte_length_range = [100, 130]; % [minimum myocyte length, maximum myocyte length]
    max_fibredir_deviation = innitial_myocyte_diameter / 6; % maximum deviation in sheetlet-normal plane per cardiomyocyte length
    
    
    % Sheetlet innitial volume fraction
    shtlt_ICV = 0.95;
    appx_target_ICV = 1.36; % This value can be >1 since at this point cardiomyocytes considered as cyllinders and are overlapping at nodes
    
    %Branching
    branching_reg_thickness = 40; % Region to generate inter sheetlet branches from common boundary between sheetets
    
    % Technical values
    randomization_breakpoint = 500000;
    
    % Generate RUC cardiomyocyte connectivity, sheetlet regions, and branches
    [RUC_con_table_branched, RUC_v_coord, RUC_con_v1_branched, RUC_con_v2_branched, Final_branches, non_shrink_vor_ruc, RUC_con_table_sheetlet_wise, Shrinked_after_inf, Shrinked_after_inf_mid] = RUC_microstructure_Toolbox.generate_ruc_myoconnectivity(RUC_w, ...
        RUC_t, RUC_l, shtlt_w, shtlt_t, ...
        shtlt_shrink_ratios, innitial_myocyte_diameter, myocyte_length_range, ...
        shtlt_ICV, appx_target_ICV, max_fibredir_deviation, branching_reg_thickness, ...
        randomization_breakpoint);
    RUC_shtletwise_cell_id_caps = RUC_microstructure_Toolbox.shtltwise_cell_ids(RUC_con_table_sheetlet_wise);

    
    % Save microstructure connectivity information
    save(append(ModelName,'_conData.mat'));

else

    % Load required parameters from previously saved workspace
    workspace_address = append(ModelName,'_conData.mat');
    
    load(workspace_address,'Final_branches', ...
        'non_shrink_vor_ruc', ...
        'RUC_con_table_branched', ...
        'RUC_con_table_sheetlet_wise', ...
        'RUC_con_v1_branched', ...
        'RUC_con_v2_branched', ...
        'RUC_v_coord', ...
        'RUC_shtletwise_cell_id_caps', ...
        'Shrinked_after_inf', ...
        'Shrinked_after_inf_mid', ...
        'RUC_w','RUC_t','RUC_l');
end



%% Perform alterations on coordinate system (Rotations and Translations of microstructure) and Cutout a cubical region of interest

% Translate the full RUC so that the minimum point of the model aligns with
% a defined reference point.

% Define target minimum coordinate after translation
min_point =[0,0,0];

% Find coordinate translation amount and overall RUC dimensions
[coord_details] = Coordinate_Toolbox.Modification_amount_and_dimension(Shrinked_after_inf,RUC_v_coord,min_point);
modify_amount = coord_details(1,:);
RUC_dim = ceil(coord_details(2,:));

% Calculate translated vertex coordinates for myocytes and sheetlets
Modified_point_list = Coordinate_Toolbox.Modify_RUC_coord_mask(RUC_v_coord,modify_amount); %type- cell
Modified_coordinates = Coordinate_Toolbox.cell_to_double_coordinates(Modified_point_list); %type- double
Modified_sheetlet_regions = Coordinate_Toolbox.Modify_shtlt_coord(Shrinked_after_inf,modify_amount);
Modified_mid_sheetlets = Coordinate_Toolbox.Modify_shtlt_coord(Shrinked_after_inf_mid,modify_amount);

% Plot translated cell network
fig3 = Plotting_Toolbox.plot_microstructure(3, Modified_point_list, RUC_con_v1_branched, RUC_con_v2_branched);

% Save translated microstructure and sheetlet information for mechanics input
save(append(ModelName,'_mechInput.mat'),"RUC_w","RUC_l","RUC_t","Modified_coordinates","Modified_sheetlet_regions","RUC_con_table_sheetlet_wise")

% Cutout a ROI - TODO
ROI_min_point = [30, 50, 80];
ROI_dimension = [300, 250, 250];
RUC_con = [RUC_con_v1_branched; RUC_con_v2_branched]';
[ROI_Pts, ROI_Conn, ROI_parentElem, ROI_tRange] = ROI_cutout_toolbox.clipLineElementsToRUC(Modified_coordinates, RUC_con, ROI_min_point, ROI_dimension, 1e-6);
[front, back, left, right, down, up] = ROI_cutout_toolbox.find_boundary_points(ROI_min_point, ROI_min_point + ROI_dimension, ROI_Pts, 0.05);

save(append(ModelName,'_mechInput_ROI.mat'),"ROI_min_point","ROI_dimension","ROI_Pts","ROI_Conn", "ROI_parentElem", "RUC_con_table_sheetlet_wise", 'RUC_shtletwise_cell_id_caps', "front", "back", "left", "right", "down", "up")
%% Generate Voxelwise mesh for vertual MR simulation

% Generate geometric cell lists for the complete network, sheetlet-wise
% myocytes, and inter-sheetlet branches.
all_cell_list = Mask_generation_Toolbox.gen_cell_detail_list(RUC_con_table_branched,Modified_point_list,12);
shtlt_wise_cell_list = Mask_generation_Toolbox.Sheetletwise_cell_details(RUC_con_table_sheetlet_wise,Modified_point_list,12);
inter_shtlt_branch_list = Mask_generation_Toolbox.gen_cell_detail_list(Final_branches,Modified_point_list,12);

% Generate a binary mask (ICV - 1; ECV - 0) and a mask with voxels assigned
% to belonging cell id (ICV - id of the myocyte which a voxel belongs to; ECV - 0)

[mask_ROI, mask_ROI_id, Reshaped_mask, Reshaped_mask_id, vol_frac] = Mask_generation_Toolbox.Generate_mask(RUC_dim, all_cell_list, shtlt_wise_cell_list, inter_shtlt_branch_list, Modified_sheetlet_regions, Modified_mid_sheetlets, ROI_min_point, ROI_dimension, RUC_con_table_sheetlet_wise, Final_branches);

print(vol_frac)

% Save .mat for diffusion simulation
filename = ppend(ModelName,'_diff.mat');
Cells_mask = mask_ROI;
Cell_params = mask_ROI;

save(filename,"Cells_mask", "Cell_params");

% Save .mat of the mask for further use
filename = ppend(ModelName,'_mask_data.mat');
save(filename, "mask_ROI", "mask_ROI_id", "Reshaped_mask", "Reshaped_mask_id", "vol_frac", "all_cell_list", "shtlt_wise_cell_list", "inter_shtlt_branch_list")