%% Define mask parameters
% Convert RUC dimensions to physical voxel size (meters)
voxel = RUC_dim.*1e-6;

% Define voxel resolution (size of each voxel in meters)
Resolution = voxel./RUC_dim;

% Number of sheetlets in the model
sheetlet_count = length(Modified_sheetlet_regions);


%% Generate mask

% This section builds a 3D voxelized mask representing intracellular regions

tic

% Initialize mask arrays
% cell_mask: binary mask (1 = intracellular, 0 = extracellular)
% cell_mask_id: stores the ID of the cell occupying each voxel

NbPixel=round(voxel./Resolution);
cell_mask = zeros(NbPixel);
cell_mask_id = zeros(NbPixel);

% Generate voxel lists for each sheetlet
Sheetletwise_voxel_list = Mask_generate_Toolbox.Gen_sheetlet_voxel_list(Modified_sheetlet_regions,NbPixel);
Sheetletwise_mid_voxel_list = Mask_generate_Toolbox.Gen_sheetlet_voxel_list(Modified_mid_sheetlets,NbPixel);


% Initialize storage for voxel-to-cell mapping (sheetlet-wise)
sheetlet_wise_voxel_cell_id = cell(sheetlet_count,1);

% Generate mask inside sheetlets (excluding inter-sheetlet branches)

cell_no_add = 0;
for stlt = 1:sheetlet_count

    % Get cell objects and voxel coordinates for current sheetlet    
    objects = shtlt_wise_cell_list{stlt,1};
    voxels = Sheetletwise_voxel_list{stlt,1};
    points = voxels.*1e-6; %to convert voxel to voxel coordinate
    
    % Perform collision detection (voxel inside any cylindrical cell)
    [In_out_list,cell_id_list] = Mask_generate_Toolbox.Collison_detection(objects,points);
    
    % Extract voxels inside intracellular volume (ICV)
    ICV_voxel_num = find(In_out_list);
    ICV_voxels = voxels(ICV_voxel_num,:);
    
    % Convert local cell IDs to global IDs
    cell_id_list_global = cell_id_list+cell_no_add;
    cell_no_add = cell_no_add + size(objects,1);
    
    % Assign voxel values in mask
    for vox=1:length(ICV_voxels)
        cell_mask(ICV_voxels(vox,1),ICV_voxels(vox,2),ICV_voxels(vox,3)) = 1;
        vox_id = ICV_voxel_num(vox,1);
        cell_mask_id(ICV_voxels(vox,1),ICV_voxels(vox,2),ICV_voxels(vox,3)) = cell_id_list_global(vox_id,1);
    end
    
    % Store voxel-to-cell mapping for current sheetlet    
    sheetlet_wise_voxel_cell_id{stlt,1} = cell_id_list;
    
    

end

%save('mask_without_branches.mat')

%% Filtering voxels
% Identify voxels not yet assigned (i.e., cleavage planes)

% Initialize mask identifier (tracks processed voxels)
cell_mask_identifier = zeros(NbPixel);

% Mark voxels inside mid-sheetlet regions (exclude from branching step)
rgn_count = size(Sheetletwise_mid_voxel_list,1);
for stlt = 1:rgn_count
    mid_voxels = Sheetletwise_mid_voxel_list{stlt,1};

    for vox=1:length(mid_voxels)
        cell_mask_identifier(mid_voxels(vox,1),mid_voxels(vox,2),mid_voxels(vox,3)) = 1;
    end
end

% Mark already assigned intracellular voxels
marked_voxels = find(cell_mask);
cell_mask_identifier(marked_voxels) = 1;

% Extract voxels not yet assigned (candidate cleavage plane voxels)
filtered_voxel_id = find(~cell_mask_identifier);

% Generate full voxel grid
[XX1,YY1,ZZ1] = ndgrid(1:1:NbPixel(1),1:1:NbPixel(2),1:1:NbPixel(3));
all_voxel_list = [XX1(:) YY1(:) ZZ1(:)];

% Extract only unassigned voxels
filtered_voxel_list = all_voxel_list(filtered_voxel_id,:); % list of voxels in the cleavage planes

%% Create mask for inter-sheetlet branches
objects = inter_shtlt_branch_list;
voxels = filtered_voxel_list;
points = voxels.*1e-6; %to convert voxel to voxel coordinate

% Collision detection for branches
[In_out_list,cell_id_list] = Mask_generate_Toolbox.Collison_detection(objects,points);
ICV_voxel_num = find(In_out_list);
ICV_voxels = voxels(ICV_voxel_num,:);

% Convert to global cell IDs
cell_id_list_global = cell_id_list+cell_no_add;

% Assign branch voxels into mask
for vox=1:length(ICV_voxels)
    cell_mask(ICV_voxels(vox,1),ICV_voxels(vox,2),ICV_voxels(vox,3)) = 1;
    vox_id = ICV_voxel_num(vox,1);
    cell_mask_id(ICV_voxels(vox,1),ICV_voxels(vox,2),ICV_voxels(vox,3)) = cell_id_list_global(vox_id,1);
end

% Store branch voxel mapping
sheetlet_wise_voxel_cell_id{end+1,1} = cell_id_list;
toc

%% Extract subdomain for diffusion simulation

% Define subdomain origin and dimensions

domain_min_point = [60,100, 100];
domain_dim = [500,500,500];

% Extract diffusion simulation mask
mask_diff_domain = Diff_domain_mask(cell_mask,domain_min_point,domain_dim);

% Compute volume fraction of intracellular region
vol_frac = length(find(mask_diff_domain))/(500^3)

%% functions

function Diff_mask = Diff_domain_mask(cell_mask,domain_min_point,domain_dim)
    % Diff_domain_mask
    % Extracts a subdomain from the full mask for diffusion simulations.
    %
    % Inputs:
    %   cell_mask        - Full 3D mask
    %   domain_min_point - Starting index of subdomain
    %   domain_dim       - Dimensions of subdomain
    %
    % Output:
    %   Diff_mask        - Extracted subdomain mask
    
    Diff_mask = zeros(domain_dim(1),domain_dim(2),domain_dim(3));
    
    
    for x_in = 1:domain_dim(1)
        for y_in = 1:domain_dim(2)
            for z_in = 1:domain_dim(3)
                Diff_mask(x_in,y_in,z_in) = cell_mask(domain_min_point(1)+x_in-1,domain_min_point(2)+y_in-1,domain_min_point(3)+z_in-1);
            end
        end
    end
end