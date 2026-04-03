%% Define mask parameters
voxel = RUC_dim.*1e-6;
Resolution = voxel./RUC_dim;
sheetlet_count = length(Modified_sheetlet_regions);


%% Generate mask

tic
%innitiate mask
NbPixel=round(voxel./Resolution);
cell_mask = zeros(NbPixel);
cell_mask_id = zeros(NbPixel);

%list voxels in each sheetlet
Sheetletwise_voxel_list = Mask_generate_Toolbox.Gen_sheetlet_voxel_list(Modified_sheetlet_regions,NbPixel);
Sheetletwise_mid_voxel_list = Mask_generate_Toolbox.Gen_sheetlet_voxel_list(Modified_mid_sheetlets,NbPixel);
%%
%list of cell_ids_that each voxel belongs to (A sheetletwise_list)
sheetlet_wise_voxel_cell_id = cell(sheetlet_count,1);

%Create mask inside sheetlets(To this point it does not include mask for
%inter sheetlet branches)

cell_no_add = 0;
for stlt = 1:sheetlet_count
    
    objects = shtlt_wise_cell_list{stlt,1};
    voxels = Sheetletwise_voxel_list{stlt,1};
    points = voxels.*1e-6; %to convert voxel to voxel coordinate
    
    [In_out_list,cell_id_list] = Mask_generate_Toolbox.Collison_detection(objects,points);
    
    ICV_voxel_num = find(In_out_list);
    ICV_voxels = voxels(ICV_voxel_num,:);
    
    %convert the cell id that voxel belongs to into global cell id
    cell_id_list_global = cell_id_list+cell_no_add;
    cell_no_add = cell_no_add + size(objects,1);
    
    
    for vox=1:length(ICV_voxels)
        cell_mask(ICV_voxels(vox,1),ICV_voxels(vox,2),ICV_voxels(vox,3)) = 1;
        vox_id = ICV_voxel_num(vox,1);
        cell_mask_id(ICV_voxels(vox,1),ICV_voxels(vox,2),ICV_voxels(vox,3)) = cell_id_list_global(vox_id,1);
    end
    
    %Generate mask with every voxel has the value of id of which cell it
    %belongs to
    
    sheetlet_wise_voxel_cell_id{stlt,1} = cell_id_list;
    
    

end

%save('mask_without_branches.mat')

%Filtering voxels
%(Identify voxels which are not labeled yet and in the cleavage planes)

%Mark voxels inside middle of the sheetlets (where inter sheetlet branches
%are not present) to remove those voxels from mask generation step for
%inter sheetlet branches

cell_mask_identifier = zeros(NbPixel);
rgn_count = size(Sheetletwise_mid_voxel_list,1);
for stlt = 1:rgn_count
    mid_voxels = Sheetletwise_mid_voxel_list{stlt,1};

    for vox=1:length(mid_voxels)
        cell_mask_identifier(mid_voxels(vox,1),mid_voxels(vox,2),mid_voxels(vox,3)) = 1;
    end
end

%Mark already marked voxels in the mask identifier
marked_voxels = find(cell_mask);
cell_mask_identifier(marked_voxels) = 1;

%List of voxels which did not marked yet and in the cleavage planes
filtered_voxel_id = find(~cell_mask_identifier);
[XX1,YY1,ZZ1] = ndgrid(1:1:NbPixel(1),1:1:NbPixel(2),1:1:NbPixel(3));
all_voxel_list = [XX1(:) YY1(:) ZZ1(:)];

filtered_voxel_list = all_voxel_list(filtered_voxel_id,:); % list of voxels in the cleavage planes

%Create mask for inter sheetlet branches
objects = inter_shtlt_branch_list;
voxels = filtered_voxel_list;
points = voxels.*1e-6; %to convert voxel to voxel coordinate

[In_out_list,cell_id_list] = Mask_generate_Toolbox.Collison_detection(objects,points);
ICV_voxel_num = find(In_out_list);
ICV_voxels = voxels(ICV_voxel_num,:);

cell_id_list_global = cell_id_list+cell_no_add;

for vox=1:length(ICV_voxels)
    cell_mask(ICV_voxels(vox,1),ICV_voxels(vox,2),ICV_voxels(vox,3)) = 1;
    vox_id = ICV_voxel_num(vox,1);
    cell_mask_id(ICV_voxels(vox,1),ICV_voxels(vox,2),ICV_voxels(vox,3)) = cell_id_list_global(vox_id,1);
end

sheetlet_wise_voxel_cell_id{end+1,1} = cell_id_list;
toc

%% Get a voxel from mask to perform diffusion simmulation

domain_min_point = [60,100, 100];
domain_dim = [500,500,500];
mask_diff_domain = Diff_domain_mask(cell_mask,domain_min_point,domain_dim);
vol_frac = length(find(mask_diff_domain))/(500^3)

%% functions

function Diff_mask = Diff_domain_mask(cell_mask,domain_min_point,domain_dim)


Diff_mask = zeros(domain_dim(1),domain_dim(2),domain_dim(3));


for x_in = 1:domain_dim(1)
    for y_in = 1:domain_dim(2)
        for z_in = 1:domain_dim(3)
            Diff_mask(x_in,y_in,z_in) = cell_mask(domain_min_point(1)+x_in-1,domain_min_point(2)+y_in-1,domain_min_point(3)+z_in-1);
        end
    end
end
end