%Mask generation toolbox

classdef Mask_generation_Toolbox

    % Mask_generation_Toolbox
    %
    % Utility class for generating voxelized myocardial microstructure
    % masks from geometric cardiomyocyte and sheetlet representations.
    %
    % Main capabilities:
    %   1. Generate voxel masks for intra-sheetlet cardiomyocytes.
    %   2. Generate voxel masks for inter-sheetlet branches.
    %   3. Assign unique cell identifiers to intracellular voxels.
    %   4. Reshape masks to improve cellular boundaries and continuity.
    %   5. Extract regions of interest (ROIs) for diffusion simulations.
    %   6. Export masks for visualization in VTK format.
    %
    % The resulting masks represent intracellular volume (ICV) and
    % extracellular volume (ECV) at voxel resolution and can be used for
    % diffusion MRI simulations or image-based analysis.
    
    methods(Static)

        function [mask_ROI, mask_ROI_id, Reshaped_mask, Reshaped_mask_id, vol_frac] = Generate_mask(RUC_dim, all_cell_list, shtlt_wise_cell_list, inter_shtlt_branch_list, Modified_sheetlet_regions, Modified_mid_sheetlets, ROI_min_point, ROI_dimention, RUC_con_table_sheetlet_wise, Final_branches, model_name)
            
            % Generate_mask
            %
            % Main workflow for converting the generated myocardial
            % microstructure into a voxelized intracellular mask.
            %
            % Workflow:
            %   Step 1:
            %       Define voxel grid dimensions and sheetlet voxel regions.
            %
            %   Step 2:
            %       Generate intracellular masks for cardiomyocytes located
            %       inside individual sheetlets.
            %
            %   Step 3:
            %       Identify unassigned voxels within cleavage-plane regions.
            %
            %   Step 4:
            %       Generate masks for inter-sheetlet branches.
            %
            %   Step 5:
            %       Reshape the mask to improve cellular continuity and
            %       prevent unrealistic voxel-level overlaps.
            %
            %   Step 6:
            %       Define cell boundaries while preserving connectivity
            %       between physically connected cardiomyocytes.
            %
            %   Step 7:
            %       Extract a user-defined region of interest (ROI) for
            %       diffusion simulations.
            %
            %   Step 8:
            %       Compute intracellular volume fraction and export masks
            %       as VTK files for visualization.
            %
            % Outputs:
            %   mask_ROI          - Binary ROI mask (ICV = 1, ECV = 0)
            %   mask_ROI_id       - ROI mask containing cell IDs
            %   Reshaped_mask     - Full binary mask after reshaping
            %   Reshaped_mask_id  - Full cell-ID mask after reshaping
            %   vol_frac          - Intracellular volume fraction within ROI

            % ____________________________________________________________%
            % Define Mask Parameters

            % Convert RUC dimensions to physical voxel size (meters)
            voxel = RUC_dim.*1e-6;
            % Define voxel resolution (size of each voxel in meters)
            Resolution = voxel./RUC_dim;
            % Number of sheetlets in the model
            sheetlet_count = length(Modified_sheetlet_regions);
            % N voxels in each dimension
            NbPixel=round(voxel./Resolution);

            % ____________________________________________________________%
            % Innitiate mask and identify voxel lists

            Sheetletwise_voxel_list = Mask_generation_Toolbox.Gen_sheetlet_voxel_list(Modified_sheetlet_regions,NbPixel);
            Sheetletwise_mid_voxel_list = Mask_generation_Toolbox.Gen_sheetlet_voxel_list(Modified_mid_sheetlets,NbPixel);

            % ____________________________________________________________%
            % Update mask with cells inside sheetlets (Excluding cleavage planes)

            [cell_mask_inshtlt, cell_mask_id_inshtlt, sheetlet_wise_voxel_cell_id, n_insheetlet_cells] = Mask_generation_Toolbox.update_inside_sheetlets(shtlt_wise_cell_list, sheetlet_count, Sheetletwise_voxel_list, NbPixel);

            % ____________________________________________________________%
            % Update mask with cells in cleavage planes (Inter-sheetlet
            % branches)
            
            % List unasigned voxels in cleavage planes (Here a further
            % shrinked mash is used to remove sheetlet voxels avoiding
            % removing voxels at the boundary which maybe belongs to
            % inter-sheetlet branches which are yet to be updated)
            cleavage_plane_voxels = Mask_generation_Toolbox.filter_unassigned_voxels(NbPixel, Sheetletwise_mid_voxel_list, cell_mask);
               
            [mask_before_reshape, mask_id_before_reshape, intershtlt_branch_id_list] = Mask_generation_Toolbox.update_inter_sheetlet_branches(inter_shtlt_branch_list, cleavage_plane_voxels, n_insheetlet_cells, cell_mask_inshtlt, cell_mask_id_inshtlt);
            
            % Store branch voxel mapping
            sheetlet_wise_voxel_cell_id{end+1,1} = intershtlt_branch_id_list;

            % ____________________________________________________________%
            % Prepare region-wise cell list and voxel list for mask
            % reshaping

            % Add column to each sheeetletwise cell list representing their global id of
            % the cell
            cell_no_add = 0;
            for stlt = 1:sheetlet_count
                for i = 1:size(shtlt_wise_cell_list{stlt,1},1)
                    shtlt_wise_cell_list{stlt,1}(i,9) = i+cell_no_add;
                end
                cell_no_add = cell_no_add + size(shtlt_wise_cell_list{stlt,1},1);
            end
            
            branch_list = inter_shtlt_branch_list(:,1:8).*1e6; %to convert into micro-meters
            for i = 1:size(branch_list,1)
                branch_list(i,9) = i+cell_no_add;
            end

            reshaped_cell_details = all_cell_list(:,1:8).* 1e6; %to convert into micro-meters
            
            reshaped1_rgnwise_cell_list = shtlt_wise_cell_list;
            reshaped1_rgnwise_cell_list{end+1} = branch_list;
            
            Rgn_wise_voxel_list = Sheetletwise_voxel_list;
            Rgn_wise_voxel_list{end+1} = cleavage_plane_voxels;
            % ____________________________________________________________%
            % Reshaping the mask

            % Reshaping the cells to fit inside sheetlets while not overlaping with
            % each other

            % Step 1 - add icv voxels to previous slice to be compatible with shaked next slice
            [shaked_mask_id] = Mask_reshape_Toolbox.reshape_mask_shake_slices(mask_before_reshape, mask_id_before_reshape, reshaped1_rgnwise_cell_list, Rgn_wise_voxel_list);

            
            % Step 2 - Define boundary between cells ----------------------

            % Identify connections between cells in regions to avoid adding
            % ecv between connected cells (e.g - branches)
            [shtlt_wise_connection_list] = Coordinate_Toolbox.find_connected_cells(RUC_con_table_sheetlet_wise,Final_branches);
            [Global_connection_list, N_cells_before_branch] = Coordinate_Toolbox.glob_con_list(shtlt_wise_connection_list);

            % Add boundary between cells
            [Reshaped_mask_id] = Mask_reshape_Toolbox.reshape_mask_define_cellboundary(shaked_mask_id, Global_connection_list, N_cells_before_branch);

            % Plot a reshaped slice for testing purpose
            test = squeeze(Reshaped_mask_id(50,:,:));
            figure
            imagesc(test)

            % ____________________________________________________________%
            
            % Define final masks and volume fractions
            mask_dim_3D = size(Reshaped_mask_id);
            Reshaped_mask = zeros(mask_dim_3D);
            
            icv_ids = find(Reshaped_mask_id);
            Reshaped_mask(icv_ids) = 1;
            
            mask_ROI = Mask_generation_Toolbox.Diff_domain_mask(Reshaped_mask, ROI_min_point, ROI_dimention);
            mask_ROI_id = Mask_generation_Toolbox.Diff_domain_mask(Reshaped_mask_id, ROI_min_point, ROI_dimention);
            
            vol_frac = length(find(mask_ROI)) / (ROI_dimention(1) * ROI_dimention(2) * ROI_dimention(3));
                       
            % ____________________________________________________________%
            % Write mask to a VTK file to visualize
            f_vtk = append(model_name, '_vf' , string(vol_frac), '_ROI.vtk');
            Mat2VTK(f_vtk,mask_ROI,'binary');
            
            
            f_vtk = append(model_name, '_vf' , string(vol_frac), '_fullmask.vtk');
            Mat2VTK(f_vtk,Reshaped_mask,'binary');


        end
        
        % ----------------------------------------------------------------%
        % Function to update mask inside sheetlets

        function [cell_mask, cell_mask_id, sheetlet_wise_voxel_cell_id, cell_no_add] = update_inside_sheetlets(shtlt_wise_cell_list, sheetlet_count, Sheetletwise_voxel_list, NbPixel)
            
            % update_inside_sheetlets
            %
            % Generates the intracellular mask for cardiomyocytes located inside
            % individual sheetlets. For each sheetlet, voxel coordinates are tested
            % against cylindrical cell segments using collision detection. Voxels
            % intersecting cardiomyocyte geometry are assigned to the corresponding
            % cell and labeled using global cell IDs.
            %
            % Inputs:
            %   shtlt_wise_cell_list        - Cell array containing sheetlet-wise
            %                                 geometric descriptions of cell segments.
            %   sheetlet_count              - Number of sheetlets in the model.
            %   Sheetletwise_voxel_list     - Cell array containing voxel coordinates
            %                                 belonging to each sheetlet region.
            %   NbPixel                     - Number of voxels in each mask dimension.
            %
            % Outputs:
            %   cell_mask                   - Binary mask where ICV = 1 and ECV = 0.
            %   cell_mask_id                - Cell-ID mask where each intracellular
            %                                 voxel stores the corresponding global
            %                                 cell ID.
            %   sheetlet_wise_voxel_cell_id - Sheetlet-wise voxel-to-cell mapping.
            %   cell_no_add                 - Total number of intra-sheetlet cells,
            %                                 used as an offset for later branch IDs.
            
            % Innitialize mask 
            cell_mask = zeros(NbPixel); % mask values - [icv = 1, ecv = 0]
            cell_mask_id = zeros(NbPixel); % mask values - [icv = cellid, ecv = 0]
            
            % Initialize storage for voxel-to-cell mapping (sheetlet-wise)
            sheetlet_wise_voxel_cell_id = cell(sheetlet_count,1);
            
            % Update masks inside sheetlets
            cell_no_add = 0;
            for stlt = 1:sheetlet_count
            
                % Get cell objects and voxel coordinates for current sheetlet    
                objects = shtlt_wise_cell_list{stlt,1};
                voxels = Sheetletwise_voxel_list{stlt,1};
                points = voxels.*1e-6; %to convert voxel to voxel coordinate
                
                % Perform collision detection (voxel inside any cylindrical cell)
                [In_out_list,cell_id_list] = Mask_generation_Toolbox.Collison_detection(objects,points);
                
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


        end

        % ----------------------------------------------------------------%
        % Function to identify inasigned voxels in cleavege planes

        function filtered_voxel_list = filter_unassigned_voxels(NbPixel, Sheetletwise_mid_voxel_list, cell_mask)
            
            % filter_unassigned_voxels
            %
            % Identifies voxels that remain unassigned after generating the
            % intra-sheetlet cardiomyocyte mask. Mid-sheetlet voxels and already
            % occupied intracellular voxels are excluded, leaving candidate voxels
            % primarily located within cleavage-plane regions. These voxels are
            % later used for inter-sheetlet branch mask generation.
            %
            % Inputs:
            %   NbPixel                    - Number of voxels in each mask dimension.
            %   Sheetletwise_mid_voxel_list - Cell array containing voxels located
            %                                 inside mid-sheetlet regions.
            %   cell_mask                  - Binary mask after intra-sheetlet cells
            %                                have been assigned.
            %
            % Output:
            %   filtered_voxel_list        - List of unassigned voxel coordinates,
            %                                representing candidate cleavage-plane
            %                                voxels.


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
        end

        % ----------------------------------------------------------------%

        function [mask_before_reshape, mask_id_before_reshape, cell_id_list] = update_inter_sheetlet_branches(inter_shtlt_branch_list, cleavage_plane_voxels, n_insheetlet_cells, cell_mask, cell_mask_id)
            
            % update_inter_sheetlet_branches
            %
            % Generates the intracellular mask associated with inter-sheetlet
            % branches by evaluating candidate voxels located within cleavage-plane
            % regions. Voxels intersecting branch segments are assigned to the
            % corresponding branch and incorporated into the mask.
            %
            % Inputs:
            %   inter_shtlt_branch_list - Geometric description of inter-sheetlet
            %                             branch segments.
            %   cleavage_plane_voxels   - List of unassigned voxels located within
            %                             cleavage-plane regions.
            %   n_insheetlet_cells      - Total number of intra-sheetlet cells.
            %                             Used to convert local branch IDs into
            %                             global cell IDs.
            %   cell_mask              - Binary intracellular mask generated from
            %                            intra-sheetlet cardiomyocytes.
            %   cell_mask_id           - Cell-ID mask corresponding to cell_mask.
            %
            % Outputs:
            %   mask_before_reshape    - Updated binary mask including
            %                            inter-sheetlet branches.
            %   mask_id_before_reshape - Updated cell-ID mask including branch IDs.
            %   cell_id_list           - Local branch IDs returned by collision
            %                            detection.
            
            % Create mask for inter-sheetlet branches
            objects = inter_shtlt_branch_list;
            voxels = cleavage_plane_voxels;
            points = voxels.*1e-6; %to convert voxel to voxel coordinate

            % Collision detection for branches
            [In_out_list,cell_id_list] = Mask_generation_Toolbox.Collison_detection(objects,points);
            ICV_voxel_num = find(In_out_list);
            ICV_voxels = voxels(ICV_voxel_num,:);

            % Convert to global cell IDs
            cell_id_list_global = cell_id_list + n_insheetlet_cells;

            mask_before_reshape = cell_mask;
            mask_id_before_reshape = cell_mask_id;
            % Assign branch voxels into mask
            for vox=1:length(ICV_voxels)
                mask_before_reshape(ICV_voxels(vox,1),ICV_voxels(vox,2),ICV_voxels(vox,3)) = 1;
                vox_id = ICV_voxel_num(vox,1);
                mask_id_before_reshape(ICV_voxels(vox,1),ICV_voxels(vox,2),ICV_voxels(vox,3)) = cell_id_list_global(vox_id,1);
            end


        end

        % ----------------------------------------------------------------%
        % Function to create cell detail list from connectivity table and coordinates
        function Cell_list = gen_cell_detail_list(con_table,point_list,max_cell_dia)

            % gen_cell_detail_list
            % Generates a structured array describing each cell segment based
            % on connectivity and coordinate data. Each cell is approximated
            % as a cylinder.
            %
            % Inputs:
            %   con_table     - Connectivity table (cell array of node indices)
            %   point_list    - Cell array of 3D coordinates
            %   max_cell_dia  - Maximum cell diameter
            %
            % Output:
            %   Cell_list     - Nx9 matrix containing:
            %       [v1x v1y v1z v2x v2y v2z radius length volume]

            Cell_list = zeros(length(con_table),9);
            cell_dia = ones(length(con_table)).*max_cell_dia;
            for cell_id = 1:length(con_table)
                %Creating array contain cell data
                % column1,2,3 - vertex1 coordinates x,y,z
                % column4,5,6 - vertex2 coordinates x,y,z
                % column7 - radius of cell
                % column8 - length of cell
                % column9 - Volume of cell

                v1_id = con_table{cell_id}(1);
                v2_id = con_table{cell_id}(2);
                Cell_list(cell_id,1) = point_list{v1_id}(1)*1e-6;  %micro_m to m
                Cell_list(cell_id,2) = point_list{v1_id}(2)*1e-6;
                Cell_list(cell_id,3) = point_list{v1_id}(3)*1e-6;
                Cell_list(cell_id,4) = point_list{v2_id}(1)*1e-6;
                Cell_list(cell_id,5) = point_list{v2_id}(2)*1e-6;
                Cell_list(cell_id,6) = point_list{v2_id}(3)*1e-6;
                Cell_list(cell_id,7) = (cell_dia(cell_id)/2)*1e-6;
                v1_coord = point_list{v1_id};
                v2_coord = point_list{v2_id};
                cell_vec = v1_coord-v2_coord;
                Cell_list(cell_id,8) = norm(cell_vec)*1e-6;
                Cell_list(cell_id,9) = pi*(Cell_list(cell_id,7)^2)*Cell_list(cell_id,8);
            end
        end

             
        % ----------------------------------------------------------------%
        % Function to generate sheetlet wise cell list
        function shtlt_wise_cell_list = Sheetletwise_cell_details(shtltwise_contable,point_list,max_cell_dia)
            
            % Sheetletwise_cell_details
            % Generates cell detail lists for each sheetlet separately.
            %
            % Inputs:
            %   shtltwise_contable - Cell array of connectivity tables per sheetlet
            %   point_list         - Cell array of node coordinates
            %   max_cell_dia       - Maximum cell diameter
            %
            % Output:
            %   shtlt_wise_cell_list - Cell array of cell detail matrices
            
            shtlt_wise_cell_list = cell(length(shtltwise_contable),1);
            
            for i = 1:length(shtltwise_contable)
                con_list = shtltwise_contable{i};
                cell_det_list = Mask_generation_Toolbox.gen_cell_detail_list(con_list,point_list,max_cell_dia);
                shtlt_wise_cell_list{i} = cell_det_list;
            end
            
            
        end
        

        % ----------------------------------------------------------------%
        % Function to generate list of each voxel inside a sheetlet
        
        function Sheetletwise_voxel_list = Gen_sheetlet_voxel_list(sheetlets,NbPixel)

            % Gen_sheetlet_voxel_list
            % Generates voxel coordinates belonging to each sheetlet region.
            %
            % Inputs:
            %   sheetlets - Cell array of polygon definitions
            %   NbPixel   - Mask dimensions [Nx, Ny, Nz]
            %
            % Output:
            %   Sheetletwise_voxel_list - Cell array of voxel coordinate lists
            
            %identify_voxels_in_each_sheetlet
            sheetlet_voxels = Mask_generation_Toolbox.Sheetlet_boundary_mask(NbPixel,sheetlets);
            
            %Create sheetletwise woxel list
            Sheetletwise_voxel_list = cell(length(sheetlet_voxels),1);

            for i = 1:length(sheetlet_voxels)
                D2_list = sheetlet_voxels{1,i};
                num_voxels = size(D2_list,1)*NbPixel(1);
                voxel_list = zeros(num_voxels,3);

                vox_num = 1;
                for vox_X = 1:NbPixel(1)
                    for pix = 1:size(D2_list,1)
                        voxel_list(vox_num,:) = [vox_X,D2_list(pix,1),D2_list(pix,2)];
                        vox_num = vox_num+1;
                    end
                end
                Sheetletwise_voxel_list{i,1} = voxel_list;
            end
            
        end


        % ----------------------------------------------------------------%
        % function to identify sheetlet boundaries in the mask
        function [sheetlet_voxels] = Sheetlet_boundary_mask(mask_dim,polygons)

            % Sheetlet_boundary_mask
            % Identifies voxel coordinates that lie inside each sheetlet
            % polygon using 2D inpolygon checks.
            %
            % Inputs:
            %   mask_dim - Mask dimensions [Nx, Ny, Nz]
            %   polygons - Cell array of polygon vertex coordinates
            %
            % Output:
            %   sheetlet_voxels - Cell array of voxel coordinates per sheetlet

            %sheetlet_mask = zeros(mask_dim(1,2:end));
            sheetlet_voxels = {};
            %sheetlet_vox_boundaries = {};
            [YY,ZZ] = ndgrid(1:1:mask_dim(2),1:1:mask_dim(3));
            pixels = [YY(:) ZZ(:)];
            
            for shtlt = 1:length(polygons)
                polyV_y = polygons{shtlt}(:,1);
                polyV_z = polygons{shtlt}(:,2);
                pix_y = pixels(:,1);
                pix_z = pixels(:,2);
                
                [in,on] = inpolygon(pix_y,pix_z,polyV_y,polyV_z);
                pixels_in_poly = [pix_y(in),pix_z(in)];
                %pixels_on_poly = [pix_y(on),pix_z(on)];
                
                sheetlet_voxels{shtlt} = pixels_in_poly;
                %sheetlet_vox_boundaries{shtlt} = pixels_on_poly;
            end
            % plot 2d sheetlet structure
            Mask_2d = zeros(mask_dim(1,2:end));
            for i1 = 1:length(sheetlet_voxels)
                shtlt2 = sheetlet_voxels{i1};

                for j = 1:size(shtlt2,1)
                    Mask_2d(shtlt2(j,1),shtlt2(j,2)) = 1;
                end
            end
            figure
            s=imagesc(Mask_2d),colormap('gray')
        end

        % ----------------------------------------------------------------%
        % function to detect collison between object list and voxel list

        function [In_out_list,id_list] = Collison_detection(objects,points)

            % Collison_detection
            % Determines whether each voxel lies inside any cylindrical cell.
            %
            % Inputs:
            %   objects - Cell list (Nx7+) describing cylinders
            %   points  - Voxel coordinate list
            %
            % Outputs:
            %   In_out_list - Logical array (1 if inside any cell)
            %   id_list     - Index of intersecting cell

            size_p = size(points,1);
            In_out_list = false(size_p,1);
            id_list = zeros(size_p,1);
            
            parfor cpt_pixel = 1:1:size_p
                pixel_coord = points(cpt_pixel,:);

                for cell_index = 1:1:size(objects,1)
                    V1_vec = [objects(cell_index,1), objects(cell_index,2), objects(cell_index,3)];
                    V2_vec = [objects(cell_index,4), objects(cell_index,5), objects(cell_index,6)];

                    %Check the distance from cell centerline to the pixel coordinate
                    dis_to_pp = Mask_generation_Toolbox.min_dis_line_point(V1_vec,V2_vec,pixel_coord);

                    %if the minimum distance from cell centerline to pixel point is
                    %less than radius, that pixel is in the intra celluler volume
                    %(Represented by 1 in IN_OUT_LIST). It will break the loop if it is
                    %inside at least one cell

                    %Pixels in outside cell (In ECV) will be represented as a 0 in IN_OUT_LIST

                    if dis_to_pp < objects(cell_index,7)
                        In_out_list(cpt_pixel,1)=1;
                        id_list(cpt_pixel,1) = cell_index;
                        break;
                    else
                        continue;
                    end

                end
            end
        end
        
        % ----------------------------------------------------------------%
        % Function to calculate minimum distance between line and a point
        
        function dis_pt = min_dis_line_point(v1_line,v2_line,Pt)

            % min_dis_line_point
            % Computes minimum distance between a point and a line segment.
            % Only valid if projection lies within segment bounds.
            %
            % Inputs:
            %   v1_line, v2_line - Segment endpoints
            %   Pt               - Query point
            %
            % Output:
            %   dis_pt           - Minimum distance

            line_seg = v2_line - v1_line;
            unitV_line = line_seg/norm(line_seg);
            line_to_p = Pt - v1_line;

            %Find the projection of the point on line and check whether the point
            %is inside the line segment

            proj = dot(line_to_p,unitV_line)*unitV_line+v1_line;

            vec1 = v1_line - proj;
            vec2 = v2_line - proj;

            %if it is on line segment, return calculated distance
            if dot(vec1,vec2) < 0
                dis_pt = norm(cross(line_seg,line_to_p)) / norm(line_seg);
            else
                dis_pt = realmax;
            end
        end

        % ----------------------------------------------------------------%
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

        % ----------------------------------------------------------------%
        
       
    end
end