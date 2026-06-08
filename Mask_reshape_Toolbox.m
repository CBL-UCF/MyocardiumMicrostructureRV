classdef Mask_reshape_Toolbox

    % Mask_reshape_Toolbox
    %
    % Utility class for post-processing voxelized myocardial masks.
    %
    % The reshaping procedure is designed to improve cellular continuity
    % and produce more realistic cell boundaries after voxelization.
    %
    % Main operations:
    %   1. Recover disconnected intracellular voxels between adjacent slices.
    %   2. Introduce controlled slice perturbations ("shaking") to avoid
    %      unrealistic voxel alignment.
    %   3. Define cellular boundaries while preserving connectivity between
    %      physically connected cardiomyocytes.
    %   4. Generate masks suitable for diffusion simulations and image-based
    %      analyses.

    methods(Static)

        % ----------------------------------------------------------------%
        % function to perform first reshaping step (shaking slces)
        function [reshaped_cell_mask_id] = reshape_mask_shake_slices(cell_mask_id, reshaped1_rgnwise_cell_list, Rgn_wise_voxel_list)
            
            % reshape_mask_shake_slices
            %
            % Performs the first mask-reshaping stage by introducing small random
            % slice-wise perturbations and restoring cellular continuity between
            % neighboring slices.
            %
            % The procedure:
            %   1. Randomly shifts the next slice within a prescribed range.
            %   2. Compares the shifted slice with the current slice.
            %   3. Identifies voxels belonging to the same cardiomyocyte that are
            %      present in the shifted slice but absent in the current slice.
            %   4. Restores those voxels to improve continuity of cellular geometry.
            %
            % Inputs:
            %   cell_mask_id               - Cell-ID mask before reshaping.
            %   reshaped1_rgnwise_cell_list
            %                              - Region-wise cell information.
            %   Rgn_wise_voxel_list        - Region-wise voxel coordinates.
            %
            % Output:
            %   reshaped_cell_mask_id      - Cell-ID mask after slice-based
            %                                continuity correction.

            mask_dim = size(cell_mask_id);
            reshaped_temp_cell_mask_id = cell_mask_id;
            reshaped_cell_mask_id = cell_mask_id;
            
            for itr =1:6
                for slice = 1:mask_dim(1)-1
                    slice_pre = squeeze(reshaped_temp_cell_mask_id(slice,:,:));
                    slice_nxt = squeeze(reshaped_temp_cell_mask_id(slice+1,:,:));
            
                    shaked_slice_nxt = Mask_reshape_Toolbox.Shake_slice(slice_nxt,1,1,2);
            
            
                    [slice_pre_reshaped] = Mask_reshape_Toolbox.fill_prev_slice(slice,slice_pre,shaked_slice_nxt,Rgn_wise_voxel_list,reshaped1_rgnwise_cell_list);
            
                    reshaped_cell_mask_id(slice,:,:) =  slice_pre_reshaped;   
                end
                reshaped_temp_cell_mask_id = reshaped_cell_mask_id;
            end

        end

        % ----------------------------------------------------------------%

        function [reshaped_mask_id] = reshape_mask_define_cellboundary(shaked_mask_id, Global_connection_list, N_cells_before_branch)
            
            % reshape_mask_define_cellboundary
            %
            % Performs the second mask-reshaping stage by enforcing boundaries
            % between neighboring cardiomyocytes.
            %
            % Adjacent voxels are examined in both vertical and horizontal
            % directions. A boundary is inserted whenever neighboring voxels belong
            % to different cells that are not physically connected according to the
            % global connectivity graph.
            %
            % This prevents unrelated cardiomyocytes from sharing voxel boundaries
            % while preserving continuity across valid cell-cell connections and
            % branch locations.
            %
            % Inputs:
            %   shaked_mask_id        - Cell-ID mask after slice perturbation.
            %   Global_connection_list
            %                         - List of physically connected cells and
            %                           branches.
            %   N_cells_before_branch - Number of cardiomyocytes prior to branch
            %                           insertion.
            %
            % Output:
            %   reshaped_mask_id      - Final reshaped mask with cell boundaries
            %                           enforced.

            mask_dim = size(shaked_mask_id);
            reshaped_mask_id = shaked_mask_id;

            %add boundaries using verticle columns of pixels
            for slice = 1:1:mask_dim(1)
                slice_mask = squeeze(reshaped_mask_id(slice,:,:));
            
                slc_dim = size(slice_mask);
                for col = 1:slc_dim(2)-1
                    col_pre = slice_mask(:,col);
                    col_nxt = slice_mask(:,col+1);
                    
                    for i = 1:slc_dim(1)
                        pix1 = min(col_pre(i,1),col_nxt(i,1));
                        pix2 = max(col_pre(i,1),col_nxt(i,1));
                        
                        cons=[];
                        if pix1 < N_cells_before_branch && pix1>0
                            cons = Global_connection_list{pix1};
                        end
                        
                        if pix1==0 || pix2 ==0
                            continue;
                        elseif pix1 == pix2
                            continue;
                        elseif ismember(pix2,cons)
                            continue;
                        else
                            col_nxt(i,1) = 0;
                        end
                    end
                    slice_mask(:,col+1) = col_nxt;
                end
                reshaped_mask_id(slice,:,:) = slice_mask;
            end 

            %add boundaries using horizontal rows of pixels

            for slice = 1:1:mask_dim(1)
                slice_mask = squeeze(reshaped_mask_id(slice,:,:));
            
                slc_dim = size(slice_mask);
                for raw = 1:slc_dim(1)-1
                    raw_pre = slice_mask(raw,:);
                    raw_nxt = slice_mask(raw+1,:);
                    
                    for i = 1:slc_dim(2)
                        pix1 = min(raw_pre(1,i),raw_nxt(1,i));
                        pix2 = max(raw_pre(1,i),raw_nxt(1,i));
                        
                        cons=[];
                        if pix1 < N_cells_before_branch && pix1>0
                            cons = Global_connection_list{pix1};
                        end
                        
                        if pix1==0 || pix2 ==0
                            continue;
                        elseif pix1 == pix2
                            continue;
                        elseif ismember(pix2,cons)
                            continue;
                        else
                            raw_nxt(1,i) = 0;
                        end
                    end
                    slice_mask(raw+1,:) = raw_nxt;
                end
                reshaped_mask_id(slice,:,:) = slice_mask;
            end 
            
        end
        
        % ----------------------------------------------------------------%
        % function to create array including intercell distances between a given list of cells
        
        function [intercell_dis] = cell_distance_array(N_cells,cell_list)
            
            % cell_distance_array
            % Computes pairwise minimum distances between cylindrical cell
            % segments represented by their endpoints.
            %
            % Inputs:
            %   N_cells   - Number of cells
            %   cell_list - Nx6 matrix containing endpoints:
            %               [x1 y1 z1 x2 y2 z2]
            %
            % Output:
            %   intercell_dis - NxN symmetric matrix of minimum distances
            %                   between cell segments

            intercell_dis = zeros(N_cells,N_cells);
            
            for cell_1 = 1:N_cells
                for cell_2 = 1:N_cells
                    if cell_2 > cell_1
                        l1p1 = cell_list(cell_1,1:3);
                        l1p2 = cell_list(cell_1,4:6);
                        l2p1 = cell_list(cell_2,1:3);
                        l2p2 = cell_list(cell_2,4:6);
                        distance = Coordinate_Toolbox.min_dis_2lines(l1p1,l1p2,l2p1,l2p2);
                        if distance == 0.0 %to avoid calculating distance between interconeected cells
                            distance = realmax;
                        end
                        
                        intercell_dis(cell_1,cell_2) = distance;
                        intercell_dis(cell_2,cell_1) = distance;
                    elseif cell_1 == cell_2 %for the same cell
                        distance = realmax;
                        intercell_dis(cell_1,cell_2) = distance;
                    end
                end
            end
            
        end
        

        % ----------------------------------------------------------------%
        % function to fill the ecv voxels of the previous slice which are filled by the same myocyte 
        % in the next slice
        
        function [slice_pre_reshaped] = fill_prev_slice(slice_no,slice_pre,slice_nxt,Rgn_wise_voxel_list,Rgn_wise_cell_list)
            
            % fill_prev_slice
            % Corrects extracellular voxels (ECV) in a previous slice using
            % information from the next slice by propagating cell IDs where
            % continuity is detected.
            %
            % Inputs:
            %   slice_no              - Current slice index
            %   slice_pre             - Previous slice mask
            %   slice_nxt             - Next slice mask
            %   Rgn_wise_voxel_list   - Region/sheetlet-wise voxel lists
            %   Rgn_wise_cell_list    - Region-wise cell data
            %
            % Output:
            %   slice_pre_reshaped    - Updated previous slice mask

            slice_change = slice_nxt-slice_pre;
            
            slice_pre_reshaped = slice_pre;
            
            sheetlet_count = length(Rgn_wise_voxel_list);
            for shtlt =1:sheetlet_count
                voxels = Rgn_wise_voxel_list{shtlt,1};
                cell_list = Rgn_wise_cell_list{shtlt,1}(:,1:9); 
                stlt_cell_id_list = cell_list(:,9);

                shtlt_slice_id = find(voxels(:,1)==slice_no);
                shtlt_slice = voxels(shtlt_slice_id,2:3);
                n_vox = size(shtlt_slice,1);
                shtlt_change = zeros(n_vox,1);
                for k=1:n_vox
                    shtlt_change(k,1) = slice_change(shtlt_slice(k,1),shtlt_slice(k,2));
                end
                unique_changes = unique(shtlt_change);
                slice_cells = intersect(stlt_cell_id_list,unique_changes);

                for cell = 1:size(slice_cells,1)
                    cell_id = slice_cells(cell,1);
                    changed_pix_ids = find(shtlt_change==cell_id);
                    changed_pix = shtlt_slice(changed_pix_ids,:);
                    n_ch_pix = size(changed_pix,1);
                    for p = 1:n_ch_pix
                        slice_pre_reshaped(changed_pix(p,1),changed_pix(p,2)) = cell_id;
                    end
                end
            end
            
        end
        % ----------------------------------------------------------------%
        % Function to shake a slice of mask   
        
        function [shaked_slice] = Shake_slice(mask_slice,horizontal,vertical,range)

            % Shake_slice
            % Applies a random spatial shift (translation) to a 2D mask slice
            % in horizontal and/or vertical directions.
            %
            % Inputs:
            %   mask_slice - 2D mask array
            %   horizontal - Boolean flag to enable horizontal shift
            %   vertical   - Boolean flag to enable vertical shift
            %   range      - Maximum translation magnitude (pixels)
            %
            % Output:
            %   shaked_slice - Translated mask slice
        
            dim = size(mask_slice);
            shaked_slice = zeros(dim);
            hshaked_slice = zeros(dim);
            trans = -1*range + randi([1,2*range],1,50);
            non_zero_trans = find(trans);
            cc = length(non_zero_trans);
            h_trans = trans(non_zero_trans(randi([1,cc],1,1)));
            v_trans = trans(non_zero_trans(randi([1,cc],1,1)));
            
            if horizontal
                if h_trans >= 0
                    for col = 1:1:dim(1)-h_trans
                        hshaked_slice(col+h_trans,:) = mask_slice(col,:);
                    end
                else
                    for col = abs(h_trans)+1:1:dim(1)
                        hshaked_slice(col+h_trans,:) = mask_slice(col,:);
                    end
                end
            else
                hshaked_slice = mask_slice;
            end
            
            if vertical
                if v_trans >= 0
                    for raw = 1:1:dim(2)-v_trans
                        shaked_slice(:,raw+v_trans) = hshaked_slice(:,raw);
                    end
                else
                    for raw = abs(v_trans)+1:1:dim(2)
                        shaked_slice(:,raw+v_trans) = hshaked_slice(:,raw);
                    end
                end
            else
                shaked_slice = hshaked_slice;
            end
            % ----------------------------------------------------------------%    
        
        
        end     
                
        
    end
    
end