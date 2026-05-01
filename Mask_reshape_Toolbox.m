classdef Mask_reshape_Toolbox

    % Mask_reshape_Toolbox
    % Utility class for post-processing voxel masks, including inter-cell
    % distance calculations, slice-based corrections, and spatial perturbations.    
    methods(Static)
        
        %% function to create array including intercell distances between a given list of cells
        
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
        

        
        %% function to fill the ecv voxels of the previous slice which are filled by the same myocyte 
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
        
     %% Function to shake a slice of mask   
        
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
                    
        
        
        end
        
       
        %%
        
        
        
        
        
        
        
        
    end
    
end