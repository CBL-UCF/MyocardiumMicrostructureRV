%Mask generate toolbox

classdef Mask_generate_Toolbox

    % Mask_generate_Toolbox
    % Utility class for converting geometric representations of cells and
    % sheetlets into voxelized masks, including collision detection and
    % voxel classification.
    
    methods(Static)
        
        %% Function to create cell detail list from connectivity table and coordinates
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
        
        
        %% Function to generate sheetlet wise cell list
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
                cell_det_list = Mask_generate_Toolbox.gen_cell_detail_list(con_list,point_list,max_cell_dia);
                shtlt_wise_cell_list{i} = cell_det_list;
            end
            
            
        end
              
        %% Function to generate list of each voxel inside a sheetlet
        
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
            sheetlet_voxels = Mask_generate_Toolbox.Sheetlet_boundary_mask(NbPixel,sheetlets);
            
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
        
        %% function to identify sheetlet boundaries in the mask
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
        
        %% function to detect collison between object list and voxel list
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
                    dis_to_pp = Mask_generate_Toolbox.min_dis_line_point(V1_vec,V2_vec,pixel_coord);

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
        
        %% Function to calculate minimum distance between line and a point
        
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
        
       
    end
end