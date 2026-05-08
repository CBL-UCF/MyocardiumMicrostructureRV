classdef Sheetlet_arrangement_Toolbox

    % Class - Sheetlet_arrangement_Toolbox
    % Static utility class for generating and processing sheetlet
    % arrangements using Voronoi/Thiessen tessellation
    
    methods(Static)

        function [Shrinked_after_inf, Shrinked_after_inf_bp, Shrinked_after_inf_mid, non_shrink_vor_ruc, Nn, Ns, seed_array, vor_vertex] = generate_sheetlet_arrangement( ...
                RUC_w, RUC_t, shtlt_w, shtlt_t, shtlt_shrink_ratios)

            % Function - generate_sheetlet_arrangement
            % Generates a sheetlet arrangement within a representative unit cell
            % using Voronoi tessellation, then shrinks and filters the resulting
            % polygons to define sheetlet boundaries and related regions.
            %
            % Inputs:
            %   RUC_w                - Representative unit cell width.
            %   RUC_t                - Representative unit cell thickness.
            %   shtlt_w              - Approximate sheetlet width.
            %   shtlt_t              - Approximate sheetlet thickness.
            %   shtlt_shrink_ratios  - Shrink ratios used to create cleavage
            %                          plane space and auxiliary sheetlet regions.
            %
            % Outputs:
            %   Shrinked_after_inf      - Shrunk bounded sheetlet polygons.
            %   Shrinked_after_inf_bp   - Shrunk polygons used for base-point
            %                             generation.
            %   Shrinked_after_inf_mid  - Shrunk polygons used for mid-sheetlet
            %                             mask construction.
            %   non_shrink_vor_ruc      - Original Voronoi cells after removing
            %                             the outer polygon layer.\
            %   Nn                      - No fo sheetlets in normal direction
            %   Ns                      - No of sheetlets in sheetlet direction

            %_____________________________________________________________%            
            % Generate seed points for voronoi tesselation 
            
            [seed_points, Nn, Ns] = Sheetlet_arrangement_Toolbox.slt_mid_points(RUC_w,RUC_t, shtlt_w, shtlt_t);
            seed_array = []; % Convert cell to array
            for i=1:length(seed_points)
                seed_array(i,1)=seed_points{i}(1);
                seed_array(i,2)=seed_points{i}(2);
            end

            %_____________________________________________________________%
            % Generate thiessen polygons            
            [Vx,Vy] = voronoi(seed_array(:,1),seed_array(:,2));
            [vor_vertex,vor_cells] = voronoin(seed_array);

            %_____________________________________________________________%
            % Shrink thiessen polygons to make room for cleavage planes

            [shrink_ratio_boundary_y, shrink_ratio_boundary_z, shrink_ratio_basepoints_y, shrink_ratio_basepoints_z, shrink_ratio_mid_y, shrink_ratio_mid_z] = Sheetlet_arrangement_Toolbox.define_shrink_ratios(shtlt_shrink_ratios);
            
            Shrinked_tess = Sheetlet_arrangement_Toolbox.shrink_vor_cells(seed_array, ...
                vor_vertex, vor_cells, shrink_ratio_boundary_y,shrink_ratio_boundary_z);

            Shrinked_tess_bp = Sheetlet_arrangement_Toolbox.shrink_vor_cells(seed_array, ...
                vor_vertex,vor_cells,shrink_ratio_basepoints_y,shrink_ratio_basepoints_z);
            
            Shrinked_tess_mid = Sheetlet_arrangement_Toolbox.shrink_vor_cells(seed_array, ...
                vor_vertex,vor_cells,shrink_ratio_mid_y,shrink_ratio_mid_z);
            
            %_____________________________________________________________%
            % Remove outer layer of thiessen polygons (to remove skewed polygons
            % on the outer layer of polygons generated from 'voronoin' function)

            Shrinked_in_ruc = Sheetlet_arrangement_Toolbox.remove_outer_polylayer(Shrinked_tess, Nn); % Include vertex coordinates
            Shrinked_in_ruc_bp = Sheetlet_arrangement_Toolbox.remove_outer_polylayer(Shrinked_tess_bp, Nn);
            Shrinked_in_ruc_mid = Sheetlet_arrangement_Toolbox.remove_outer_polylayer(Shrinked_tess_mid, Nn);
            non_shrink_vor_ruc = Sheetlet_arrangement_Toolbox.remove_outer_polylayer(vor_cells, Nn);
            
            %_____________________________________________________________%
            % Remove polygons with 'infinity' as a vertex in case if
            % remained after removing outer layer
            Shrinked_after_inf = Sheetlet_arrangement_Toolbox.remove_unbounded_polygons(Shrinked_in_ruc);
            Shrinked_after_inf_bp = Sheetlet_arrangement_Toolbox.remove_unbounded_polygons(Shrinked_in_ruc_bp);
            Shrinked_after_inf_mid = Sheetlet_arrangement_Toolbox.remove_unbounded_polygons(Shrinked_in_ruc_mid);

        end

        % ----------------------------------------------------------------%
        function [mid_point_list, Nn, Ns] = slt_mid_points(RUC_w,RUC_t,Slt_w, Slt_t)
            
            % Function - slt_mid_points
            % Generates approximate sheetlet seed/mid-point coordinates within
            % the representative unit cell. These points are later used as seeds
            % for the Voronoi tessellation.
            %
            % Inputs:
            %   RUC_w  - Representative unit cell width.
            %   RUC_t  - Representative unit cell thickness.
            %   Slt_w  - Approximate sheetlet width.
            %   Slt_t  - Approximate sheetlet thickness.
            %
            % Outputs:
            %   mid_point_list - Cell array containing [y,z] seed coordinates.
            %   Nn             - Number of seed points in the normal direction.
            %   Ns             - Number of seed points in the sheetlet direction.

            mid_point_list = {};
            
            Ns = round(RUC_w/(Slt_w*1.1))+2;
            Nn = round(RUC_t/(Slt_t*1.1))+2;          
            
            for pts = 1:Ns
                for ptn = 1:Nn
                    if rem(pts,2)==1
                        z=0;
                    else
                        z = 0;%0.02*Slt_t; %0.05*Slt_t;
                    end
                    if rem(ptn,2)==1
                        y=0;
                    else
                        y = 0; %0.02*Slt_w; %0.5*Slt_t;
                    end
                    if rem(ptn,2)==1
                        y=y+(pts+0.46+0.08*rand())*Slt_w;
                        z=z+(ptn+0.48+0.04*rand())*Slt_t;
                    else
                        y=y+(0+pts+0.46+0.08*rand())*Slt_w; %add 0.5 of sheetlet width to shift the layer
                        z=z+(ptn+0.48+0.04*rand())*Slt_t;
                    end
                    point = [y,z];
                    mid_point_list{end+1} = point;
                end
            end
        end

        % ----------------------------------------------------------------%

        function [shrink_ratio_boundary_y, shrink_ratio_boundary_z, shrink_ratio_basepoints_y, shrink_ratio_basepoints_z, shrink_ratio_mid_y, shrink_ratio_mid_z] = define_shrink_ratios(shtlt_shrink_ratios)
            
            % Function - define_shrink_ratios
            % Defines shrink ratios for multiple sheetlet-related regions.
            % These ratios control how Voronoi polygons are contracted toward
            % their corresponding seed points in the y and z directions.
            %
            % Input:
            %   shtlt_shrink_ratios - Base shrink ratios for the y and z
            %                         directions.
            %
            % Outputs:
            %   shrink_ratio_boundary_y    - y-direction shrink ratio for
            %                                sheetlet boundaries.
            %   shrink_ratio_boundary_z    - z-direction shrink ratio for
            %                                sheetlet boundaries.
            %   shrink_ratio_basepoints_y  - y-direction shrink ratio for
            %                                cardiomyocyte base-point regions.
            %   shrink_ratio_basepoints_z  - z-direction shrink ratio for
            %                                cardiomyocyte base-point regions.
            %   shrink_ratio_mid_y         - y-direction shrink ratio for
            %                                mid-sheetlet regions.
            %   shrink_ratio_mid_z         - z-direction shrink ratio for
            %                                mid-sheetlet regions.
            r = shtlt_shrink_ratios;

            % Define sheetlet boundary ratios
            shrink_ratio_boundary_y = [r(1), 1-r(1)]; 
            shrink_ratio_boundary_z = [r(2), 1-r(2)];

            % Define boundary ratios for seeding cardiomyocyte mid points
            shrink_ratio_basepoints_y = [r(1)-0.050, 1-r(1)+0.050]; 
            shrink_ratio_basepoints_z = [r(2)-0.075, 1-r(2)+0.075];

            % Defnie mid sheetlet boundary rations for mask construction procedure
            shrink_ratio_mid_y = [r(1)-0.250, 1-r(1)+0.250]; 
            shrink_ratio_mid_z = [r(2)-0.250, 1-r(2)+0.250];

        end

        % ----------------------------------------------------------------%

        function shrinked_tess_coord = shrink_vor_cells(seed_array,v,vor_cells,shrink_ratio_y,shrink_ratio_z)
            
            % Function - shrink_vor_cells
            % Shrinks each Voronoi polygon toward its corresponding seed point
            % using independent shrink ratios in the y and z directions.
            %
            % Inputs:
            %   seed_array      - Array of seed point coordinates, where each
            %                     row contains [y,z].
            %   v               - Voronoi vertex coordinate array returned by
            %                     voronoin.
            %   vor_cells       - Cell array of Voronoi cells returned by
            %                     voronoin.
            %   shrink_ratio_y  - Two-component shrink ratio in the y direction.
            %   shrink_ratio_z  - Two-component shrink ratio in the z direction.
            %
            % Output:
            %   shrinked_tess_coord - Cell array containing shrunk polygon
            %                         vertex coordinates.
            
            shrinked_tess_coord ={};
        
            for vor = 1:length(seed_array)
                
                seed_y = seed_array(vor,1);
                seed_z = seed_array(vor,2);
                
                ex_poly = vor_cells{vor};
                ex_poly_coord = [];
                
                for i=1:length(ex_poly)
                    ex_poly_coord(i,1) = v(ex_poly(i),1);
                    ex_poly_coord(i,2) = v(ex_poly(i),2);
                end
                
                ex_poly_coord;
                aa_y = shrink_ratio_y(1);
                bb_y = shrink_ratio_y(2);
                aa_z = shrink_ratio_z(1);
                bb_z = shrink_ratio_z(2);
                nw_poly_coord =[];
                
                for point_v = 1:length(ex_poly_coord)
                    nw_poly_coord(point_v,1) = (bb_y*seed_array(vor,1)+aa_y*ex_poly_coord(point_v,1))/(aa_y+bb_y);
                    nw_poly_coord(point_v,2) = (bb_z*seed_array(vor,2)+aa_z*ex_poly_coord(point_v,2))/(aa_z+bb_z);
                end
                
                shrinked_tess_coord{end+1} = nw_poly_coord;
            end
        end

        % ----------------------------------------------------------------%

        function filtered_polygons = remove_outer_polylayer(all_polygons, Nn)
            
            % Function - remove_outer_polylayer
            % Removes the outer layer of Voronoi polygons. This is used to
            % exclude skewed or boundary-affected polygons produced by the
            % Voronoi tessellation.
            %
            % Inputs:
            %   all_polygons - Cell array containing polygon data.
            %   Nn           - Number of polygons in the normal direction.
            %
            % Output:
            %   filtered_polygons - Cell array after removing the outer polygon
            %                       layer.

            filtered_polygons = {};
            for p_no = 1: length(all_polygons)
                if p_no < Nn
                    continue
                elseif p_no > length(all_polygons)-Nn
                    continue
                elseif rem(p_no,Nn) == 0
                    continue
                elseif rem(p_no,Nn) == 1
                    continue
                else 
                   filtered_polygons{end+1} =  all_polygons{p_no};
                end
            end

        end

        % ----------------------------------------------------------------%

        function filtered_polygons = remove_unbounded_polygons(all_polygons)
            
            % Function - remove_unbounded_polygons
            % Removes polygons that contain infinite coordinates. These
            % correspond to unbounded Voronoi cells and are not suitable for
            % finite representative unit cell construction.
            %
            % Input:
            %   all_polygons - Cell array containing polygon coordinate arrays.
            %
            % Output:
            %   filtered_polygons - Cell array containing only bounded polygons.

            filtered_polygons = {};
            for poly = 1:length(all_polygons)
                coords_poly = all_polygons{poly};
                inf_indexes = isinf(coords_poly); %return 'true' for corresponding matrix index if the element is inf'
                N_inf = nonzeros(inf_indexes);
                %considered double have infinity coordinates if N_inf has any entities
                
                if ~isempty(N_inf)
                    continue
                else
                    filtered_polygons{end+1} = all_polygons{poly};
                end        
            end

        end

        % ----------------------------------------------------------------%
    end

end