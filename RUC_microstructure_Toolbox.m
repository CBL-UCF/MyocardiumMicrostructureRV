classdef RUC_microstructure_Toolbox

    methods(Static)

        function ruc_con = generate_ruc_myoconnectivity(RUC_w, RUC_t, RUC_l, shtlt_w, shtlt_t, shtlt_shrink_ratios, fiber_diameter, fiber_length_range, shtlt_ICV, max_fibredir_deviation, randomization_breakpoint)
            
            % ____________________________________________________________%
            % Generate sheetlet arrangement on first face and plot to
            % monitor microstructure constructuin
            ruc_con = 0;
            [Shrinked_after_inf, Shrinked_after_inf_bp, Shrinked_after_inf_mid, non_shrink_vor_ruc, Nn, Ns] = Sheetlet_arrangement_Toolbox.generate_sheetlet_arrangement(RUC_w, RUC_t, shtlt_w, shtlt_t, shtlt_shrink_ratios);

            fig1 = Plotting_Toolbox.monitor_sheetlet_boundaries_in3D(1, 'Monitor Microstructure Construction', Shrinked_after_inf, RUC_l);

            sheetlet_face1 = Shrinked_after_inf;
            sheetlet_face1_bp = Shrinked_after_inf_bp;
            
            % ____________________________________________________________%
            % Record coordinates of each vertex in the RUC and connectivity

            stlt_count = length(sheetlet_face1);
            rng_count = 0;

            coordinate_array = {};
            connectivity_array = {};
            connectivity_V1_array = {};
            connectivity_V2_array = {};
            face1_point_array = {};
            linewise_point_array = {};
            
            for stlt = 1: stlt_count
                
                stlt_vertex = sheetlet_face1{1,stlt};
                min_Y = min(stlt_vertex(:,1));
                min_Z = min(stlt_vertex(:,2));
                max_Y = max(stlt_vertex(:,1));
                max_Z = max(stlt_vertex(:,2));
            
                L_shtlt = RUC_l;        %length in micrometers(Along fibre dirrection)
                
                %Srink voronoi in 1:9 ratio to place base points of the extending
                %cardiomyocyte chains inside the voronoi
                
                stlt_bp_boundary = sheetlet_face1_bp{1,stlt};
                stlt_bp_boundary(end+1,:) = stlt_bp_boundary(1,:);
                Y_bp_vertex = transpose(stlt_bp_boundary(:,1));
                Z_bp_vertex = transpose(stlt_bp_boundary(:,2));
                
                target_ICV_sheetlet = 0.95; %TODO
                
                %c_sheetlet_innitial
                [face1_pt_list, connected_p_list, con_table_branched_final, graph_v1_branched, graph_v2_branched, linewise_point_list] = Sheetlet_microstructure_Toolbox.generate_insheetlet_microstructure(L_shtlt, ...
                min_Y, min_Z, max_Y, max_Z, fiber_diameter, fiber_length_range, shtlt_ICV, ...
                Y_bp_vertex, Z_bp_vertex, max_fibredir_deviation, target_ICV_sheetlet, randomization_breakpoint, rng_count);
                
                % Record microstructure inside each sheetlet
                
                coordinate_array{end+1} = connected_p_list;
                connectivity_array{end+1} = con_table_branched_final;
                connectivity_V1_array{end+1} = graph_v1_branched;
                connectivity_V2_array{end+1} = graph_v2_branched;
                face1_point_array{end+1} = face1_pt_list;
                linewise_point_array{end+1} = linewise_point_list;

                rng_count = rng_count+1;

            end







        end

    end

end