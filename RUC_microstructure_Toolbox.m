classdef RUC_microstructure_Toolbox

    % RUC_microstructure_Toolbox
    %
    % Utility class for generating representative unit cell (RUC)
    % myocardial microstructures. The workflow includes:
    %   - Sheetlet arrangement generation
    %   - Intra-sheetlet cardiomyocyte connectivity generation
    %   - Global connectivity assembly
    %   - Neighbouring sheetlet identification
    %   - Common-edge detection
    %   - Inter-sheetlet branch generation
    %   - Branch filtering to achieve target intracellular volume fraction
    %
    % The generated microstructure can subsequently be voxelized and used
    % for diffusion simulations or mechanical modeling.

    methods(Static)

        function [RUC_con_table_branched, RUC_v_coord, RUC_con_v1_branched, RUC_con_v2_branched, Final_branches, non_shrink_vor_ruc, RUC_con_table_sheetlet_wise, Shrinked_after_inf, Shrinked_after_inf_mid] = generate_ruc_myoconnectivity(RUC_w, RUC_t, RUC_l, shtlt_w, shtlt_t, shtlt_shrink_ratios, fiber_diameter, fiber_length_range, shtlt_ICV, appx_target_ICV, max_fibredir_deviation, branching_reg_thickness, randomization_breakpoint)
            
            % generate_ruc_myoconnectivity
            %
            % Main workflow for generating a complete myocardial RUC.
            %
            % Steps:
            %   1. Generate sheetlet arrangement.
            %   2. Generate cardiomyocyte connectivity within each sheetlet.
            %   3. Assemble sheetlet-wise connectivity into a global graph.
            %   4. Identify neighbouring sheetlets.
            %   5. Detect common sheetlet boundaries.
            %   6. Generate potential inter-sheetlet branches.
            %   7. Select branches to achieve the target volume fraction.
            %   8. Return the final branched connectivity representation.

            % ____________________________________________________________%
            % Generate sheetlet arrangement on first face and plot to
            % monitor microstructure constructuin
            ruc_con = 0;
            [Shrinked_after_inf, Shrinked_after_inf_bp, Shrinked_after_inf_mid, non_shrink_vor_ruc, Nn, Ns, seed_array, vor_vertex] = Sheetlet_arrangement_Toolbox.generate_sheetlet_arrangement(RUC_w, ...
                RUC_t, shtlt_w, shtlt_t, shtlt_shrink_ratios);

            fig1 = Plotting_Toolbox.monitor_sheetlet_boundaries_in3D(1, 'Monitor Microstructure Construction', Shrinked_after_inf, RUC_l);

            sheetlet_face1 = Shrinked_after_inf;
            sheetlet_face1_bp = Shrinked_after_inf_bp;
            
            % ____________________________________________________________%
            % Generate microstructure inside each sheetlet region
            
            % Generate sheetletwise connectivity
            [coordinate_array, connectivity_array, connectivity_V1_array, connectivity_V2_array, face1_point_array, linewise_point_array] = RUC_microstructure_Toolbox.generate_sheetlet_wise_connectivity(sheetlet_face1, ...
                sheetlet_face1_bp, RUC_l, shtlt_ICV, fiber_diameter, ...
                fiber_length_range, max_fibredir_deviation, randomization_breakpoint);
            
            % Get RUC cardiomyocyte connectivity to a single graph
            % definition
            [RUC_v_coord, RUC_con_table, RUC_con_table_sheetlet_wise, RUC_con_v1, RUC_con_v2] = RUC_microstructure_Toolbox.Define_single_connectivity_for_RUC(coordinate_array, connectivity_array);

            % Plot sheetlet connectivity before inter sheetlet branches
            fig2 = Plotting_Toolbox.plot_microstructure(2, RUC_v_coord, RUC_con_v1, RUC_con_v2);

            % ____________________________________________________________%
            % Generate inter sheetlet branches

            % Identify neighbouring sheetlets for each sheetlet
            [Neghbour_stlt__rucid] = RUC_microstructure_Toolbox.identify_neighbouring_sheetlets(seed_array, Nn, Ns);

            % Generate common edge array for RUC
            min_edgelength = 30; % TODO
            [common_edge_array_ruc] = RUC_microstructure_Toolbox.find_common_edges(Neghbour_stlt__rucid, non_shrink_vor_ruc, min_edgelength, vor_vertex);
            
            % Define inter-sheetlet branches
            [Possible_Edge_wise_branches] = RUC_microstructure_Toolbox.find_potential_intersheetlet_branches(RUC_v_coord, ...
                common_edge_array_ruc, face1_point_array, coordinate_array, connectivity_array, branching_reg_thickness, linewise_point_array, fiber_diameter/2);

            % Concaternate potential branch list for RUC
            Possible_ruc_branches = {};
            for edge_no = 1:length(Possible_Edge_wise_branches)
                b_s = Possible_Edge_wise_branches{edge_no};
                for br = 1:length(b_s)
                    Possible_ruc_branches{end+1} = b_s{br};
                end
            end

            %_____________________________________________________________%

            % Compute volume fractions before inter-sheetlet branching

            Total_area = RUC_microstructure_Toolbox.Poly_area_calc(non_shrink_vor_ruc,vor_vertex);
            Total_volume = Total_area*RUC_l;
            Sheetlet_myocyte_vol = Sheetlet_microstructure_Toolbox.Myo_vol_calc(RUC_con_table,RUC_v_coord,fiber_diameter/2);
            Vf_without_branching = Sheetlet_myocyte_vol/Total_volume;
            Possible_branch_vol = Sheetlet_microstructure_Toolbox.Myo_vol_calc(Possible_ruc_branches,RUC_v_coord,fiber_diameter/2);
            max_vol_frac = (Sheetlet_myocyte_vol+Possible_branch_vol)/Total_volume;


            %_____________________________________________________________%\
            % Select branches randomly from potential branches until
            % approximate volume fraction is achieved

            Final_branches = RUC_microstructure_Toolbox.filter_final_branches(appx_target_ICV, max_vol_frac, Total_volume, Vf_without_branching, Possible_ruc_branches, RUC_v_coord, fiber_diameter/2);

            %_____________________________________________________________%
            % Define branched connectivity table

            RUC_con_v1_branched = RUC_con_v1;
            RUC_con_v2_branched = RUC_con_v2;
            RUC_con_table_branched = RUC_con_table;  %connectivity table (type-cell)
            
            
            % for branch_no1 = 1:length(Final_branches)
            %     RUC_con_v1_branched(end+1) = Final_branches{branch_no1}(1);
            %     RUC_con_v2_branched(end+1) = Final_branches{branch_no1}(2);
            % end
            % for branch_no2 = 1:length(Final_branches)
            %     RUC_con_table_branched{end+1} = Final_branches{branch_no2};
            % end
        
        end

        % ----------------------------------------------------------------%
        function [coordinate_array, connectivity_array, connectivity_V1_array, connectivity_V2_array, face1_point_array, linewise_point_array] = generate_sheetlet_wise_connectivity(sheetlet_face1, ...
                sheetlet_face1_bp, RUC_l, shtlt_ICV, fiber_diameter, fiber_length_range, max_fibredir_deviation, randomization_breakpoint)
            
            % generate_sheetlet_wise_connectivity
            %
            % Generates cardiomyocyte microstructure independently inside
            % each sheetlet region. Connectivity, graph representation,
            % base points, and line-wise point chains are stored separately
            % for each sheetlet.

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
                
                
                
                %c_sheetlet_innitial
                [face1_pt_list, connected_p_list, con_table_branched_final, graph_v1_branched, graph_v2_branched, linewise_point_list] = Sheetlet_microstructure_Toolbox.generate_insheetlet_microstructure(L_shtlt, ...
                min_Y, min_Z, max_Y, max_Z, fiber_diameter, fiber_length_range, shtlt_ICV, ...
                Y_bp_vertex, Z_bp_vertex, max_fibredir_deviation, randomization_breakpoint, rng_count);
                
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

        %-----------------------------------------------------------------%

        function [RUC_v_coord, RUC_con_table, RUC_con_table_sheetlet_wise, RUC_con_v1, RUC_con_v2] = Define_single_connectivity_for_RUC(coordinate_array, connectivity_array)
            
            % Define_single_connectivity_for_RUC
            %
            % Concatenates all sheetlet-wise coordinate and connectivity
            % information into a single global graph representation of the
            % representative unit cell.
            %
            % Outputs:
            %   RUC_v_coord              - Global node coordinate list.
            %   RUC_con_table            - Global connectivity table.
            %   RUC_con_table_sheetlet_wise
            %                           - Sheetlet-wise connectivity.
            %   RUC_con_v1, RUC_con_v2  - Graph edge definitions.

            % Concaternate coordinate list
            RUC_v_coord = {};

            for stlt1 = 1:length(coordinate_array)
                coords_stlt = coordinate_array{1,stlt1};
                for coorda = 1:length(coords_stlt)
                    RUC_v_coord{end+1} = coords_stlt{coorda};
                end
            end

            % redefine connectivity table to comply with concaternated
            % point coordinates

            RUC_con_table = {};
            RUC_con_table_sheetlet_wise = {};
            
            
            for stlt2 = 1:length(connectivity_array)
            if stlt2 == 1
                node_num_addition = 0;
            else
                node_num_addition = node_num_addition + length(coordinate_array{1,stlt2-1});
            end
                
            con_stlt = connectivity_array{1,stlt2};
            for_sheetlet = {};
            for con = 1:length(con_stlt)
                 RUC_con_table{end+1} = con_stlt{con} + node_num_addition;
                 for_sheetlet{end+1} = con_stlt{con} + node_num_addition;
            end
            RUC_con_table_sheetlet_wise{stlt2} = for_sheetlet;
            end
                
            RUC_con_v1 = [];
            RUC_con_v2 = [];
            
            for myo = 1:length(RUC_con_table)
             RUC_con_v1(end+1) = RUC_con_table{myo}(1);
             RUC_con_v2(end+1) = RUC_con_table{myo}(2);
            end
       
        end

        % ----------------------------------------------------------------%

        function [Neghbour_stlt__rucid] = identify_neighbouring_sheetlets(seed_array, Nn, Ns)
            
            % identify_neighbouring_sheetlets
            %
            % Uses Delaunay triangulation of sheetlet seed points to
            % identify neighbouring sheetlets within the representative
            % unit cell.
            %
            % Neighbor relationships are re-indexed to correspond to the
            % sheetlet numbering used after removal of outer and unbounded
            % Voronoi regions.

            DT_ruc_before = delaunayTriangulation(seed_array);

            seed_con = DT_ruc_before.ConnectivityList;
    
            all_seed_ids = (1:1:length(seed_array));
            
            % Filter sheetlet ids from innitial voronoi tesselation which
            % are inside the RUC
            RUC_seed_ids =[];
            for p_id = 1: length(all_seed_ids)
                if p_id < Nn
                    continue
                elseif p_id > length(all_seed_ids)-Nn
                    continue
                elseif rem(p_id,Nn) == 0
                    continue
                elseif rem(p_id,Nn) == 1
                    continue
                else 
                   RUC_seed_ids(end+1) =  all_seed_ids(p_id);
                end
            end

            % Find neighbours for sheetlet inside the RUC
            Neighbour_stlt_array = {};
    
            for stlt = 1:length(RUC_seed_ids)
                neigbours =[];
                sheetlet_no = RUC_seed_ids(stlt);
                corr_p_ids = find(seed_con == sheetlet_no); % ids of the considered sheetlet seed points in the delunui connectivity array
                
                uu = length(seed_con);
                for inin = 1:length(corr_p_ids)
                    id1 = corr_p_ids(inin);
                    if id1<= uu
                        neigbours(end+1) = seed_con(id1+uu);
                        neigbours(end+1) = seed_con(id1+2*uu);
                    elseif id1> uu && id1<= 2*uu
                        neigbours(end+1) = seed_con(id1-uu);
                        neigbours(end+1) = seed_con(id1+uu);
                    elseif id1 > 2*uu
                        neigbours(end+1) = seed_con(id1-uu);
                        neigbours(end+1) = seed_con(id1-2*uu);
                    end
                end
                Neighbour_stlt_array{end+1} = neigbours;
            end

            % Filter neighbouring sheetlets which are included in RUC
            Neghbour_stlt_array_ruc = {}; %array of sheetlet ids only in RUC
    
            %From connected points to each vertexes, select only points that are in RUC
            %and neglect others
            for stlt = 1:length(Neighbour_stlt_array)
                ruc_neighbour = intersect(Neighbour_stlt_array{stlt},RUC_seed_ids);
                Neghbour_stlt_array_ruc{end+1} = ruc_neighbour;
            end

            %Above array indexes are correseponding to seed point id before removing
            %outer layer of sheetlets and unbounded Thiessen polygons.
            %Now, we redefine the same array using sheetlet ids used in RUC
            %connectivity 
            
            Neghbour_stlt__rucid = {};
            for stlt6 = 1:length(Neghbour_stlt_array_ruc)
                ruc_neighbour1 = [];
                for vv = 1:length(Neghbour_stlt_array_ruc{stlt6})
                    ruc_neighbour1(end+1) = find(RUC_seed_ids==Neghbour_stlt_array_ruc{stlt6}(vv));
                end
                Neghbour_stlt__rucid{end+1} = ruc_neighbour1;
            end

        end

        %-----------------------------------------------------------------%

        function  [common_edge_array_ruc] = find_common_edges(Neghbour_stlt__rucid, non_shrink_vor_ruc, min_edgelength, vor_vertex)
            
            % find_common_edges
            %
            % Determines shared boundaries between neighbouring sheetlets.
            %
            % Outputs an array containing:
            %   - Sheetlet pair IDs
            %   - Shared edge coordinates
            %   - Edge length
            %
            % Edges shorter than a specified threshold are discarded.

            %identifying shared edges
            common_edges_ruc = {};
            for stlt7 = 1: length(Neghbour_stlt__rucid)
                common_edge={};
                n_cells = Neghbour_stlt__rucid{stlt7};
                for n_cell_id = 1:length(n_cells)
                    common_edge{end+1} = intersect(non_shrink_vor_ruc{stlt7},non_shrink_vor_ruc{n_cells(n_cell_id)});
                end
                common_edges_ruc{end+1} = common_edge;
            end

            %Creating an array that include all the common edges inside the RUC and
            %their corresponding sheetlet numbers
            %Arrangement:-
            % shtlt1 | shtlt2 | commonedge1(Y) | commonedge1(Z) | commonedge2(Y) | commonedge2(Z) | edge length
            
            common_edge_array_stlt_wise =[];
            edge_no = 1;
            for stlt8 = 1: length(common_edges_ruc)
                com_b = common_edges_ruc{stlt8};
                for edgeno = 1:length(com_b)
                    if length(com_b{edgeno}) < 2 %to avoid diagonal neighbours which have only one common vertex(not a edge)
                        continue;
                    else
                        common_edge_array_stlt_wise(edge_no,1) = stlt8;
                        common_edge_array_stlt_wise(edge_no,2) = Neghbour_stlt__rucid{stlt8}(edgeno);
            
                        vertex1_id = com_b{edgeno}(1);
                        common_edge_array_stlt_wise(edge_no,3) = vor_vertex(vertex1_id,1);
                        common_edge_array_stlt_wise(edge_no,4) = vor_vertex(vertex1_id,2);
            
                        vertex2_id = com_b{edgeno}(2);
                        common_edge_array_stlt_wise(edge_no,5) = vor_vertex(vertex2_id,1);
                        common_edge_array_stlt_wise(edge_no,6) = vor_vertex(vertex2_id,2); 
            
                        edge_no = edge_no + 1;
                    end
                    
                end
                
            end

            %from abouve array,remove edges that are repeating since it is sheetletwise
            % array 
            %Also remove edges that are shorter than 50 micro meters
            %Added length of each edge to this array in 7th column
            
            common_edge_array_ruc = [];
            cea_size = size(common_edge_array_stlt_wise);
            for ed = 1:cea_size(1)
                if common_edge_array_stlt_wise(ed,1)<common_edge_array_stlt_wise(ed,2)
                    v1_edge = [common_edge_array_stlt_wise(ed,3),common_edge_array_stlt_wise(ed,4)];
                    v2_edge = [common_edge_array_stlt_wise(ed,5),common_edge_array_stlt_wise(ed,6)];
                    
                    edge_length = norm(v1_edge-v2_edge);
                    % Remove any egdes shorter than given length (used 30
                    % micrometer in paper)
                    if edge_length > min_edgelength
                        common_edge_array_ruc = [common_edge_array_ruc;common_edge_array_stlt_wise(ed,:),edge_length];
                    else 
                        continue
                    end
                else
                    continue;
                end
            end

        end

        %-----------------------------------------------------------------%

        function [Possible_Edge_wise_branches] = find_potential_intersheetlet_branches(RUC_v_coord, common_edge_array_ruc, face1_point_array, coordinate_array, connectivity_array, branching_reg_thickness, linewise_point_array, f_rad)
            
            % find_potential_intersheetlet_branches
            %
            % Identifies candidate branches between neighbouring sheetlets.
            %
            % Candidate branches are generated only from cardiomyocyte
            % chains located within a prescribed distance from a common
            % sheetlet boundary.
            
            Possible_Edge_wise_branches = {};
    
            no_of_nodes = length(RUC_v_coord);
            node_branch_count = zeros(1,no_of_nodes);
            
            for com_edge = 1:length(common_edge_array_ruc(:,1))
                vtx1 = [0,common_edge_array_ruc(com_edge,3),common_edge_array_ruc(com_edge,4)];
                vtx2 = [0,common_edge_array_ruc(com_edge,5),common_edge_array_ruc(com_edge,6)];
                shtlt1 = common_edge_array_ruc(com_edge,1);
                shtlt2 = common_edge_array_ruc(com_edge,2);
                
                points_stlt1_f1 = face1_point_array{shtlt1};
                points_stlt2_f1 = face1_point_array{shtlt2};
                
                points_stlt1 = coordinate_array{shtlt1};
                points_stlt2 = coordinate_array{shtlt2};
                 
                %Sort coordinates of sheetlet 1 face1 points according to the distance from 
                %the common edge
                stlt1_sorted_p_list = RUC_microstructure_Toolbox.sort_pts_from_common_face(points_stlt1_f1,vtx1,vtx2);

                % Filter cardiomyocyte chains which are included within
                % 'branching_reg_thickness' from common edge
                bp_check = 1;
                outer_stlt_branches_bp_wise = {};
                stlt1_effective_bp_list = {}; %base points that can create non-overlapping branches
                
                eff_b = 1;
                for bp_in = 1:length(stlt1_sorted_p_list)
                    if stlt1_sorted_p_list{3,bp_in} <= branching_reg_thickness  %Considering basepoints only within 40 micrometer from the edge
                        stlt1_effective_bp_list{1,eff_b} = stlt1_sorted_p_list{1,bp_in};
                        stlt1_effective_bp_list{2,eff_b} = stlt1_sorted_p_list{2,bp_in};
                        eff_b=eff_b+1;
                    else
                        continue
                    end
                end

                for bp_no1 = 1: length(stlt1_effective_bp_list)

                    [Branches_for_BP1, node_branch_count] = RUC_microstructure_Toolbox.find_potential_branches_cellchain(bp_no1, ...
                    shtlt1, shtlt2, stlt1_effective_bp_list, points_stlt2_f1, linewise_point_array, coordinate_array, ...
                    points_stlt1, points_stlt2, connectivity_array, node_branch_count, f_rad);


                    outer_stlt_branches_bp_wise{end+1} = Branches_for_BP1;
                end

                %Create possible branches for all face 1 basepoint chains and extract
                %the first N_branch_edge number of basepoints with highest possible
                %branches in their chains andd append to the branch table
                
                branch_counts_edge = [];
                for bpn = 1:length(outer_stlt_branches_bp_wise)
                    branch_counts_edge(end+1) = length(outer_stlt_branches_bp_wise{bpn});
                end
                %Take the indexes of basepoints that can create branches for
                %considering common edge
                possible_bp = find(branch_counts_edge);
                
                effective_edge_branches = {};
                for add = 1:length(possible_bp)
                    add_id = possible_bp(add);
                    branches_bp = outer_stlt_branches_bp_wise{add_id};
                    for i3 = 1:length(branches_bp)
                        effective_edge_branches{end+1} = branches_bp{i3};
                    end
                end
                
                Possible_Edge_wise_branches{end+1} = effective_edge_branches;

            end
        end

        %-----------------------------------------------------------------%

        function [Branches_for_BP1, node_branch_count] = find_potential_branches_cellchain(bp_no1, ...
                shtlt1, shtlt2, stlt1_effective_bp_list, points_stlt2_f1, linewise_point_array, coordinate_array, ...
                points_stlt1, points_stlt2, connectivity_array, node_branch_count, f_rad)

            % find_potential_branches_cellchain
            %
            % Generates candidate branches between two neighboring
            % cardiomyocyte chains.
            %
            % Candidate branches must satisfy:
            %   - No duplicate connectivity
            %   - No overlap with existing cardiomyocytes
            %   - Node branch-count restrictions
            %
            % Valid branches are returned for later filtering.

            BP1_coord = stlt1_effective_bp_list{1,bp_no1};
            dis_bp1_bp2 = [];
            %Check the distances with all other base points in the adjacent
            %sheetlet
            for pid = 1:length(points_stlt2_f1)
                dis_bp1_bp2(end+1) = norm(points_stlt2_f1{pid}-BP1_coord);
            end
            [dis_bps,BP2_eff_list] = sort(dis_bp1_bp2); %get the closest sorted cardiomyocyte chains from selected BP1 chain
            
            Branches_for_BP1 = {};
            for itra = 1: 4
                I=BP2_eff_list(itra);
                BP2_coord = points_stlt2_f1{I}; %This is the closest base point in sheetlet 2(face1) to the selected basepoint in sheetlet 1
                
                linewise_BP1 = linewise_point_array{1,shtlt1}{1,stlt1_effective_bp_list{2,bp_no1}};
                linewise_BP2 = linewise_point_array{1,shtlt2}{1,I};
                
                %Find all possible branches that can be created between these two
                %chains adhering to overlapping and dirrection constraints
                
        
                %select a vertex from extended chain of base point 1, connect it to
                %the next face vertex of the adjacent sheetlet cardiomyocyte chain
                %and check whether if it is possible to have that branch
                
                max_vertex_id = min(length(linewise_BP1),length(linewise_BP2))-1;
                
                % Iterate through vertexes in cell chain
                for itr1 = 1:max_vertex_id
                
                    cb_v1_id = itr1 ;             %Cleavage branch vertex 1 is from the sheetlet 1
                    cb_v2_id = cb_v1_id + 1;      %vertex 2 from the adjacent myocyte chain is selected from next face 
                                                  %from V1 to maintain preferential
                                                  %dirrection of cariomyocytes
                    
                    cb_v1 = linewise_BP1{cb_v1_id};
                    cb_v2 = linewise_BP2{cb_v2_id};
                    
                    %Check whether the branch is already exist in the chain

                    cb_v1_graph_id = RUC_microstructure_Toolbox.find_coord_id(coordinate_array,cb_v1,shtlt1);
                    cb_v2_graph_id = RUC_microstructure_Toolbox.find_coord_id(coordinate_array,cb_v2,shtlt2);
                    branch_exist = 0;
                    for cell1 = 1:length(Branches_for_BP1)
                        if [cb_v1_graph_id,cb_v2_graph_id] == Branches_for_BP1{cell1} | [cb_v2_graph_id,cb_v1_graph_id] == Branches_for_BP1{cell1}
                            branch_exist = 1;
                            break;
                        end
                    end


                    if ~branch_exist
                        branches_overlap = 0;
                        %Should check whether this branch is overlaping with any other
                        %cardiomyocyte
        
                        shtlt1_myocytes = connectivity_array{shtlt1};
                        shtlt2_myocytes = connectivity_array{shtlt2};
                        
                        branch_clearance = RUC_microstructure_Toolbox.find_branch_clearance(cb_v1, cb_v2, shtlt1_myocytes, shtlt2_myocytes, points_stlt1, points_stlt2);
                        
                        if min(branch_clearance) < 2*f_rad %from sheetlet innitial
                           branches_overlap = 1;
                           continue;
                        end
                    end

                    %cHECK WHETHER THE SAME NODE HAS MORE THAN TWO BRANCHES
                    branch_count_per_node_exeed = 0;
                    if node_branch_count(1,cb_v1_graph_id)>0 || node_branch_count(1,cb_v2_graph_id)>0
                        branch_count_per_node_exeed = 1;
                    end
                    
                    % if the branch is not already available and branch is not
                    % overlapping with any other cardiomyocytes in the ruc, branch
                    % is added to a cardiomyocyte chainwise list
                    
                    if ~branches_overlap && ~branch_exist && ~branch_count_per_node_exeed
                        node_branch_count(cb_v1_graph_id) = node_branch_count(cb_v1_graph_id)+1;
                        node_branch_count(cb_v2_graph_id) = node_branch_count(cb_v2_graph_id)+1;
                        
                        Branches_for_BP1{end+1} = [cb_v1_graph_id,cb_v2_graph_id];
                    end

                end
            end

        end

        % ----------------------------------------------------------------%

        function branch_clearance = find_branch_clearance(cb_v1, cb_v2, shtlt1_myocytes, shtlt2_myocytes, points_stlt1, points_stlt2)
            
            % find_branch_clearance
            %
            % Computes minimum distances between a candidate branch and all
            % cardiomyocytes in two neighbouring sheetlets.
            %
            % Used to reject branches that would geometrically overlap
            % existing microstructure.

            branch_clearance = [];
            for cell2 = 1: length(shtlt1_myocytes)
                ch_v1_in = shtlt1_myocytes{cell2}(1);
                ch_v2_in = shtlt1_myocytes{cell2}(2);
                %neglect checking with cardiomyocytes include the same
                %node
                D_cells = Coordinate_Toolbox.DistBetween2Segment(cb_v1,cb_v2,points_stlt1{ch_v1_in},points_stlt1{ch_v2_in});
                if D_cells==0
                    branch_clearance(end+1) = 1000000;
                else
                    branch_clearance(end+1) = D_cells; 
                end
                
                
            end
            for cell2 = 1: length(shtlt2_myocytes)
                ch_v1_in = shtlt2_myocytes{cell2}(1);
                ch_v2_in = shtlt2_myocytes{cell2}(2);
                %neglect checking with cardiomyocytes include the same
                %node
                D_cells = Coordinate_Toolbox.DistBetween2Segment(cb_v1,cb_v2,points_stlt2{ch_v1_in},points_stlt2{ch_v2_in});
                if D_cells==0
                    branch_clearance(end+1) = 1000000;
                else
                    branch_clearance(end+1) = D_cells; 
                end
            end

        end


        %-----------------------------------------------------------------%

        function sorted_point_list = sort_pts_from_common_face(coord_list,line_vertex1,line_vertex2)
            
            % sort_pts_from_common_face
            %
            % Sorts cardiomyocyte base points according to their distance
            % from a common boundary between two neighbouring sheetlets.
            %
            % Used to identify chains most likely to generate
            % inter-sheetlet branches.

            point_dis_from_common_edge = [];
            
            % finding the seed points from sheetlet 1, within 15 micrometers from
            % the common face of two sheetlets
            
            for pts = 1:length(coord_list)
                point_f1_coord = coord_list{pts};
                point_f1_coord(1)=0;
                point_dis_from_common_edge(end+1) = Coordinate_Toolbox.min_dis_point_to_line(point_f1_coord,line_vertex1, line_vertex2);
            end
            
            [sorted_dis,sorted_id] = sort(point_dis_from_common_edge);
            
            sorted_point_list = {};
            for order = 1:length(sorted_id)
                sorted_point_list{1,order} = coord_list{sorted_id(order)};
                sorted_point_list{2,order} = sorted_id(order);
                sorted_point_list{3,order} = sorted_dis(order);
                
            end
            
        end

        %-----------------------------------------------------------------%

        function coordinate_id = find_coord_id(coordinate_array,p_cord,sheetlet_id)
        
            % find_coord_id
            %
            % Returns the global node index corresponding to a coordinate
            % within a specified sheetlet.
            %
            % Used when converting local sheetlet connectivity into the
            % global RUC graph representation.

            coordinate_id = 0;
        
            if sheetlet_id == 1
                coordinate_id = coordinate_id;
            else
                for stlt9 = 1:sheetlet_id-1
                    coordinate_id = coordinate_id + length(coordinate_array{1,stlt9});
                end
            end
            
            points_in_stlt = coordinate_array{sheetlet_id};
            for cd_no = 1:length(points_in_stlt)
                if p_cord == points_in_stlt{cd_no}
                    coordinate_id = coordinate_id + cd_no;
                else
                    continue
                end
            end

        end

        %-----------------------------------------------------------------%

        function Poly_area = Poly_area_calc(poly_connectivity,coordinates)
            
            % Poly_area_calc
            %
            % Computes the total cross-sectional area occupied by all
            % sheetlet polygons in the representative unit cell.
            %
            % Used for volume fraction calculations.

            Poly_area = 0;
            for poly = 1:length(poly_connectivity)
                polygon = poly_connectivity{poly};
                Vertex_coord_Y = [];
                Vertex_coord_Z = [];
                for vno = 1:length(polygon)
                    Vertex_coord_Y(end+1) = coordinates(polygon(vno),1);
                    Vertex_coord_Z(end+1) = coordinates(polygon(vno),2);
                end
                Add_area = polyarea(Vertex_coord_Y,Vertex_coord_Z);
                Poly_area = Poly_area + Add_area;
            end
        end

        %-----------------------------------------------------------------%

        function Final_branches = filter_final_branches(target_ICV, max_vol_frac, Total_volume, Vf_without_branching, Possible_ruc_branches, RUC_v_coord, f_rad)
            
            % filter_final_branches
            %
            % Randomly selects branches from the pool of candidate
            % inter-sheetlet branches until the target intracellular
            % volume fraction is approximately achieved.
            %
            % The selected branches form the final inter-sheetlet
            % connectivity network.

            if max_vol_frac >= target_ICV
                target_vol_frac = max_vol_frac;
            else
                target_vol_frac = target_ICV;
            end

            target_br_vol_frac = target_vol_frac - Vf_without_branching;
            target_br_vol = target_br_vol_frac * Total_volume;
            
            % Filter branches

            random_branch_id = randperm(length(Possible_ruc_branches));

            Final_branches = {};
            for br_id = 1:length(random_branch_id)
                Final_branches{end+1} = Possible_ruc_branches{random_branch_id(br_id)};
                cur_br_vol = Sheetlet_microstructure_Toolbox.Myo_vol_calc(Final_branches,RUC_v_coord,f_rad);
                if cur_br_vol>=target_br_vol
        
        
                    break
                else
                    continue
                end
            end


        end

        %-----------------------------------------------------------------%

        function shtletwise_cell_id_caps = shtltwise_cell_ids(RUC_con_table_sheetlet_wise)

            shtletwise_cell_id_caps = [];
            cap = 0;

            for stlt = 1:length(RUC_con_table_sheetlet_wise)

                cap = cap + length(RUC_con_table_sheetlet_wise{stlt});
                shtletwise_cell_id_caps(end+1) = cap;
            end

        end


        

    end

end