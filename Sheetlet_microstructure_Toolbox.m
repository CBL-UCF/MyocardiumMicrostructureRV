classdef Sheetlet_microstructure_Toolbox

    % Class - Sheetlet_microstructure_Toolbox
    % Static utility class for generating cardiomyocyte-like microstructure
    % inside individual sheetlets, including fiber chains, connectivity, and
    % possible intra-sheetlet branching.

    methods(Static)



        function [face1_pt_list, connected_p_list, con_table_branched_final, graph_v1_branched, graph_v2_branched, linewise_point_list] = generate_insheetlet_microstructure(L_shtlt, ...
                min_Y, min_Z, max_Y, max_Z, fiber_diameter, fiber_length_range, shtlt_ICV, ...
                Y_bp_vertex, Z_bp_vertex, max_fibredir_deviation, target_ICV_sheetlet, randomization_breakpoint, rng_count)

            % Function - generate_insheetlet_microstructure
            % Generates an intra-sheetlet cardiomyocyte network by placing
            % base points, extruding fibers along the fiber direction, creating
            % connectivity, identifying possible branches, and selecting final
            % branches to approach a target intra-cellular volume fraction.
            %
            % Inputs:
            %   L_shtlt                  - Sheetlet length in the fiber direction.
            %   min_Y, min_Z             - Minimum y and z bounds of the sheetlet.
            %   max_Y, max_Z             - Maximum y and z bounds of the sheetlet.
            %   fiber_diameter           - Cardiomyocyte/fiber diameter.
            %   fiber_length_range       - Minimum and maximum fiber segment length.
            %   shtlt_ICV                - Initial intra-cellular volume fraction.
            %   Y_bp_vertex, Z_bp_vertex - Polygon vertices defining the valid
            %                              base-point region in the sheetlet.
            %   max_fibredir_deviation   - Maximum lateral deviation from the
            %                              preferred fiber direction.
            %   target_ICV_sheetlet      - Target intra-cellular volume fraction
            %                              after adding branches.
            %   randomization_breakpoint - Maximum randomization attempts before
            %                              stopping point placement.
            %   rng_count                - Random seed value.
            %
            % Outputs:
            %   face1_pt_list            - Base points generated on the first face.
            %   connected_p_list         - Unique points included in the graph.
            %   con_table_branched_final - Final connectivity table after branching.
            %   graph_v1_branched        - First vertex index list for graph edges.
            %   graph_v2_branched        - Second vertex index list for graph edges.
            %   linewise_point_list      - Ordered point lists for each fiber chain.

            fiber_radius = fiber_diameter/2;
            % Approximate firber length projected in X direction
            appx_fiber_length = (fiber_length_range(1) + fiber_length_range(2)) / 2;
            
            % Volume fractions inside sheetlet
            Vf_myocyte = shtlt_ICV;
            Vf_ecm = 1 - shtlt_ICV;
            
            % Approximate number of fibre per chain along fibre dirrection
            Nf_l = round(L_shtlt/appx_fiber_length);
            
            % Approximate number of fibre per sheetlet cross section
            vor_area = polyarea(Y_bp_vertex,Z_bp_vertex);
            Nf_f =  round((0.95*(vor_area)/(pi*fiber_radius^2)),0); % Assume a larger value to iterate through

            % Generate myocyte basepoints in RUC face 1
            [face1_pt_list, Nf_f, area_frac] = Sheetlet_microstructure_Toolbox.define_fiberbasepoints_face1(Nf_f, min_Y, min_Z, max_Y, max_Z, Y_bp_vertex, Z_bp_vertex, fiber_radius, vor_area, randomization_breakpoint, rng_count);
    
            % Extrude points in fiber direction with slight deviations from
            % preferential direction

            rng_count = rng_count+5;
            facevise_p_list = Sheetlet_microstructure_Toolbox.extrude_points_in_fiberdir(Nf_l, face1_pt_list, max_fibredir_deviation, fiber_length_range, fiber_radius,  randomization_breakpoint, rng_count);

            % Construct the connectivity table for sheetlet
            [connected_p_list, linewise_table, con_table, graph_v1, graph_v2] = Sheetlet_microstructure_Toolbox.construct_connectivity_table(Nf_f, Nf_l, facevise_p_list, fiber_radius);
            
            % Approximate fiber volume fraction before branching 
            Vf_before_branch = Sheetlet_microstructure_Toolbox.fibreVolfrac(linewise_table, fiber_radius, L_shtlt, vor_area);

            % Creating branches inside sheetlet
            [linewise_point_list, possible_branch_table, branch_table_coord] = Sheetlet_microstructure_Toolbox.generate_inshtlt_branches(con_table, ...
                linewise_table, fiber_radius, face1_pt_list, connected_p_list);         
            

            % Add bracnhes selected from possible branches upto approximate
            % target volume fraction

            Total_volume = L_shtlt * vor_area;

            [con_table_branched_final, graph_v1_branched, graph_v2_branched] = Sheetlet_microstructure_Toolbox.select_final_branches(target_ICV_sheetlet, ...
                Vf_before_branch, fiber_radius, Total_volume, branch_table_coord, possible_branch_table, graph_v1, graph_v2, connected_p_list);

            % Monitor generate microstructure in 3D plot
            fig1 = Plotting_Toolbox.monitor_insheetlet_microstructure(1, connected_p_list, graph_v1_branched,graph_v2_branched);

        end

        % ----------------------------------------------------------------%

        function [face1_basepoints, Nf_f, area_frac] = define_fiberbasepoints_face1(Nf_f, ...
                min_Y, min_Z, max_Y, max_Z, Y_bp_vertex, Z_bp_vertex, fiber_radius, ...
                vor_area, rndm_br_point, rng_count)

            % Function - define_fiberbasepoints_face1
            % Randomly generates non-overlapping (circles) fiber base points on the first
            % face of a sheetlet while ensuring that all accepted points lie
            % inside the specified polygonal region.
            %
            % Inputs:
            %   Nf_f             - Target number of fiber base points.
            %   min_Y, min_Z     - Minimum y and z bounds of the sampling box.
            %   max_Y, max_Z     - Maximum y and z bounds of the sampling box.
            %   Y_bp_vertex      - y-coordinates of the valid polygon boundary.
            %   Z_bp_vertex      - z-coordinates of the valid polygon boundary.
            %   fiber_radius     - Fiber radius used for overlap checking.
            %   vor_area         - Area of the sheetlet polygon.
            %   rndm_br_point    - Maximum number of randomization attempts.
            %   rng_count        - Random seed value.
            %
            % Outputs:
            %   face1_basepoints - Accepted base-point coordinates on face 1.
            %   Nf_f             - Final number of accepted base points.
            %   area_frac        - Area fraction occupied by accepted fibers.


            face1_basepoints = {};

            W_shtlt = max_Y-min_Y;  % bounding box width in micrometers (Along Sheetlet dirrection)
            T_shtlt = max_Z-min_Z;  % bounding box thickness in micrometers (Along normal dirrection)

            % Function to calculate euclidean distance between two points
            Dis = @(a1,a2,b1,b2,c1,c2) sqrt(((a1-a2)^2+(b1-b2)^2+(c1-c2)^2));

            p_face=0;
            br_point=0;

            rng(rng_count) % to innitialize the random number generating global stream

            while p_face<Nf_f

                P_x = 0;
                P_y = min_Y + rand(1,1)*W_shtlt;
                P_z = min_Z + rand(1,1)*T_shtlt;

                % make sure that the points are inside the Thiessen polygon 
                % defined for fiber basepoints                    
                inout_vor = inpolygon(P_y,P_z,Y_bp_vertex,Z_bp_vertex);
                if ~inout_vor
                    continue
                end
                                
                br_point=br_point+1; 

                % make sure distance between basepoints are more than cell diameter
                D_list=[];
                for c = 1:length(face1_basepoints)
                    D = Dis(P_x,face1_basepoints{c}(1),P_y,face1_basepoints{c}(2),P_z,face1_basepoints{c}(3));
                    D_list(end+1) = D;
                end

                %Introducing a breakpoint if number of cardiomyocytes across the space
                %cannot reach 50% of the area
            
                if br_point==rndm_br_point
                    fprintf('Randomization iterations reached %i and generated points = %i / %i \n', rndm_br_point, length(face1_basepoints),Nf_f)
                    break
                end

                %Make sure cardiomycytes in the region do not overlap each other  
                if min(D_list)<2*fiber_radius
                    continue    
                else 
                    face1_basepoints{end+1} = [P_x,P_y,P_z];
                    p_face=p_face+1;
                end

            end

            Nf_f = length(face1_basepoints);
            area_frac = (Nf_f*pi*fiber_radius^2)/vor_area;         

        end

        % ----------------------------------------------------------------%

        function facevise_p_list = extrude_points_in_fiberdir(Nf_l, face1_pt_list, ...
                max_fibredir_deviation, fiber_length_range, fiber_radius,  randomization_breakpoint, rng_count)

            % Function - extrude_points_in_fiberdir
            % Extrudes each base point along the fiber direction to generate
            % longitudinal fiber chains. New points are accepted only when they
            % satisfy the minimum spacing requirement relative to existing points.
            %
            % Inputs:
            %   Nf_l                     - Number of fiber segments per chain.
            %   face1_pt_list            - Base points on the first sheetlet face.
            %   max_fibredir_deviation   - Maximum lateral deviation in y and z.
            %   fiber_length_range       - Minimum and maximum segment length.
            %   fiber_radius             - Fiber radius for spacing checks.
            %   randomization_breakpoint - Maximum number of placement attempts.
            %   rng_count                - Random seed value.
            %
            % Output:
            %   facevise_p_list          - Cell array of points arranged by face
            %                              along the fiber direction.

            %Consider each point in the first face of the sheetlet and extrude those
            %points Nf_l times, checking minimum distance of fibre diameter with any other point
            %in the whole point list

            Point_list = face1_pt_list;

            appx_fiber_length = (fiber_length_range(1) + fiber_length_range(2)) / 2;
            fiber_l_deviation = fiber_length_range(2) - fiber_length_range(1);

            facevise_p_list = cell(1, Nf_l + 1); % create empty cells in a cell array for each plane along x dirrection
            facevise_p_list{1} = face1_pt_list;

            rng(rng_count)

            Dis = @(a1,a2,b1,b2,c1,c2) sqrt(((a1-a2)^2+(b1-b2)^2+(c1-c2)^2));

            for pt=1:length(face1_pt_list) % Iterate through points in first face and extrude
                Base_point = face1_pt_list{pt};
                
                ln=1;
                                
                br_point=0;
                while ln <= Nf_l % Iterate though number of cardiomyocyres per chain
                    Px_new = Base_point(1)+ ln*appx_fiber_length -(fiber_l_deviation/2) + fiber_l_deviation*rand(1,1);
                    Py_new = Base_point(2)-(max_fibredir_deviation/2)+max_fibredir_deviation*rand(1,1);     
                    Pz_new = Base_point(3)-(max_fibredir_deviation/2)+max_fibredir_deviation*rand(1,1);
                    
                    D_list=[];
                    for c = 1:length(Point_list)
                        D = Dis(Px_new,Point_list{c}(1),Py_new,Point_list{c}(2),Pz_new,Point_list{c}(3));
                        D_list(end+1) = D;
                    end
                    
                    br_point = br_point+1;
                    
                    if br_point==randomization_breakpoint
                        fprintf('Randomization iterations reached %i and extruded points = %i / %i for %i base point \n', randomization_breakpoint, ln,Nf_l, pt)
                        for it = 1:Nf_l-ln+1
                            Point_list{end+1} = [NaN, NaN, NaN];
                            facevise_p_list{ln+1}{end+1} = [NaN, NaN, NaN]; %Add points to the facevise point list
                            %p_face = p_face + 1;
                            ln=ln+1;
                        end
                        break
                    end
                    
                   if min(D_list) < 2*fiber_radius
                        continue    
                   else 
                        Point_list{end+1} = [Px_new,Py_new,Pz_new];
                        facevise_p_list{ln+1}{end+1} = [Px_new,Py_new,Pz_new]; %Add points to the facevise point list
                        %p_face = p_face + 1;
                        ln = ln + 1;
                   end
                end
            end
        end

        % ----------------------------------------------------------------%

        function [connected_p_list, linewise_table, con_table, graph_v1, graph_v2] = construct_connectivity_table(Nf_f, Nf_l, facevise_p_list, fiber_radius)
            
            % Function - construct_connectivity_table
            % Constructs fiber-chain connectivity from face-wise point lists.
            % The function connects consecutive points along each fiber chain,
            % checks spacing from existing line segments, and stores unique graph
            % vertices and edge connectivity.
            %
            % Inputs:
            %   Nf_f            - Number of fiber chains.
            %   Nf_l            - Number of fiber segments per chain.
            %   facevise_p_list - Cell array of points arranged by face.
            %   fiber_radius    - Fiber radius used for line-spacing checks.
            %
            % Outputs:
            %   connected_p_list - Unique point coordinates included in the graph.
            %   linewise_table   - Connectivity stored separately for each chain.
            %   con_table        - Global connectivity table.
            %   graph_v1         - First vertex index list for graph edges.
            %   graph_v2         - Second vertex index list for graph edges.

            Dis = @(a1,a2,b1,b2,c1,c2) sqrt(((a1-a2)^2+(b1-b2)^2+(c1-c2)^2));

            con_table = {};
            linewise_table = cell(1, Nf_f);
            in_v = 0;

            %innitiate connectivity of vertex lists of the graph
            graph_v1=[];
            graph_v2=[];

            %innitiate a list of points that actually got a connection after all checks
            connected_p_list = {};

            break_extending_line = 0; 
            min_line_dis=[];


            for bp=1:Nf_f           %id of the basepoint that is extending
                for v = 1:Nf_l      %id of the cell number in lingitudinal dirrection(x) coreespoinding correspoinding to basepoint
                    v1=facevise_p_list{v}{bp};
                    v2=facevise_p_list{v+1}{bp};
                    
                    %Checking minimum distance of the line connecting v1 and v2 with 
                    %all other lines already cerated
                    
                    line_dis_list=[1000];
                    
                    for bp_c = 1:bp-1

                        for ln_c = 1:length(linewise_table{bp_c})
                            v1_c=linewise_table{bp_c}{ln_c}{1};
                            v2_c=linewise_table{bp_c}{ln_c}{2};
                            
                            % ddd = min_dis_2lines(v1,v2,v1_c,v2_c); %Calling function to calculate distance
                            ddd = Coordinate_Toolbox.DistBetween2Segment(v1,v2,v1_c,v2_c);
                            line_dis_list(end+1) = ddd;
                            break_extending_line = 0;
                            if ddd < 2*fiber_radius
                                v2=v2_c;    %Since always v2_c is the one have higher 
                                            %x coordinate
                                            break_extending_line = 1;
                                break;   %%% Check weather you need to check the 
                                         %%%minimum distance for this line as well
                            end
                        end

                        if break_extending_line==1
                                break;
                        end
                    end
                        
                    aaaa=min(line_dis_list);
                    min_line_dis(end+1)= aaaa;
                   
                    linewise_table{bp}{end+1}={v1,v2};  %cell which have the all connected point pairs
                    
                    v1_connected=false;                 %Check whether vertex 1 is already in the graph
                    for in_v = 1:length(connected_p_list)  % If vertex 1 is already in the grapg, get the index to create connectivity table
                        if connected_p_list{in_v}==v1
                            v1_in = in_v;
                            v1_connected=true;
                            break
                        end
                    end
                    if not(v1_connected)                %if vertex 1 is not in the add it to both corrdinates list and connectivity list
                        connected_p_list{end+1}=v1;
                        v1_in=length(connected_p_list);
                        
                    end
                    
                    v2_connected=false;                %Check whether vertex 2 is already in the graph
                    for in_v = 1:length(connected_p_list)   % If vertex 2 is already in the grapg, get the index to create connectivity table
                        if connected_p_list{in_v}==v2
                            v2_in=in_v;
                            v2_connected=true;
                            break
                        end
                    end
                    if not(v2_connected)                 %if vertex 2 is not in the add it to both corrdinates list and connectivity list
                        connected_p_list{end+1}=v2;
                        v2_in=length(connected_p_list);
                    end
                    
                    graph_v1(end+1)= v1_in;
                    graph_v2(end+1)= v2_in;
                    
                    con_table{end+1}=[v1_in,v2_in];
                    
                    if break_extending_line == 1
                        break;
                    end
                   
                end
                
            end
        end
        % ----------------------------------------------------------------%

        function [linewise_point_list, possible_branch_table, branch_table_coord] = generate_inshtlt_branches(con_table, ...
                linewise_table, fiber_radius, face1_pt_list, connected_p_list)

            % Function - generate_inshtlt_branches
            % Identifies possible intra-sheetlet branches between neighboring
            % fiber chains. Candidate branches are accepted only if they do not
            % already exist, do not overlap existing fibers, and do not exceed
            % the permitted number of branches per node.
            %
            % Inputs:
            %   con_table        - Existing fiber connectivity table.
            %   linewise_table   - Fiber-chain connectivity organized by chain.
            %   fiber_radius     - Fiber radius used for overlap checking.
            %   face1_pt_list    - Base points for all fiber chains.
            %   connected_p_list - Unique graph point coordinates.
            %
            % Outputs:
            %   linewise_point_list   - Ordered point list for each fiber chain.
            %   possible_branch_table - Candidate branch connectivity pairs.
            %   branch_table_coord    - Candidate branch endpoint coordinates.

            linewise_point_list = Sheetlet_microstructure_Toolbox.create_linewise_pointlist(linewise_table);
            

            con_table_branched = con_table;

            possible_branch_table = {};
            branch_table_coord = {};  

            

            branch_bp_count = length(face1_pt_list);  % No of cell chains

            no_of_nodes = length(connected_p_list);
            node_branch_count = zeros(1,no_of_nodes);

            for Bp = 1:branch_bp_count
                Bp1_coord = face1_pt_list{1,Bp};
                 
                                 
                dis_bp1_bp2 = [];
                for pid = 1:length(face1_pt_list)
                    dis = norm([0,face1_pt_list{pid}(2),face1_pt_list{pid}(3)]-[0,Bp1_coord(2),Bp1_coord(3)]);
                    if dis == 0
                       dis_bp1_bp2(end+1) = 100000; % Define maximum distance for own point to avoid mistaken as minimum distance
                    else
                       dis_bp1_bp2(end+1) = dis;
                    end
                end

                [dis_bps,Bp2_eff_list] = sort(dis_bp1_bp2); %get the closest sorted cardiomyocyte chains from selected BP1 chain
                 
                Branches_Bp1 = {};
                 
                Eff_bp2_no = min(10,length(face1_pt_list)); %No of chains considered as potential chains for a selected BP1 chain to create branches
                
                for itra = 1: Eff_bp2_no

                    I = Bp2_eff_list(itra);
                     %[M,I] = min(dis_bp1_bp2);
                    BP2_coord = face1_pt_list{I}; % This is the I'th closest base point in sheetlet 2(face1) to the selected basepoint in sheetlet 1
            
                    linewise_Bp1 = linewise_point_list{Bp};
                    linewise_Bp2 = linewise_point_list{I};
                    
                    %Find all possible branches that can be created between these two
                    %chains adhering to overlapping and dirrection constraints
            
            
                    %select a vertex from extended chain of base point 1, connect it to
                    %the next face vertex of the adjacent sheetlet cardiomyocyte chain
                    %and check whether if it is possible to have that branch
            
                    max_vertex_id = min(length(linewise_Bp1),length(linewise_Bp2))-1;
            
                    %create iterations for number of branches in cardiomyocyte chain
                    %for that specific basepoint
            
            
                    for itr1 = 1:max_vertex_id
                        cb_v1_id = itr1 ;             %Cleavage branch vertex 1 is from the sheetlet 1
                        cb_v2_id = cb_v1_id + 1;      %vertex 2 from the adjacent myocyte chain is selected from next face 
                                                        %wrt V1 to maintain preferential
                                                        %dirrection of cariomyocytes
            
                        cb_v1 = linewise_Bp1{cb_v1_id};
                        cb_v2 = linewise_Bp2{cb_v2_id};
            
                        %Check whether the branch is already exist in the chain
                        
                        cb_v1_graph_id = Sheetlet_microstructure_Toolbox.find_coord_id_shtlt(connected_p_list,cb_v1);
                        cb_v2_graph_id = Sheetlet_microstructure_Toolbox.find_coord_id_shtlt(connected_p_list,cb_v2);
                        
                        branch_exist = 0;
                        for cell1 = 1:length(con_table_branched)
                            if [cb_v1_graph_id,cb_v2_graph_id] == con_table_branched{cell1} | [cb_v2_graph_id,cb_v1_graph_id] == con_table_branched{cell1} | cb_v1_graph_id==cb_v2_graph_id
                                branch_exist = 1;
                                break;
                            end
                        end
                          
                        if ~branch_exist
                            branches_overlap = 0;
                                %Should check whether this branch is overlaping with any other
                                %cardiomyocyte
            
                            cell_spcaes = [];
                            for cell2 = 1: length(con_table_branched)
                                 ch_v1_in = con_table_branched{cell2}(1);
                                 ch_v2_in = con_table_branched{cell2}(2);
                                 %neglect checking with cardiomyocytes include the same
                                 %node
                                 
                                 
                                 D_cells = Coordinate_Toolbox.DistBetween2Segment(cb_v1,cb_v2,connected_p_list{ch_v1_in},connected_p_list{ch_v2_in});
                                 %if D_cells < 1e-10
                                 if cb_v1_graph_id == ch_v1_in || cb_v2_graph_id == ch_v1_in || cb_v1_graph_id == ch_v2_in || cb_v2_graph_id == ch_v2_in
                                     cell_spcaes(end+1) = 100000;
                                 else
                                     cell_spcaes(end+1) = D_cells;
                                 end
            
                            end
                            
                            if min(cell_spcaes)<2*fiber_radius
                                 branches_overlap = 1;
                                 continue;
                            end
                        
                        end
                        
                        %Check whether the same node has more than two
                        %branches
                        branch_count_per_node_exeed = 0;
                        if node_branch_count(1,cb_v1_graph_id)>1 || node_branch_count(1,cb_v2_graph_id)>1
                            branch_count_per_node_exeed = 1;
                        end
                            
                            
                            
                            
                        if  ~branch_exist && ~branches_overlap && ~branch_count_per_node_exeed
                            
                            node_branch_count(cb_v1_graph_id) = node_branch_count(cb_v1_graph_id)+1;
                            node_branch_count(cb_v2_graph_id) = node_branch_count(cb_v2_graph_id)+1;
                            
                            Branches_Bp1{end+1} = [cb_v1_graph_id,cb_v2_graph_id];
                            possible_branch_table{end+1} = [cb_v1_graph_id,cb_v2_graph_id];
                            branch_table_coord{end+1} = {cb_v1,cb_v2};
                            con_table_branched{end+1} = [cb_v1_graph_id,cb_v2_graph_id];

                        end
                    end
               end
           
            end

        end

        % ----------------------------------------------------------------%

        function linewise_point_list = create_linewise_pointlist(linewise_table)

            % Function - create_linewise_pointlist
            % Converts a line-wise connectivity table into ordered point lists
            % for each fiber chain.
            %
            % Input:
            %   linewise_table      - Connectivity table organized by fiber chain.
            %
            % Output:
            %   linewise_point_list - Ordered point coordinates for each chain.

            linewise_point_list = {};
            for lpi=1:length(linewise_table)
                pointss = {};
                for lpi1 = 1:length(linewise_table{lpi})
                    pointss{end+1} = linewise_table{lpi}{lpi1}{1};
                    if lpi1 == length(linewise_table{lpi})
                        pointss{end+1} = linewise_table{lpi}{lpi1}{2};
                    end
                end
                linewise_point_list{end+1} = pointss;
            end      

        end

        % ----------------------------------------------------------------%

        function Vf = fibreVolfrac(con_list,fibre_r,RUC_l,vor_area)  %Calculating total fibre volume fraction in the RUC acoording to connectivity table
            
            % Function - fibreVolfrac
            % Calculates the total fiber volume fraction from a chain-wise
            % connectivity list by summing segment lengths and assuming each
            % segment is a cylinder.
            %
            % Inputs:
            %   con_list - Connectivity list organized by fiber chain.
            %   fibre_r  - Fiber radius.
            %   RUC_l    - Representative unit cell length.
            %   vor_area - Cross-sectional polygon area.
            %
            % Output:
            %   Vf       - Fiber volume fraction, expressed as a percentage.

            total_fibre_length = 0;
            
            
            for bp_vol = 1:length(con_list)
                for line = 1:length(con_list{bp_vol})
                    total_fibre_length = total_fibre_length + norm(con_list{bp_vol}{line}{2}-con_list{bp_vol}{line}{1});
                end
            end
            
            total_fibre_volume = total_fibre_length*pi*fibre_r^2;
            RUC_total_vol = RUC_l*vor_area;
            Vf = (total_fibre_volume/RUC_total_vol)*100;
            
        end 

        % ----------------------------------------------------------------%

        function coordinate_id = find_coord_id_shtlt(coordinate_array,p_cord)

            % Function - find_coord_id_shtlt
            % Finds the index of a specified coordinate within a cell array of
            % stored coordinates.
            %
            % Inputs:
            %   coordinate_array - Cell array containing coordinate vectors.
            %   p_cord           - Coordinate vector to locate.
            %
            % Output:
            %   coordinate_id    - Index of the matching coordinate.

            for cd_no = 1:length(coordinate_array)
                if p_cord == coordinate_array{cd_no}
                    coordinate_id = cd_no;
                else
                    continue
                end
            end
        
        end

        % ----------------------------------------------------------------%

        function [con_table_branched_final, graph_v1_branched, graph_v2_branched] = select_final_branches(Required_V_frac, ...
                Vf_before_branch, fiber_radius, Total_volume, branch_table_coord, possible_branch_table, graph_v1, graph_v2, connected_p_list)
            
            % Function - select_final_branches
            % Randomly selects a subset of possible branches to approach the
            % required fiber volume fraction after branching. The selected
            % branches are appended to the existing graph connectivity.
            %
            % Inputs:
            %   Required_V_frac      - Required final volume fraction.
            %   Vf_before_branch     - Fiber volume fraction before branching.
            %   Total_volume         - Total volume of the sheetlet region.
            %   branch_table_coord   - Endpoint coordinates for possible branches.
            %   possible_branch_table- Connectivity pairs for possible branches.
            %   graph_v1             - Existing first vertex index list.
            %   graph_v2             - Existing second vertex index list.
            %
            % Outputs:
            %   con_table_branched_final - Final connectivity table including
            %                              selected branches.
            %   graph_v1_branched        - Updated first vertex index list.
            %   graph_v2_branched        - Updated second vertex index list.

            % Approximate fiber volume fraction after branching

            Max_branch_vol = 0;
            for lines = 1:length(branch_table_coord)
                Max_branch_vol = Max_branch_vol + pi * fiber_radius^2 *norm(branch_table_coord{lines}{2}-branch_table_coord{lines}{1}); % Assumed cells as cylinders
            end
            
            Max_Vf_branched= Vf_before_branch + 100 * Max_branch_vol/Total_volume;

            % Select final branches considering input approximate volume
            % fraction
            random_branch_id = randperm(length(possible_branch_table));
            
            graph_v1_branched = graph_v1;
            graph_v2_branched = graph_v2;

            if Max_Vf_branched >= Required_V_frac
                target_vol_frac = Max_Vf_branched;
            else
                target_vol_frac = Required_V_frac;
            end

            target_br_vol_frac = target_vol_frac - Vf_before_branch;
            target_br_vol = target_br_vol_frac*Total_volume;
            
            % Pick final branches to reach given volume fracktion
            Final_branches = {};
            for br_id = 1:length(random_branch_id)
                Final_branches{end+1} = possible_branch_table{random_branch_id(br_id)};
                cur_br_vol = Sheetlet_microstructure_Toolbox.Myo_vol_calc(Final_branches,connected_p_list,fiber_radius);
                if cur_br_vol>=target_br_vol
                    break
                else
                    continue
                end
            end

            for branch_no1 = 1:length(Final_branches)
                graph_v1_branched(end+1) = Final_branches{branch_no1}(1);
                graph_v2_branched(end+1) = Final_branches{branch_no1}(2);
            end
             
            con_table_branched_final = {};
            for cell_no = 1:length(graph_v1_branched)
                con_table_branched_final{end+1} = [graph_v1_branched(cell_no),graph_v2_branched(cell_no)];
            end
        end

        % ----------------------------------------------------------------%

        function cell_vol = Myo_vol_calc(cell_array,cell_coord,cell_radius)
            % Myo_vol_calc
            % Calculates the total myocyte/cell volume represented by a connectivity
            % array. Each connected cell segment is approximated as a cylinder whose
            % length is the distance between two connected node coordinates.
            %
            % Inputs:
            %   cell_array  - Cell array containing connectivity pairs. Each entry
            %                 stores the indices of two connected cell coordinates.
            %   cell_coord  - Cell array containing coordinate vectors for each node.
            %   cell_radius - Radius of the cylindrical cell segment.
            %
            % Output:
            %   cell_vol    - Total volume of all connected cylindrical cell segments.


            cell_vol = 0;
            for cell = 1:length(cell_array)
                cell_v1 = cell_array{cell}(1);
                cell_v2 = cell_array{cell}(2);
                
                cell_len = norm(cell_coord{cell_v1}-cell_coord{cell_v2});
                cell_area = pi*(cell_radius^2);
                cell_vol_add = cell_len*cell_area;
                cell_vol = cell_vol + cell_vol_add;
            end
        end

    end
end

            