%% Inputs

% Dimensions of RUC
% Approximate sheetlet width(250 microm) and thickness(100 microm) used to
% create voronoi patterns - Sheetlet_voronoi.m

%b_sheetlet_arrangement
[Shrinked_after_inf, Shrinked_after_inf_bp, Shrinked_after_inf_mid, non_shrink_vor_ruc, Nn, Ns] = Sheetlet_arrangement_Toolbox.generate_sheetlet_arrangement(RUC_w, RUC_t, shtlt_w, shtlt_t, shtlt_shrink_ratios);

fig1 = Plotting_Toolbox.monitor_sheetlet_boundaries_in3D(1, 'Monitor Microstructure Construction', Shrinked_after_inf, RUC_l);

sheetlet_face1 = Shrinked_after_inf;
sheetlet_face1_bp = Shrinked_after_inf_bp;

stlt_count = length(sheetlet_face1);
rng_count = 0;



%%
% Record coordinates of each vertex in the RUC and connectivity
    
coordinate_array = {};
connectivity_array = {};
connectivity_V1_array = {};
connectivity_V2_array = {};
face1_point_array = {};
linewise_point_array = {};
area_frac_list = [];

for stlt = 1: stlt_count
    
    stlt_vertex = sheetlet_face1{1,stlt};
    min_Y = min(stlt_vertex(:,1));
    min_Z = min(stlt_vertex(:,2));
    max_Y = max(stlt_vertex(:,1));
    max_Z = max(stlt_vertex(:,2));

    SLT_L = RUC_l;        %length in micrometers(Along fibre dirrection)
    SLT_W = max_Y-min_Y;  %Width in micrometers (Along Sheetlet dirrection)
    SLT_T = max_Z-min_Z;  %Thickness in micrometers (Along normal dirrection)
    
    %Srink voronoi in 1:9 ratio to place base points of the extending
    %cardiomyocyte chains inside the voronoi
    
    stlt_bp_boundary = sheetlet_face1_bp{1,stlt};
    stlt_bp_boundary(end+1,:) = stlt_bp_boundary(1,:); % Add first point to close polygon
    Y_bp_vertex = transpose(stlt_bp_boundary(:,1)); % Y coordinates of basepoint polygons
    Z_bp_vertex = transpose(stlt_bp_boundary(:,2)); % Z coordinates of basepoint polygons
       
    
    c_sheetlet_innitial
    rng_count = rng_count+1;
end

%% Getting full cardiomyocyte network into a single graph definition

RUC_v_coord = {};

for stlt1 = 1:length(coordinate_array)
    coords_stlt = coordinate_array{1,stlt1};
    for coorda = 1:length(coords_stlt)
        RUC_v_coord{end+1} = coords_stlt{coorda};
    end
end
        
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
        
  %% Graphical representation of connectivity network in RUC

figure(1)

hold on;
XX_ruc=[];
YY_ruc=[];
ZZ_ruc=[];

for pt_g_ruc=1:length(RUC_v_coord)
    XX_ruc(end+1)=RUC_v_coord{pt_g_ruc}(1);
    YY_ruc(end+1)=RUC_v_coord{pt_g_ruc}(2);
    ZZ_ruc(end+1)=RUC_v_coord{pt_g_ruc}(3);
end

con_graph_ruc=graph(RUC_con_v1,RUC_con_v2);
pp_ruc=plot(con_graph_ruc);
pp_ruc.XData=XX_ruc;
pp_ruc.YData=YY_ruc;
pp_ruc.ZData=ZZ_ruc;      
  


 %% Identifying neghbouring sheetlets
 
 if Mesoscale == 1
     N_stlt = Ns-2;
     N_normal = Nn-2;
     
     Neghbour_stlt_ids = {};
      
     % Indexes of neighbouring sheetlets  for each sheetlet in RUC 
     for stlt3 = 1:length(face1_point_array)
         if stlt3 == 1
             Neghbour_stlt = [stlt3+1, stlt3+N_normal];
         elseif stlt3 == N_normal
             Neghbour_stlt = [stlt3-1, stlt3+N_normal];
         elseif stlt3 == N_normal*N_stlt
             Neghbour_stlt = [stlt3-1, stlt3-N_normal];
         elseif stlt3 == N_normal*(N_stlt-1)+1
             Neghbour_stlt = [stlt3+1, stlt3-N_normal];
         elseif stlt3 < N_normal
             Neghbour_stlt = [stlt3-1, stlt3+1, stlt3+N_normal];
         elseif stlt3 > N_normal*(N_stlt-1)+1
             Neghbour_stlt = [stlt3-1, stlt3+1, stlt3-N_normal];
         elseif rem(stlt3,N_normal) == 0
             Neghbour_stlt = [stlt3-N_normal, stlt3-1, stlt3+N_normal];
         elseif rem(stlt3,N_normal) == 1
             Neghbour_stlt = [stlt3-N_normal, stlt3+1, stlt3+N_normal];
         else
             Neghbour_stlt = [stlt3-N_normal, stlt3+1, stlt3+N_normal,stlt3-1];
         end
         Neghbour_stlt_ids{end+1} =Neghbour_stlt; 
     end
    
    
    DT_ruc_before = delaunayTriangulation(seed_array);
    % DT_ruc = delaunayTriangulation(seeds_ruc);
    
    
    seed_con = DT_ruc_before.ConnectivityList;
    
    all_seed_ids = (1:1:length(seed_points));
    
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
    
    Neighbour_stlt_array = {};
    
    for stlt4 = 1:length(RUC_seed_ids)
        neigbours =[];
        sheetlet_no = RUC_seed_ids(stlt4);
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
    
    Neghbour_stlt_array_ruc = {}; %array of sheetlet ids only in RUC
    
    %From connected points to each vertexes, select only points that are in RUC
    %and neglect others
    for stlt5 = 1:length(Neighbour_stlt_array)
        ruc_neighbour = intersect(Neighbour_stlt_array{stlt5},RUC_seed_ids);
        Neghbour_stlt_array_ruc{end+1} = ruc_neighbour;
    end
    
    %Above array indexes are correseponding to seed point id before removing
    %outer layer of sheetlets and voronois containing infinity as a vertex.
    %Now, we redefine the same array using sheetlet ids only inside the RUC
    
    Neghbour_stlt__rucid = {};
    for stlt6 = 1:length(Neghbour_stlt_array_ruc)
        ruc_neighbour1 = [];
        for vv = 1:length(Neghbour_stlt_array_ruc{stlt6})
            ruc_neighbour1(end+1) = find(RUC_seed_ids==Neghbour_stlt_array_ruc{stlt6}(vv));
        end
        Neghbour_stlt__rucid{end+1} = ruc_neighbour1;
    end
    
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
    % shtlt1 | shtlt2 | commonedge1(Y) | commonedge1(Z) | commonedge2(Y) | commonedge2(Z) |
    
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
            % Remove any egdes shorter than 50 micrometers
            if edge_length > 30
                common_edge_array_ruc = [common_edge_array_ruc;common_edge_array_stlt_wise(ed,:),edge_length];
            else 
                continue
            end
        else
            continue;
        end
    end
 end

%% Create between sheetlet cardiomyocyte branches

% extract closest face 1 points for each common edge
%Assuming face 1 is in x=0 plane, finding points near to cboundaries

if Mesoscale == 1

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
        
        stlt1_sorted_p_list = sorted_pts_from_common_face(points_stlt1_f1,vtx1,vtx2);
        
        %Select the number of branches per common boundary based on edge length
        %Selected #branches is 1/50 of the edge length (myocyte D = 15)
        %N_branch_edge = round(common_edge_array_ruc(com_edge,7)/50);
        
        
        
        %Iterate the simmulation until create branches in N_branch # base
        %points
        bp_check = 1;
        
        outer_stlt_branches_bp_wise = {};
        
        stlt1_effective_bp_list = {}; %base points that can create non-overlapping branches
        
        eff_b = 1;
        for bp_in = 1:length(stlt1_sorted_p_list)
            if stlt1_sorted_p_list{3,bp_in}<=40  %Considering basepoints only within 40 micrometer from the edge
                stlt1_effective_bp_list{1,eff_b} = stlt1_sorted_p_list{1,bp_in};
                stlt1_effective_bp_list{2,eff_b} = stlt1_sorted_p_list{2,bp_in};
                eff_b=eff_b+1;
            else
                continue
            end
        end
        
        
        
        
        for bp_no1 = 1: length(stlt1_effective_bp_list)
            BP1_coord = stlt1_effective_bp_list{1,bp_no1}; %From sheetlet1_face 1 points
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
                
                %create iterations for number of branches in cardiomyocyte chain
                %for that specific basepoint
                
                
                
                
                for itr1 = 1:max_vertex_id
                    
                    cb_v1_id = itr1 ;             %Cleavage branch vertex 1 is from the sheetlet 1
                    cb_v2_id = cb_v1_id + 1;      %vertex 2 from the adjacent myocyte chain is selected from next face 
                                                  %from V1 to maintain preferential
                                                  %dirrection of cariomyocytes
                    
                    cb_v1 = linewise_BP1{cb_v1_id};
                    cb_v2 = linewise_BP2{cb_v2_id};
                    
                    %Check whether the branch is already exist in the chain
                    
                                
                    cb_v1_graph_id = find_coord_id(coordinate_array,cb_v1,shtlt1);
                    cb_v2_graph_id = find_coord_id(coordinate_array,cb_v2,shtlt2);
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
                        
                        branch_clearance = [];
                        for cell2 = 1: length(shtlt1_myocytes)
                            ch_v1_in = shtlt1_myocytes{cell2}(1);
                            ch_v2_in = shtlt1_myocytes{cell2}(2);
                            %neglect checking with cardiomyocytes include the same
                            %node
                            D_cells = min_dis_2lines(cb_v1,cb_v2,points_stlt1{ch_v1_in},points_stlt1{ch_v2_in});
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
                            D_cells = min_dis_2lines(cb_v1,cb_v2,points_stlt2{ch_v1_in},points_stlt2{ch_v2_in});
                            if D_cells==0
                                branch_clearance(end+1) = 1000000;
                            else
                                branch_clearance(end+1) = D_cells; 
                            end
                        end
                        
                        if min(branch_clearance)<2*f_rad %from sheetlet innitial
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

%%
if Mesoscale == 1
    Possible_ruc_branches = {};
    for edge_no = 1:length(Possible_Edge_wise_branches)
        b_s = Possible_Edge_wise_branches{edge_no};
        for br = 1:length(b_s)
            Possible_ruc_branches{end+1} = b_s{br};
        end
    end
end
    
%% Available volume fraction calculations
Total_area = Poly_area_calc(non_shrink_vor_ruc,vor_vertex);
Total_volume = Total_area*RUC_l;

%f_rad imported from sheetlet innitial
Sheetlet_myocyte_vol = Myo_vol_calc(RUC_con_table,RUC_v_coord,f_rad);

Vf_without_branching = Sheetlet_myocyte_vol/Total_volume;

if Mesoscale == 1
    Possible_branch_vol = Myo_vol_calc(Possible_ruc_branches,RUC_v_coord,f_rad);
else
    Possible_branch_vol = 0;
end

Max_vol_frac = (Sheetlet_myocyte_vol+Possible_branch_vol)/Total_volume;


%% Add branches to the RUC untill the required volume fraction is achieved

if Mesoscale == 1
    random_branch_id = randperm(length(Possible_ruc_branches));
    %to pick branches shuffled from all possible branches list
    
    %Adjusted to make use half of total branches temporaryly
    Required_V_frac = 1.36;
    if Max_vol_frac >= Required_V_frac
        target_vol_frac = Max_vol_frac;
    else
        target_vol_frac = Required_V_frac;
    end
    
    
    target_br_vol_frac = target_vol_frac - Vf_without_branching;
    target_br_vol = target_br_vol_frac*Total_volume;
    
    Final_branches = {};
    for br_id = 1:length(random_branch_id)
        Final_branches{end+1} = Possible_ruc_branches{random_branch_id(br_id)};
        cur_br_vol = Myo_vol_calc(Final_branches,RUC_v_coord,f_rad);
        if cur_br_vol>=target_br_vol
            break
        else
            continue
        end
    end
end


    

%% connectivity graph including inter sheetlet branching
 RUC_con_v1_branched = RUC_con_v1;
 RUC_con_v2_branched = RUC_con_v2;
 RUC_con_table_branched = RUC_con_table;  %connectivity table (type-cell)
 
 if Mesoscale ==1
     for branch_no1 = 1:length(Final_branches)
         RUC_con_v1_branched(end+1) = Final_branches{branch_no1}(1);
         RUC_con_v2_branched(end+1) = Final_branches{branch_no1}(2);
     end
     for branch_no2 = 1:length(Final_branches)
         RUC_con_table_branched{end+1} = Final_branches{branch_no2};
     end
 end

 %%
figure(2)

hold on;
XX_ruc=[];
YY_ruc=[];
ZZ_ruc=[];

for pt_g_ruc=1:length(RUC_v_coord)
    XX_ruc(end+1)=RUC_v_coord{pt_g_ruc}(1);
    YY_ruc(end+1)=RUC_v_coord{pt_g_ruc}(2);
    ZZ_ruc(end+1)=RUC_v_coord{pt_g_ruc}(3);
end

con_graph_ruc_branched=graph(RUC_con_v1_branched,RUC_con_v2_branched);
%con_graph_ruc_branched=graph(RUC_con_v1,RUC_con_v2);
pp_ruc_branched=plot(con_graph_ruc_branched);
pp_ruc_branched.XData=XX_ruc;
pp_ruc_branched.YData=YY_ruc;
pp_ruc_branched.ZData=ZZ_ruc; 

%%
save(append(ModelName,'_conData'));

%% functions used

function cell_vol = Myo_vol_calc(cell_array,cell_coord,cell_radius)
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

function Poly_area = Poly_area_calc(poly_connectivity,coordinates)
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

function d = min_dis_point_to_line(pt, v1, v2)
      a = v1 - v2;
      b = pt - v2;
      
      p = v2+a.*(dot(a,b)/(norm(a)^2));
      c1 = v2-p;
      c2 = v1-p;
      
      
      %Return correct distance if the point is located in the range of line
      %segment and return higher distance(5000) if the points are outside the
      %range of line segment
      if dot(c1,c2)>0
          d = 5000;
      else
        d = norm(cross(a,b)) / norm(a);
      end
end

function dis_pt = min_dis_line_point(v1_line,v2_line,Pt)

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

function dis_min = min_dis_2lines(l1p1,l1p2,l2p1,l2p2)

%To find the minimum distance between two lines connecting 4
%points

%Variable name pattern: l1p2 = line 1 point 2

%Consider line 1 and line 2 in in parametric form l1p1+s1v1 and l2p1+s2v2
%(v1 and v2 are vectors)

%v= (v1 crossproducr v2) is dirrection vector for line M which is
%orthogonal to both line1 and line2

%get projections of points l1p1(m1) and l2p1(m2) into M line and calculate the
%distance between those two points will give the minimum distance between
%two lines



v1=l1p2-l1p1;
v2=l2p2-l2p1;

mag_v1=norm(v1);
mag_v2=norm(v2);

U_v1 = v1/mag_v1;
U_v2 = v2/mag_v2;

vcross = cross(v1,v2);

%make closest point to the origin is first point of the line
if norm(l1p1)<norm(l1p2)
    v1_1=l1p1;
    v1_2=l1p2;
else
    v1_1=l1p2;
    v1_2=l1p1;
end

if norm(l2p1)<norm(l2p2)
    v2_1=l2p1;
    v2_2=l2p2;
else
    v2_1=l2p2;
    v2_2=l2p1;
end
%If two lines are parallel,check weather they overlap
v1 = v1_2-v1_1;
v2 = v2_2-v2_1;

denom=norm(vcross)^2;

if denom==0   %Two lines are parallel condition
    d1=dot(U_v1,v1_1-v2_1);
    d2=dot(U_v1,v1_1-v2_2);
    
    if d1<=0 && d2<=0  %line 2 before line 1 condition
        if abs(d1)<abs(d2)
            dis_min = norm(v1_1-v2_1);
        else
            dis_min = norm(v1_1-v2_2);
        end 
        
    elseif d1>=mag_v1 && d2>=mag_v2  %line2 after line1 condition
        if abs(d1)<abs(d2)
            dis_min = norm(v1_2-v2_1);
        else
            dis_min = norm(v1_2-v2_2);
        end
    else  %two line segments overlap condition
        v3=v1_1+d1*U_v1-v2_1;
        dis_min = norm(v3);
    end
    
else %Two lines are not parallel condition
    %find unit vectors to express any point in v1 in the form of v1_1+t1*Unit(V1)
    
    uu =  v2_1-v1_1;
    
    det1 = det([uu;U_v2;vcross]);
    det2 = det([uu;U_v1;vcross]);
    
    % Used reference - Sheetlet model/Resources
    v22 = dot(v2,v2);
    v11 = dot(v1,v1);
    v21 = dot(v2,v1);
    v21_1 = dot(uu,v1);
    v21_2 = dot(uu,v2);
    denom1 = v21*v21 - v22*v11;
    
    if isalmost(denom,0.0, 1e-4)
        t1=0;
        t2=(v11*t1 - v21_1)/v21;
    else
        t1 = (v21_2*v21-v22*v21_1)/denom1;
        t2 = (-v21_1*v21+v11*v21_2)/denom1;
    end
    
    t1 = max(min(t1,1),0);
    t2 = max(min(t2,1),0);
    
    t1=t1*mag_v1;
    t2=t2*mag_v2;
        


    
    %find closest points on two lines m1,m2
    %t1 = dot(cross(U_v2,vcross),(-v1_1+v2_1))/dot(vcross,vcross); %References in model folder
    %t1 = det1/denom
    m1 = v1_1+t1*U_v1;
    
    %t2 = dot(cross(U_v1,vcross),(-v1_1+v2_1))/dot(vcross,vcross);  %References in model folder
    %t2 = mag_v2*det2/denom
    m2 = v2_1+t2*U_v2;
    
    %Check whether m1 and m2 on line segments
    if dot((v1_1-m1),(v1_2-m1))<0 && dot((v2_1-m2),(v2_2-m2))<0
        %both points on line segments condition
        dis_min = norm(m1-m2);
    elseif dot((v1_1-m1),(v1_2-m1))>0 && dot((v2_1-m2),(v2_2-m2))<0
        %Only m2 on linesegment condition
        if dot((v1_1-m2),U_v1)<dot((v1_2-m2),U_v1)
            dis_min = norm(v1_1-m2);
        else
            dis_min = norm(v1_2-m2);
        end
    elseif dot((v1_1-m1),(v1_2-m1))<0 && dot((v2_1-m2),(v2_2-m2))>0
        %Only m1 on line segement condition
        if dot((v2_1-m1),U_v2)<dot((v2_2-m1),U_v2)
            dis_min = norm(v2_1-m1);
        else
            dis_min = norm(v2_2-m1);
        end
    else
        %Both points not on line segment condition
        if t1<0
            m1=v1_1;
        elseif t1>mag_v1
            m1=v1_2;
        end
        if t2<0
            m2=v2_1;
        elseif t2>mag_v2
            m2=v2_2;
        end
        dis_min = norm(m1-m2);
    end
end
end


function test = isalmost(a,b,tol)
%
% usage: test = isalmost(a,b,tol)
%
% tests if matrix a is approximately equal to b within a specified 
% tolerance interval (b-tol <= a <= b+tol)
%
% note:  if b is given as a scalar, all values in a are compared against
%        the scalar value b
%
% calls: none
%
% inputs:
%
% a(nr,nc) = matrix of data values to test
% b(nr,nc) = matrix of data values for comparison (or a single scalar value)
%      tol = tolerance used in computation
%
% outputs:
%
% test(nr,nc) = matrix of test results:
%
%        test(i,j) = 0 -> a(i,j) is not equal to b(i,j) (or is NaN)
%        test(i,j) = 1 -> a(i,j) is approximately equal to b(i,j)
%
%   author : James Crawford 
%   created 01/08/2007
%
%   history: v1.0 (01/08/2007)
%
% get length of input matrix a
[nr,nc] = size(a);
% check input for consistency
if ~all(size(a) == size(b))
   if all(size(b) == [1 1])
      % convert scalar value b to a matrix of size(a)
      b = b*ones(size(a));
   else
      disp('error: input arguments are inconsistent (isalmost.m)')
      disp('(b) must be a matrix of same size as (a) or a single value')
   end
end
one = ones(size(b));
% perform test
test = (a <= b+tol*one)&(a >= b-tol*one);
%
end


function coordinate_id = find_coord_id(coordinate_array,p_cord,sheetlet_id)

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

function sorted_point_list = sorted_pts_from_common_face(coord_list,line_vertex1,line_vertex2)

    point_dis_from_common_edge = [];
    
    % finding the seed points from sheetlet 1, within 15 micrometers from
    % the common face of two sheetlets
    
    for pts = 1:length(coord_list)
        point_f1_coord = coord_list{pts};
        point_f1_coord(1)=0;
        point_dis_from_common_edge(end+1) = min_dis_point_to_line(point_f1_coord,line_vertex1, line_vertex2);
    end
    
    [sorted_dis,sorted_id] = sort(point_dis_from_common_edge);
    
    sorted_point_list = {};
    for order = 1:length(sorted_id)
        sorted_point_list{1,order} = coord_list{sorted_id(order)};
        sorted_point_list{2,order} = sorted_id(order);
        sorted_point_list{3,order} = sorted_dis(order);
        
    end
    
end





