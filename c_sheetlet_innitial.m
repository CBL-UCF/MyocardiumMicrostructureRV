% reset(obj)
% release(obj)


%Sheetlet model
%% Define inputs

% Dimensions of innitial box to get myocyte network inside sheetlets

L_SLT = SLT_L; 
W_SLT = SLT_W;
T_SLT = SLT_T;



%Cardiomyocyte dimensions
%d_f = 15;
%d_f = 10 + (20-10)*betarnd(3,3,1,500);  %Creating beta randomly distributed
                                        %array for cardiomyocite diameters 
                                        %between 10-20  micrometers (3,3 
                                        %gives normal distributionbut not 
                                        %smooth)
f_rad = 10/2;

% figure (1)
% histogram(d_f,50)

l_f = 100 + (150-100)*betarnd(3,3,1,500); %Creating beta randomly distribut
                                          %ed array for cardiomyocite 
                                          %lengths between 100-150 
                                          %micrometers (3,3 gives normal 
                                          %distributionbut not smooth)

%Volume fractions 

Vf_myocite = 0.3;
Vf_ecm =0.7;


%Sheetlet Orientation

ha = 0; %Helical angle in radian
ta = 0; %Transverse angle in radian





%% Create points on face 1 (Without layers)

%Number of fibre along fibre dirrection and within one cross section
Nf_l = round(L_SLT/100);  %Considered 100 nanometer innitial distance between two myocites in fibre dirrection

vor_area = polyarea(Y_bp_v,Z_bp_v);
Nf_f =  round((0.95*(vor_area)/(pi*f_rad^2)),0);



%Generate points on surface 1
Point_list = {};
face1_pt_list ={};
Dis = @(a1,a2,b1,b2,c1,c2) sqrt(((a1-a2)^2+(b1-b2)^2+(c1-c2)^2));


p_face=0;
br_point=0;


rng(rng_count) % to innitialize the random number generating global stream
while p_face<Nf_f 

    %P_x = -15+30*rand(1,1); % keeping variation of -15 to 30 micrometers for each vertex plane along x dirrection
    P_x = 0;
    P_y = min_Y + rand(1,1)*W_SLT;
    P_z = min_Z + rand(1,1)*T_SLT;
    
    %make sure that the points are inside the voronoi cell
    
    inout_vor = inpolygon(P_y,P_z,Y_bp_v,Z_bp_v);
    if ~inout_vor
        continue
    end


    br_point=br_point+1; 

    D_list=[];
    for c = 1:length(Point_list)
        D = Dis(P_x,Point_list{c}(1),P_y,Point_list{c}(2),P_z,Point_list{c}(3));
        D_list(end+1) = D;

    end

    %Introducing a breakpoint if number of cardiomyocytes across the space
    %cannot reach 50% of the area

    if br_point==400000
        fprintf('Randomization iterations reached 100000 and generated points = %i / %i \n',length(Point_list),Nf_f)
        break
    end

    %Make sure cardiomycytes in the region do not overlap each other  
    if min(D_list)<2*f_rad
        continue    
    else 
        Point_list{end+1} = [P_x,P_y,P_z];
        face1_pt_list{end+1} = [P_x,P_y,P_z];
        p_face=p_face+1;
    end
end
Nf_f = length(face1_pt_list);
area_frac = (Nf_f*pi*f_rad^2)/vor_area;

area_frac_list(end+1) = area_frac;
 
%% Extruding points

% Nf_l - Number of vertex planes(faces) along x dirrection

facevise_p_list={face1_pt_list};

% create empty cells in a cell array for each plane along x dirrection 
for cl=1:Nf_l
    empty_cell={};
    facevise_p_list{end+1}=empty_cell;
end

%Consider each point in the first face of the sheetlet and extrude those
%points Nf_l times, checking minimum distance of fibre diameter with any other point
%in the whole point list


max_var_yz = f_rad/3; %f_rad/3; %maximum variation in x and y coodinates of cell center line along a line = fibre diameter/3
% max_var_yz = f_rad*1.5;

rng_count = rng_count+5;
rng(rng_count)

for pt=1:length(face1_pt_list)
    Base_point = face1_pt_list{pt};
    
    ln=1;
    
    
    br_point=0;
    while ln<=Nf_l
        Px_new = Base_point(1)+ln*100-30+60*rand(1,1);
        Py_new = Base_point(2)-(max_var_yz/2)+max_var_yz*rand(1,1);     
        Pz_new = Base_point(3)-(max_var_yz/2)+max_var_yz*rand(1,1);
        
        D_list=[];
        for c = 1:length(Point_list)
            D = Dis(Px_new,Point_list{c}(1),Py_new,Point_list{c}(2),Pz_new,Point_list{c}(3));
            D_list(end+1) = D;
        end
        
        br_point = br_point+1;
        
        if br_point==400000
            fprintf('Randomization iterations reached 400000 and extruded points = %i / %i for %i base point \n',ln,Nf_l, pt)
            for it = 1:Nf_l-ln+1
                Point_list{end+1} = [NaN,NaN,NaN];
                facevise_p_list{ln+1}{end+1} = [NaN,NaN,NaN]; %Add points to the facevise point list
                p_face=p_face+1;
                ln=ln+1;
            end
            break
        end
        
       if min(D_list)<2*f_rad
            continue    
       else 
            Point_list{end+1} = [Px_new,Py_new,Pz_new];
            facevise_p_list{ln+1}{end+1} = [Px_new,Py_new,Pz_new]; %Add points to the facevise point list
            p_face=p_face+1;
            ln=ln+1;
       end
    end
end

%% construct the connectivity table

con_table = {};

linewise_table ={}; %Array of  arrays containing a array of nodes for every continuing line of fibres
%create empty cells in a cell array for each extruding basepoint line 
for c2=1:Nf_f
    empty_cell={};
    linewise_table{end+1}=empty_cell;
end

%innitiate variable to give index to each graph vertex
in_v = 0;

%innitiate connectivity of vertex lists of the graph
graph_v1=[];
graph_v2=[];

%innitiate a list of points that actually got a connection after all checks
connected_p_list = {};

break_extending_line = 0;
min_line_dis=[];
for bp=1:Nf_f    %id of the basepoint that is extending
    for v = 1:Nf_l   %id of the cell number in lingitudinal dirrection(x) coreespoinding correspoinding to basepoint
        v1=facevise_p_list{v}{bp};
        v2=facevise_p_list{v+1}{bp};
        
        %Checking minimum distance of the line connecting v1 and v2 with 
        %all other lines already cerated
        
        line_dis_list=[1000];
        
        for bp_c = 1:bp-1
            for ln_c = 1:length(linewise_table{bp_c})
                v1_c=linewise_table{bp_c}{ln_c}{1};
                v2_c=linewise_table{bp_c}{ln_c}{2};
                
                ddd = min_dis_2lines(v1,v2,v1_c,v2_c); %Calling function to calculate distance
                %ddd = DistBetween2Segment(v1,v2,v1_c,v2_c);
                line_dis_list(end+1)=ddd;
                break_extending_line = 0;
                if ddd < 2*f_rad
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

%% Creating branchs


con_table_branched = con_table;

graph_v1_branched = graph_v1;
graph_v2_branched = graph_v2;

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
    
possible_branch_table = {};
branch_table_coord = {};

Vf_before_branch = fibreVolfrac(linewise_table,f_rad,L_SLT,vor_area);

 
branch_bp_count = round(1*length(face1_pt_list));  %Create branches in 

no_of_nodes = length(connected_p_list);
node_branch_count = zeros(1,no_of_nodes);

for Bp = 1:branch_bp_count
    Bp1_coord = face1_pt_list{1,Bp};
     
    
     
    dis_bp1_bp2 = [];
    for pid = 1:length(face1_pt_list)
        dis = norm([0,face1_pt_list{pid}(2),face1_pt_list{pid}(3)]-[0,Bp1_coord(2),Bp1_coord(3)]);
        if dis == 0
           dis_bp1_bp2(end+1) = 100000; % to keep the same base point at the end on the array
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
        BP2_coord = face1_pt_list{I}; %This is the I'th closest base point in sheetlet 2(face1) to the selected basepoint in sheetlet 1

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
            

            cb_v1_graph_id = find_coord_id_shtlt(connected_p_list,cb_v1);
            cb_v2_graph_id = find_coord_id_shtlt(connected_p_list,cb_v2);
            
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
                     
                     %D_cells = min_dis_2lines(cb_v1,cb_v2,connected_p_list{ch_v1_in},connected_p_list{ch_v2_in});
                     D_cells = DistBetween2Segment(cb_v1,cb_v2,connected_p_list{ch_v1_in},connected_p_list{ch_v2_in});
                     %if D_cells < 1e-10
                     if cb_v1_graph_id == ch_v1_in || cb_v2_graph_id == ch_v1_in || cb_v1_graph_id == ch_v2_in || cb_v2_graph_id == ch_v2_in
                         cell_spcaes(end+1) = 100000;
                     else
                         cell_spcaes(end+1) = D_cells;
                     end

                end
                % N_current_branches = length(possible_branch_table);
                % if N_current_branches > 0
                %     for cell2 = 1: length(possible_branch_table)
                %         ch_v1_in = possible_branch_table{cell2}(1);
                %         ch_v2_in = possible_branch_table{cell2}(2);
                %         %D_cells = min_dis_2lines(cb_v1,cb_v2,connected_p_list{ch_v1_in},connected_p_list{ch_v2_in});
                %         D_cells = DistBetween2Segment(cb_v1,cb_v2,connected_p_list{ch_v1_in},connected_p_list{ch_v2_in});
                %         %if D_cells < 1e-10
                %         if cb_v1_graph_id == ch_v1_in || cb_v2_graph_id == ch_v1_in || cb_v1_graph_id == ch_v2_in || cb_v2_graph_id == ch_v2_in
                %             cell_spcaes(end+1) = 100000;
                % 
                %         else
                %             cell_spcaes(end+1) = D_cells;
                %         end
                %     end
                % end
                if min(cell_spcaes)<2*f_rad
                     branches_overlap = 1;
                     continue;
                end
            
            end
            
            %cHECK WHETHER THE SAME NODE HAS MORE THAN TWO BRANCHES
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

 %% Calculate new volume fraction
 
 Max_branch_vol = 0;
 for lines = 1:length(branch_table_coord)
     Max_branch_vol = Max_branch_vol+pi*f_rad^2*norm(branch_table_coord{lines}{2}-branch_table_coord{lines}{1});
 end
 Total_volume = L_SLT*vor_area;
 Max_Vf_branched= 100*Max_branch_vol/Total_volume;
 
Max_Vf_branched = Vf_before_branch+Max_Vf_branched
 
%% Add branches to the RUC untill the required volume fraction is achieved

random_branch_id = randperm(length(possible_branch_table));
%to pick branches shuffled from all possible branches list

%Adjusted to make use half of total branches temporaryly
Required_V_frac = 0.97;
if Max_Vf_branched >= Required_V_frac
    target_vol_frac = Max_Vf_branched;
else
    target_vol_frac = Required_V_frac;
end


target_br_vol_frac = target_vol_frac - Vf_before_branch;
target_br_vol = target_br_vol_frac*Total_volume;

Final_branches = {};
for br_id = 1:length(random_branch_id)
    Final_branches{end+1} = possible_branch_table{random_branch_id(br_id)};
    cur_br_vol = Myo_vol_calc(Final_branches,connected_p_list,f_rad);
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


%% Graphical representation of connectivity network

figure(1)

hold on;
XX=[];
YY=[];
ZZ=[];

for pt_g=1:length(connected_p_list)
    XX(end+1)=connected_p_list{pt_g}(1);
    YY(end+1)=connected_p_list{pt_g}(2);
    ZZ(end+1)=connected_p_list{pt_g}(3);
end

con_graph=graph(graph_v1_branched,graph_v2_branched);
pp=plot(con_graph);
pp.XData=XX;
pp.YData=YY;
pp.ZData=ZZ;

%% Record coordinates of each vertexes in the RUC and connectivity  

%addidng coordinates of vertexes in each sheetlet as cells

 
coordinate_array{end+1} = connected_p_list;
connectivity_array{end+1} = con_table_branched_final;
connectivity_V1_array{end+1} = graph_v1_branched;
connectivity_V2_array{end+1} = graph_v2_branched;
face1_point_array{end+1} = face1_pt_list;
linewise_point_array{end+1} = linewise_point_list;
        
        
 
%% Functions

function Vf = fibreVolfrac(con_list,fibre_r,RUC_l,vor_area)  %Calculating total fibre volume fraction in the RUC acoording to connectivity table

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
    
    if isalmost(denom,0.0, 1e-6)
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

function coordinate_id = find_coord_id_shtlt(coordinate_array,p_cord)


for cd_no = 1:length(coordinate_array)
    if p_cord == coordinate_array{cd_no}
        coordinate_id = cd_no;
    else
        continue
    end
end

end

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

function [distance varargout] = DistBetween2Segment(p1, p2, p3, p4)
    u = p1 - p2;
    v = p3 - p4;
    w = p2 - p4;
    
    a = dot(u,u);
    b = dot(u,v);
    c = dot(v,v);
    d = dot(u,w);
    e = dot(v,w);
    D = a*c - b*b;
    sD = D;
    tD = D;
    
    SMALL_NUM = 0.00000001;
    
    % compute the line parameters of the two closest points
    if (D < SMALL_NUM)  % the lines are almost parallel
        sN = 0.0;       % force using point P0 on segment S1
        sD = 1.0;       % to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    else                % get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0)   % sc < 0 => the s=0 edge is visible       
            sN = 0.0;
            tN = e;
            tD = c;
        elseif (sN > sD)% sc > 1 => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        end
    end
    
    if (tN < 0.0)            % tc < 0 => the t=0 edge is visible
        tN = 0.0;
        % recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        elseif (-d > a)
            sN = sD;
        else
            sN = -d;
            sD = a;
        end
    elseif (tN > tD)       % tc > 1 => the t=1 edge is visible
        tN = tD;
        % recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        elseif ((-d + b) > a)
            sN = sD;
        else 
            sN = (-d + b);
            sD = a;
        end
    end
    
    % finally do the division to get sc and tc
    if(abs(sN) < SMALL_NUM)
        sc = 0.0;
    else
        sc = sN / sD;
    end
    
    if(abs(tN) < SMALL_NUM)
        tc = 0.0;
    else
        tc = tN / tD;
    end
    
    % get the difference of the two closest points
    dP = w + (sc * u) - (tc * v);  % = S1(sc) - S2(tc)
    distance = norm(dP);
    outV = dP;
    
    varargout(1) = {outV};      % vector connecting the closest points
    varargout(2) = {p2+sc*u};   % Closest point on object 1 
    varargout(3) = {p4+tc*v};   % Closest point on object 2
    
end


