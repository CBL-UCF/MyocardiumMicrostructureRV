global RUC_w RUC_t RUC_l


%Thickness of considered region (> Sheeetlet thickness)
N = 30;        %number of sheetlets in YZ plane

%Creating mid points of the sheetlets, maintaining characteristic sizes in
%sheetlet dirrections and normal dirrections

if Mesoscale ==1
    slt_w_vor = 200;
    slt_t_vor = 45;
else
    slt_w_vor = RUC_w;
    slt_t_vor = RUC_t;
end

% Creating seeds for sheetlets. Here, Seeds were located in a way that create 
% voronoi cells with lesset thickness than width it makes cleavage planes to
% be approximately continuous representing the real structure

%seed_points = slt_mid_points(RUC_w,RUC_t,slt_w_vor,slt_t_vor);
seed_points = slt_mid_points_square(RUC_w,RUC_t,slt_w_vor,slt_t_vor);

seed_array = [];

for i=1:length(seed_points)
    seed_array(i,1)=seed_points{i}(1);
    seed_array(i,2)=seed_points{i}(2);
end
[Vx,Vy] = voronoi(seed_array(:,1),seed_array(:,2));
%tess = voronoi(seed_array(:,1),seed_array(:,2));
[vor_vertex,vor_cells] = voronoin(seed_array);




%Shrink voronoi cells to make space for cleavage planes
% shrink_ratio_boundary_y = [9.92,0.08]; %for sheetlet boundaries [9,1];%aligned latest [9.9,0.1]
% shrink_ratio_boundary_z = [8 12];%[13.8,6.2];
% shrink_ratio_basepoints_y = [17 3]; %to innitiate basepoints [18 2];%
% shrink_ratio_basepoints_z = [7 13];
% shrink_ratio_mid_y = [15 5]; %to use in generating mask 
% shrink_ratio_mid_z = [5 15];

% values used in 'Meso_sclhor_lclver_d10
shrink_ratio_boundary_y = [9.55,0.45]; 
shrink_ratio_boundary_z = [16.5 3.5];
shrink_ratio_basepoints_y = [19.4 0.6]; 
shrink_ratio_basepoints_z = [15 5];
shrink_ratio_mid_y = [14 6]; %to use in generating mask 
shrink_ratio_mid_z = [11 9];

Shrinked_tess = shrink_vor_cells(seed_array,vor_vertex,vor_cells,shrink_ratio_boundary_y,shrink_ratio_boundary_z);

%Shrink further the same voronoi cells to define area for myocyte vertexes
%in first face (Othervise when those are extruded in x dirrections,
%vertexes can move out of sheetlets)

Shrinked_tess_bp = shrink_vor_cells(seed_array,vor_vertex,vor_cells,shrink_ratio_basepoints_y,shrink_ratio_basepoints_z);
Shrinked_tess_mid = shrink_vor_cells(seed_array,vor_vertex,vor_cells,shrink_ratio_mid_y,shrink_ratio_mid_z);

%% select voronoi cells for the RUC 

%  Remove outer layer of voronoi cells (Shrinked)
%  Create an array of the data related to voronoi tesselations inside the
%  RUC domain

Shrinked_in_ruc = {}; % Include vertex coordinates
Shrinked_in_ruc_bp = {};
Shrinked_in_ruc_mid = {};

Nn = round(RUC_t/(slt_t_vor*1.1))+2;
Ns = round(RUC_w/(slt_w_vor*1.1))+2;

for p_no = 1: length(Shrinked_tess)
    if p_no < Nn
        continue
    elseif p_no > length(Shrinked_tess)-Nn
        continue
    elseif rem(p_no,Nn) == 0
        continue
    elseif rem(p_no,Nn) == 1
        continue
    else 
       Shrinked_in_ruc{end+1} =  Shrinked_tess{p_no};
    end
end

%  Remove outer layer of voronoi cells (non-shinked) (for use to identify
%  neighbouring sheetlets and shared boundarie)
%  Create an array of the data related to voronoi tesselations inside the
%  RUC domain

non_shrink_vor_ruc = {};
for p_no = 1: length(vor_cells)
    if p_no < Nn
        continue
    elseif p_no > length(vor_cells)-Nn
        continue
    elseif rem(p_no,Nn) == 0
        continue
    elseif rem(p_no,Nn) == 1
        continue
    else 
       non_shrink_vor_ruc{end+1} =  vor_cells{p_no};
    end
end

% %Get the corresponding vertex coordinates of non-shrinked vor cells in RUC
% non_shrinked_in_ruc = {};
% 
% for sss = 1: length(non_shrink_vor_ruc)
%     cell3 = non_shrink_vor_ruc{sss};
%     cell_v_coord = {};
%     for v5 = 1:length(cell3)
%        cell_v_coord{end+1} = [vor_vertex(v5,1),vor_vertex(v5,2)];
%     end
%     non_shrinked_in_ruc{end+1} = cell_v_coord;
% end

% repeat the same procedure for areas of seeding myocyte vertexes for the
% first face of RUC
for p_no = 1: length(Shrinked_tess_bp)
    if p_no < Nn
        continue
    elseif p_no > length(Shrinked_tess_bp)-Nn
        continue
    elseif rem(p_no,Nn) == 0
        continue
    elseif rem(p_no,Nn) == 1
        continue
    else 
       Shrinked_in_ruc_bp{end+1} =  Shrinked_tess_bp{p_no};
    end
end

for p_no = 1: length(Shrinked_tess_mid)
    if p_no < Nn
        continue
    elseif p_no > length(Shrinked_tess_mid)-Nn
        continue
    elseif rem(p_no,Nn) == 0
        continue
    elseif rem(p_no,Nn) == 1
        continue
    else 
       Shrinked_in_ruc_mid{end+1} =  Shrinked_tess_mid{p_no};
    end
end
%   Remove Shrinked polygons contain inf vertex coordinates (Since we cannot give
%   coordinates for the voronoi plots with inf as a coordinate value)

Shrinked_after_inf = {};
Shrinked_after_inf_bp = {};
Shrinked_after_inf_mid = {};

for poly = 1:length(Shrinked_in_ruc)
    coords_poly = Shrinked_in_ruc{poly};
    inf_indexes = isinf(coords_poly); %return 'true' for corresponding matrix index if the element is inf'
    N_inf = nonzeros(inf_indexes);
    
    %considered double have infinity coordinates if N_inf has any entities
    
    if ~isempty(N_inf)
        continue
    else
        Shrinked_after_inf{end+1} = Shrinked_in_ruc{poly};
    end        
end

for poly = 1:length(Shrinked_in_ruc_bp)
    coords_poly = Shrinked_in_ruc_bp{poly};
    inf_indexes = isinf(coords_poly); %return 'true' for corresponding matrix index if the element is inf'
    N_inf = nonzeros(inf_indexes);
    
    %considered double have infinity coordinates if N_inf has any entities
    
    if ~isempty(N_inf)
        continue
    else
        Shrinked_after_inf_bp{end+1} = Shrinked_in_ruc_bp{poly};
    end        
end

for poly = 1:length(Shrinked_in_ruc_mid)
    coords_poly = Shrinked_in_ruc_mid{poly};
    inf_indexes = isinf(coords_poly); %return 'true' for corresponding matrix index if the element is inf'
    N_inf = nonzeros(inf_indexes);
    
    %considered double have infinity coordinates if N_inf has any entities
    
    if ~isempty(N_inf)
        continue
    else
        Shrinked_after_inf_mid{end+1} = Shrinked_in_ruc_mid{poly};
    end        
end

%%  Creating 3D plot
figure(1)
title('Cardiac Meso-structure')

p_no = 1;

for vor_cell = 1:length(Shrinked_after_inf)
    v1_con =[];
    v2_con =[];
    
    vertex_X_3d = [];
    vertex_Y_3d = [];
    vertex_Z_3d = [];
    
    
    pt_corner = 1;
    for plane =1:2
        if plane == 1
            x_coord_vor = 0;
        else
            x_coord_vor = RUC_l;
        end
        
        no_of_v = length(Shrinked_after_inf{vor_cell});
        
        for pt = 1: no_of_v
            vertex_X_3d(end+1) = x_coord_vor;
            vertex_Y_3d(end+1) = Shrinked_after_inf{vor_cell}(pt,1);
            vertex_Z_3d(end+1) = Shrinked_after_inf{vor_cell}(pt,2);
            if pt_corner == no_of_v || pt_corner == 2*no_of_v
                v1_con(end+1) = pt_corner;
                if plane == 1
                    v2_con(end+1) = 1;
                else
                    v2_con(end+1) = no_of_v+1;
                end
            else
                v1_con(end+1) = pt_corner;
                v2_con(end+1) = pt_corner+1;
            end
            pt_corner = pt_corner+1;
        end
    end
    
    for edge_v = 1:no_of_v
        v1_con(end+1) = edge_v;
        v2_con(end+1) = edge_v + no_of_v;
    end
        
    
    
    XX=vertex_X_3d;
    YY=vertex_Y_3d;
    ZZ=vertex_Z_3d;



    g = graph(v1_con,v2_con);
    pp=plot(g);
    pp.XData=XX;
    pp.YData=YY;
    pp.ZData=ZZ;
    hold on;
    
end
axis equal

%% Functions

function mid_point_list = slt_mid_points(RUC_w,RUC_t,Slt_w, Slt_t)

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

function shrinked_tess_coord = shrink_vor_cells(seed_array,v,vor_cells,shrink_ratio_y,shrink_ratio_z)

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

function mid_point_list = slt_mid_points_square(RUC_w,RUC_t,Slt_w, Slt_t)

mid_point_list = {};

Ns = round(RUC_w/(Slt_w*1.1))+2;
Nn = round(RUC_t/(Slt_t*1.1))+2;



for pts = 1:Ns
    for ptn = 1:Nn
        if rem(pts,2)==1
            z=0;
        else
            z = 0; %0.5*Slt_t;
        end
        if rem(ptn,2)==1
            y=(pts+0.5)*Slt_w;
            z=z+(ptn+0.5)*Slt_t;
        else
            y=(0+pts+0.5)*Slt_w; %add 0.5 of sheetlet width to shift the layer
            z=z+(ptn+0.5)*Slt_t;
        end
        point = [y,z];
        mid_point_list{end+1} = point;
    end
end

% for pts = 1:Ns
%     for ptn = 1:Nn
%         if rem(pts,2)==1
%             z=0;
%         else
%             z = 0*Slt_t;
%         end
%         if rem(ptn,2)==1
%             y=(pts+0.5)*Slt_w;
%             z=z+(ptn+0.5)*Slt_t;
%         else
%             y=(0+pts+0.5)*Slt_w; %add 0.5 of sheetlet width to shift the layer
%             z=z+(ptn+0.5)*Slt_t;
%         end
%         point = [y,z];
%         mid_point_list{end+1} = point;
%     end
% end
        
end