clc;clear;


global RUC_w RUC_t RUC_l


%Thickness of considered region (> Sheeetlet thickness)
N = 30;        %number of sheetlets in YZ plane

%Creating mid points of the sheetlets, maintaining characteristic sizes in
%sheetlet dirrections and normal dirrections

slt_w_vor = 250;
slt_t_vor = 100;

% Creating seeds for sheetlets. Here, Seeds were located in a way that create 
% voronoi cells with lesset thickness than width it makes cleavage planes to
% be approximately continuous representing the real structure

seed_points = slt_mid_points(RUC_w,RUC_t,slt_w_vor,slt_t_vor);
seed_array = [];

for i=1:length(seed_points)
    seed_array(i,1)=seed_points{i}(1);
    seed_array(i,2)=seed_points{i}(2);
end
[Vx,Vy] = voronoi(seed_array(:,1),seed_array(:,2));
%tess = voronoi(seed_array(:,1),seed_array(:,2));
[v,vor_cells] = voronoin(seed_array);



%Shrink voronoi cells to make space for cleavage planes
shrink_ratio = [9.5,0.5];
Shrinked_tess = shrink_vor_cells(seed_array,v,vor_cells,shrink_ratio);


%  Remove outer layer of voronoi cells

Shrinked_in_ruc = {};
Nn = round(RUC_t/(slt_t_vor*1.1))+2;


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


%   Remove Shrinked polygons contain inf vertex coordinates (Since we cannot give
%   coordinates for the voronoi plots with inf as a coordinate value)

Shrinked_after_inf = {};

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





%%  Creating plot
figure(3)
title('Cardiac Meso -  RUC')

p_no = 1;

for vor_cell = 1:length(Shrinked_after_inf)
    v1_con =[];
    v2_con =[];
    
    vertex_Y = [];
    vertex_Z = [];
    
    pt_corner = 1;
    
    for pt = 1: length(Shrinked_after_inf{vor_cell})
        vertex_Y(end+1) = Shrinked_after_inf{vor_cell}(pt,1);
        vertex_Z(end+1) = Shrinked_after_inf{vor_cell}(pt,2);
        if pt_corner == length(Shrinked_after_inf{vor_cell})
            v1_con(end+1) = pt_corner;
            v2_con(end+1) = 1;
        else
            v1_con(end+1) = pt_corner;
            v2_con(end+1) = pt_corner+1;
        end
        pt_corner = pt_corner+1;
    end
    
    
    XX=zeros(1,length(vertex_Y));
    YY=vertex_Y;
    ZZ=vertex_Z;



    g = graph(v1_con,v2_con);
    pp=plot(g);
    pp.XData=XX;
    pp.YData=YY;
    pp.ZData=ZZ;
    hold on;
    
end
% axis equal


%%  Creating 3D plot
figure(4)
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