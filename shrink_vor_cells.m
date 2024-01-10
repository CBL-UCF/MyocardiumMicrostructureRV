function shrinked_tess_coord = shrink_vor_cells(seed_array,v,vor_cells,shrink_ratio)

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
    aa = shrink_ratio(1);
    bb = shrink_ratio(2);
    nw_poly_coord =[];
    for point_v = 1:length(ex_poly_coord)
        nw_poly_coord(point_v,1) = (bb*seed_array(vor,1)+aa*ex_poly_coord(point_v,1))/(aa+bb);
        nw_poly_coord(point_v,2) = (bb*seed_array(vor,2)+aa*ex_poly_coord(point_v,2))/(aa+bb);
    end
    
    shrinked_tess_coord{end+1} = nw_poly_coord;
end
    