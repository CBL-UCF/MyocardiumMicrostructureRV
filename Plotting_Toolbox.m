classdef Plotting_Toolbox

    methods(Static)

        function a = monitor_sheetlet_boundaries_in3D(fig_no, fig_title, polygon_list, RUC_l)
            a = true;

            figure(fig_no)
            title(fig_title)            
            
            p_no = 1;            
            for vor_cell = 1:length(polygon_list)
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
                    
                    no_of_v = length(polygon_list{vor_cell});
                    
                    for pt = 1: no_of_v
                        vertex_X_3d(end+1) = x_coord_vor;
                        vertex_Y_3d(end+1) = polygon_list{vor_cell}(pt,1);
                        vertex_Z_3d(end+1) = polygon_list{vor_cell}(pt,2);
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

        end

        % ----------------------------------------------------------------%

        function fig1 = monitor_insheetlet_microstructure(fig_no, connected_p_list, graph_v1_branched,graph_v2_branched)
            
            fig1 = 0;
            figure(fig_no)

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
        end

    end

end