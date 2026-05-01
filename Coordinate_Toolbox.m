classdef Coordinate_Toolbox

    % Coordinate_Toolbox
    % Utility class providing geometric operations for coordinate handling,
    % transformations, distance calculations, and connectivity processing
    % for sheetlet-based microstructural models.
    
    methods(Static)
        
        %% Function to get X,Y,Z coordinate list from a cell
        
        function [coordinates] = cell_to_double_coordinates(cell)

            % cell_to_double_coordinates
            % Converts a cell array of 3D coordinate vectors into a numeric
            % Nx3 matrix for easier processing.
            %
            % Input:
            %   cell        - Cell array where each entry is [x,y,z]
            %
            % Output:
            %   coordinates - Nx3 numeric matrix of coordinates

            XX_coord=[];
            YY_coord=[];
            ZZ_coord=[];

            for pt=1:length(cell)
                XX_coord(end+1)=cell{pt}(1);
                YY_coord(end+1)=cell{pt}(2);
                ZZ_coord(end+1)=cell{pt}(3);
            end
            coordinates(:,1) = XX_coord;
            coordinates(:,2) = YY_coord;
            coordinates(:,3) = ZZ_coord;
        end
        
        %% function to find coordinate modification amounts to translate

        %%coordinates to 'min point' befor develop mask
        
        function [coord_details] = Modification_amount_and_dimension(polygons,point_list,min_point)

            % Modification_amount_and_dimension
            % Computes translation offsets required to shift all coordinates
            % so that the minimum point aligns with a given reference point.
            % Also computes domain dimensions.
            %
            % Inputs:
            %   polygons   - Voronoi polygon list (cell array of vertices)
            %   point_list - List of 3D coordinates
            %   min_point  - Desired minimum coordinate [x,y,z]
            %
            % Output:
            %   coord_details:
            %       row 1 -> [dx, dy, dz] translation amounts
            %       row 2 -> [Lx, Ly, Lz] domain dimensions

            voronoi_vertex_y = [];
            voronoi_vertex_z = [];
            for i =1:length(polygons)
                list = polygons{i};
                for j = 1:size(list,1)
                    voronoi_vertex_y(end+1) = list(j,1);
                    voronoi_vertex_z(end+1) = list(j,2);
                end
            end

            XX_ruc_mask=[];
%             YY_ruc_mask=[];
%             ZZ_ruc_mask=[];

            for pt_g_ruc=1:length(point_list)
                XX_ruc_mask(end+1)=point_list{pt_g_ruc}(1);
%                 YY_ruc_mask(end+1)=point_list1{pt_g_ruc}(2);
%                 ZZ_ruc_mask(end+1)=point_list1{pt_g_ruc}(3);
            end

            modify_x = min_point(1)-min(XX_ruc_mask);
            modify_y = min_point(2)-min(voronoi_vertex_y);
            modify_z = min_point(3)-min(voronoi_vertex_z);
            
            domain_x= max(XX_ruc_mask)-min(XX_ruc_mask);
            domain_y =max(voronoi_vertex_y)-min(voronoi_vertex_y);
            domain_z =max(voronoi_vertex_z)-min(voronoi_vertex_z);
            
                       
            coord_details(1,:) = [modify_x,modify_y,modify_z];
            coord_details(2,:) = [domain_x,domain_y,domain_z];
        end
        
        %% Modify all coordinates in a list to a minimum point
        function [Modified_RUC_coord] =  Modify_RUC_coord_mask(point_list1, modify_amount)

            % Modify_RUC_coord_mask
            % Applies a translation to all coordinates in a list so that the
            % domain is shifted to a desired reference position.
            %
            % Inputs:
            %   point_list1   - Cell array of 3D coordinates
            %   modify_amount - Translation vector [dx, dy, dz]
            %
            % Output:
            %   Modified_RUC_coord - Translated coordinate list

            XX_ruc_mask=[];
            YY_ruc_mask=[];
            ZZ_ruc_mask=[];

            for pt_g_ruc=1:length(point_list1)
                XX_ruc_mask(end+1)=point_list1{pt_g_ruc}(1);
                YY_ruc_mask(end+1)=point_list1{pt_g_ruc}(2);
                ZZ_ruc_mask(end+1)=point_list1{pt_g_ruc}(3);
            end

            X_modify = modify_amount(1);
            Y_modify = modify_amount(2);
            Z_modify = modify_amount(3);

            XX_ruc_modified = XX_ruc_mask+X_modify;
            YY_ruc_modified = YY_ruc_mask+Y_modify;
            ZZ_ruc_modified = ZZ_ruc_mask+Z_modify;


            Modified_RUC_coord ={};
            for i = 1: length(XX_ruc_modified)
                Modified_RUC_coord{end+1} = [XX_ruc_modified(i) YY_ruc_modified(i) ZZ_ruc_modified(i)];
            end
        end
        
        
        %% function to modify sheetlet coordinates to compatible with mask 
        function [Modified_sheetlets] = Modify_shtlt_coord(shtlt_list,modify_amount)

            % Modify_shtlt_coord
            % Applies translation to 2D sheetlet polygon coordinates to align
            % them with the mask coordinate system.
            %
            % Inputs:
            %   shtlt_list    - Cell array of polygon coordinates (y,z)
            %   modify_amount - Translation vector [dx, dy, dz]
            %
            % Output:
            %   Modified_sheetlets - Translated polygons

            Modified_sheetlets = cell(1,length(shtlt_list));
            
            for i = 1:length(shtlt_list)
                shtlt1 = shtlt_list{i};
                Modified_sheetlets{i} = [shtlt1(:,1)+modify_amount(2),shtlt1(:,2)+modify_amount(3)];
            end
        end
        
        %% function to calculate minimum distance between line segment and a point in 3D space
        function dis_pt = min_dis_line_point(v1_line,v2_line,Pt)

            % min_dis_line_point
            % Computes minimum distance between a 3D point and a line segment.
            % Returns distance only if projection lies within segment bounds,
            % otherwise returns a large value.
            %
            % Inputs:
            %   v1_line, v2_line - Endpoints of the segment
            %   Pt               - Query point
            %
            % Output:
            %   dis_pt           - Minimum distance

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
        
        
        %% function to calculate minimum distance between two line segments in 3D space
        function dis_min = min_dis_2lines(l1p1,l1p2,l2p1,l2p2)

            % min_dis_2lines
            % Computes the minimum distance between two line segments in 3D,
            % accounting for parallel and skew configurations.
            %
            % Inputs:
            %   l1p1, l1p2 - Line 1 endpoints
            %   l2p1, l2p2 - Line 2 endpoints
            %
            % Output:
            %   dis_min    - Minimum distance between segments



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

                if Matrix_Toolbox.isalmost(denom,0.0, 1e-4)
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
                        %return the distance between closes points
                    end
                end
        end

        %% function to get the arrays of sheetletwise interceonnected cells
        
        function [shtlt_wise_connection_list] = find_connected_cells(RUC_con_table_sheetlet_wise,Final_branches)
            
            % find_connected_cells
            % Identifies neighboring/connected cells within each sheetlet based on
            % shared vertices, and also includes connections to branch elements.
            %
            % Inputs:
            %   RUC_con_table_sheetlet_wise - Cell array containing sheetlet-wise
            %                                  connectivity tables.
            %   Final_branches              - Cell array containing final branch
            %                                  connectivity pairs.
            %
            % Output:
            %   shtlt_wise_connection_list  - Sheetlet-wise list of connected cell
            %                                  and branch indices.

            
            glob_cell_add =0;
            n_shtlt = length(RUC_con_table_sheetlet_wise);
            shtlt_wise_connection_list = cell(n_shtlt,1);

            n_branch = length(Final_branches);
            branch_array = zeros(n_branch,2);
            for b=1:n_branch
                branch_array(b,1) = Final_branches{b}(1);
                branch_array(b,2) = Final_branches{b}(2);
            end

            sheetlet_all_cells = 0;
            for stlt = 1:n_shtlt
                sheetlet_all_cells=sheetlet_all_cells+length(RUC_con_table_sheetlet_wise{stlt});
            end


            for shtlt = 1:length(RUC_con_table_sheetlet_wise)
                conn_list = RUC_con_table_sheetlet_wise{shtlt};
                cell_c = length(conn_list);
                shtlt_con_array = zeros(cell_c,2);
                for i = 1:cell_c
                    shtlt_con_array(i,1) = conn_list{i}(1);
                    shtlt_con_array(i,2) = conn_list{i}(2);
                end

                connected_cells = cell(cell_c,1);
                for j = 1:cell_c
                    vertexes = shtlt_con_array(j,:);
                    v1_con_ids1 = find(shtlt_con_array(:,1) == vertexes(1));
                    v1_con_ids2 = find(shtlt_con_array(:,2) == vertexes(1));
                    v1_con_cells = union(v1_con_ids1,v1_con_ids2);
                    v1_con_cells_glob = v1_con_cells +glob_cell_add;

                    v2_con_ids1 = find(shtlt_con_array(:,1) == vertexes(2));
                    v2_con_ids2 = find(shtlt_con_array(:,2) == vertexes(2));
                    v2_con_cells = union(v2_con_ids1,v2_con_ids2);
                    v2_con_cells_glob = v2_con_cells +glob_cell_add;

                    v1_br_id1 = find(branch_array(:,1) == vertexes(1));
                    v1_br_id2 = find(branch_array(:,2) == vertexes(1));
                    v1_con_brs = union(v1_br_id1,v1_br_id2);
                    v1_con_brs_glob = v1_con_brs+sheetlet_all_cells;


                    v2_br_id1 = find(branch_array(:,1) == vertexes(2));
                    v2_br_id2 = find(branch_array(:,2) == vertexes(2));
                    v2_con_brs = union(v2_br_id1,v2_br_id2);
                    v2_con_brs_glob = v2_con_brs+sheetlet_all_cells;


                    all_con_cells = union(v1_con_cells_glob,v2_con_cells_glob);
                    all_con_branch = union(v1_con_brs_glob,v2_con_brs_glob);
                    all_con = union(all_con_cells,all_con_branch);

                    connected_cells{j,1} = all_con;
                end
                shtlt_wise_connection_list{shtlt,1} = connected_cells;
                glob_cell_add = glob_cell_add + cell_c;
            end            

        
        end
        
        %% function to generate global connectivity list
        
        function [Global_connection_list, N_cells_before_branch] = glob_con_list(shtlt_wise_connection_list)

            % glob_con_list
            % Combines sheetlet-wise connection lists into one global connection
            % list and counts the number of cells before adding branch elements.
            %
            % Input:
            %   shtlt_wise_connection_list - Sheetlet-wise cell connectivity list.
            %
            % Outputs:
            %   Global_connection_list     - Global list of connected cell indices.
            %   N_cells_before_branch      - Total number of non-branch cells.

            N_shtlt = length(shtlt_wise_connection_list);
            N_cells_before_branch = 0;
            Global_connection_list = {};
            for i =1:N_shtlt
                shtlt_cells = shtlt_wise_connection_list{i,1};
                cell_count = length(shtlt_cells);
                N_cells_before_branch = N_cells_before_branch+cell_count;
                for j = 1:cell_count
                    Global_connection_list{end+1} = shtlt_cells{j};
                end
            end
        end
        
        %% Function to rotate voronoi polygons by an angle around x axis

        function [Rotated_poly_list] = Rotate_2D_poly(poly_list,angle)

            % Rotate_2D_poly
            % Rotates a list of 2D polygon coordinates by a specified angle using
            % the standard 2D rotation matrix.
            %
            % Inputs:
            %   poly_list         - Cell array of polygon vertex coordinates.
            %   angle             - Rotation angle in radians.
            %
            % Output:
            %   Rotated_poly_list - Cell array of rotated polygon coordinates.

            R = [cos(angle), -sin(angle);sin(angle),cos(angle)];
            
            n_poly = length(poly_list);
            Rotated_poly_list = cell(1,n_poly);
            
            for i = 1:n_poly
                poly  = poly_list{i};
                nvec = size(poly,1);
                rotated_poly = [];
                for j=1:nvec
                    vec = poly(j,:)';
                    rot_vec = R*vec;
                    rotated_poly(j,:) = rot_vec;
                end
                Rotated_poly_list{1,i} = rotated_poly;
            end

         
        end
        
        %% Function to create stripe regions (vertical) from a square
        function [Regions] = CreateDevidedRegions(Domain, n_stripes)

                % CreateDevidedRegions
                % Divides a rectangular domain into equal vertical stripe regions.
                %
                % Inputs:
                %   Domain    - Domain dimensions [width, height].
                %   n_stripes - Number of vertical regions to create.
                %
                % Output:
                %   Regions   - Cell array of polygon coordinates for each stripe.

            Regions =cell(1,n_stripes);

            for i = 1:n_stripes
                w_rgn = Domain(1)/n_stripes;
                t = Domain(2);
                r_temp = [w_rgn*(i-1),0;
                    w_rgn*i,0;
                    w_rgn*i,t;
                    w_rgn*(i-1),t];
                Regions{1,i} = r_temp;
            end
        end


        %% Min distance between two segments

        function [distance, varargout] = DistBetween2Segment(p1, p2, p3, p4)

            % DistBetween2Segment
            % Computes minimum distance between two line segments and optionally
            % returns closest points and connecting vector.
            %
            % Inputs:
            %   p1, p2 - Segment 1 endpoints
            %   p3, p4 - Segment 2 endpoints
            %
            % Outputs:
            %   distance - Minimum distance
            %   varargout:
            %       (1) vector between closest points
            %       (2) closest point on segment 1
            %       (3) closest point on segment 2

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
        
        
        
        
        
        
    end
end