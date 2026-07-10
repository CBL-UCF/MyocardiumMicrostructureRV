classdef ROI_cutout_toolbox

    methods(Static)

        function [newPts, newConn, parentElem, tRange] = clipLineElementsToRUC(points, conn, rucMin, rucDim, tol)
            % clipLineElementsToRUC
            % Clips line elements against an axis-aligned rectangular RUC box.
            %
            % Inputs
            % ------
            % points : N x 3 array
            %     Original point coordinates.
            %
            % conn : M x 2 array
            %     Connectivity table of line elements. Each row contains two point IDs.
            %
            % rucMin : 1 x 3 array
            %     Minimum corner of the RUC: [xmin ymin zmin].
            %
            % rucDim : 1 x 3 array
            %     RUC dimensions: [Lx Ly Lz].
            %
            % tol : scalar, optional
            %     Numerical tolerance for point comparison.
            %
            % Outputs
            % -------
            % newPts : K x 3 array
            %     New point list after clipping.
            %
            % newConn : P x 2 array
            %     New line-element connectivity after clipping.
            %
            % parentElem : P x 1 array
            %     Original element ID corresponding to each new clipped element.
            %
            % tRange : P x 2 array
            %     Parametric interval [tEnter, tExit] of each clipped line segment.
            %
            % The original line is defined as:
            %     p(t) = p1 + t*(p2 - p1),  0 <= t <= 1
            
            if nargin < 5
                tol = 1e-9;
            end
        
            % Define RUC box
            rucMax = rucMin + rucDim;
        
            % In case dimensions are accidentally given negative
            boxMin = min(rucMin, rucMax);
            boxMax = max(rucMin, rucMax);
        
            % Initialize outputs
            newPts = zeros(0, 3);
            newConn = zeros(0, 2);
            parentElem = zeros(0, 1);
            tRange = zeros(0, 2);
        
            nElem = size(conn, 1);
        
            for e = 1:nElem
        
                % Original point IDs
                id1 = conn(e, 1);
                id2 = conn(e, 2);
        
                % Original coordinates
                p1 = points(id1, :);
                p2 = points(id2, :);
        
                % Clip this line segment against the RUC box
                [keep, q1, q2, tEnter, tExit] = ROI_cutout_toolbox.clipSegmentToBox(p1, p2, boxMin, boxMax, tol);
        
                if ~keep
                    continue;
                end
        
                % Ignore extremely short clipped elements
                if norm(q2 - q1) < tol
                    continue;
                end
        
                % Add/reuse points in the new point list
                newId1 = addUniquePoint(q1);
                newId2 = addUniquePoint(q2);
        
                % Avoid zero-length connectivity
                if newId1 == newId2
                    continue;
                end
        
                % Add new connectivity
                newConn(end+1, :) = [newId1, newId2];
                parentElem(end+1, 1) = e;
                tRange(end+1, :) = [tEnter, tExit];
        
            end
        
            % ----------------------------------------------------------------%
            % Nested function: add point only if it does not already exist
            
            function idx = addUniquePoint(q)
        
                % Snap points very close to RUC faces exactly onto the faces
                q = ROI_cutout_toolbox.snapPointToBox(q, boxMin, boxMax, tol);
        
                if isempty(newPts)
                    newPts(1, :) = q;
                    idx = 1;
                    return;
                end
        
                dist = sqrt(sum((newPts - q).^2, 2));
                existingIdx = find(dist <= tol, 1, 'first');
        
                if isempty(existingIdx)
                    newPts(end+1, :) = q;
                    idx = size(newPts, 1);
                else
                    idx = existingIdx;
                end
            end
        
        end
        
        
        % ----------------------------------------------------------------%
        % Helper function: clip one line segment against an axis-aligned box
        
        function [keep, q1, q2, tEnter, tExit] = clipSegmentToBox(p1, p2, boxMin, boxMax, tol)
        
            d = p2 - p1;
        
            tEnter = 0;
            tExit = 1;
        
            keep = true;
        
            for k = 1:3
        
                if abs(d(k)) < tol
                    % The segment has almost no variation in coordinate k.
                    % Therefore, it is parallel to the two box faces normal to direction k.
                    % If this constant coordinate lies outside the slab, the segment cannot intersect the box.
                    if p1(k) < boxMin(k) - tol || p1(k) > boxMax(k) + tol
                        keep = false;
                        q1 = [];
                        q2 = [];
                        return;
                    end
                else
                    % Compute intersection parameters with lower and upper planes
                    t1 = (boxMin(k) - p1(k)) / d(k);
                    t2 = (boxMax(k) - p1(k)) / d(k);
        
                    tLow = min(t1, t2);
                    tHigh = max(t1, t2);
        
                    % Update valid interval
                    tEnter = max(tEnter, tLow);
                    tExit = min(tExit, tHigh);
        
                    % No overlap means no intersection
                    if tEnter > tExit + tol
                        keep = false;
                        q1 = [];
                        q2 = [];
                        return;
                    end
                end
            end
        
            % Clamp because small numerical errors may push values slightly outside
            tEnter = max(0, min(1, tEnter));
            tExit = max(0, min(1, tExit));
        
            % Compute clipped endpoints
            q1 = p1 + tEnter * d;
            q2 = p1 + tExit * d;
        
            % Snap points close to faces
            q1 = ROI_cutout_toolbox.snapPointToBox(q1, boxMin, boxMax, tol);
            q2 = ROI_cutout_toolbox.snapPointToBox(q2, boxMin, boxMax, tol);
        
        end
        
        
        % ----------------------------------------------------------------%
        % Helper function: snap points very close to box faces
        
        function q = snapPointToBox(q, boxMin, boxMax, tol)
        
            for k = 1:3
                if abs(q(k) - boxMin(k)) <= tol
                    q(k) = boxMin(k);
                elseif abs(q(k) - boxMax(k)) <= tol
                    q(k) = boxMax(k);
                end
        
                % Also clamp tiny numerical overshoots
                if q(k) < boxMin(k) && q(k) > boxMin(k) - tol
                    q(k) = boxMin(k);
                elseif q(k) > boxMax(k) && q(k) < boxMax(k) + tol
                    q(k) = boxMax(k);
                end
            end

        end
        
        % ----------------------------------------------------------------%
        % point ids in RUC faces
        function [front, back, left, right, down, up] = find_boundary_points(box_min, box_max, points, tol)

            front = find(points(:,1) >= box_min(1)-tol & points(:,1) <= box_min(1)+tol);
            back = find(points(:,1) >= box_max(1)-tol & points(:,1) <= box_max(1)+tol);

            left = find(points(:,2) >= box_min(2)-tol & points(:,2) <= box_min(2)+tol);
            right = find(points(:,2) >= box_max(2)-tol & points(:,2) <= box_max(2)+tol);

            down = find(points(:,3) >= box_min(3)-tol & points(:,3) <= box_min(3)+tol);
            up = find(points(:,3) >= box_max(3)-tol & points(:,3) <= box_max(3)+tol);


        end




    end

end