%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PathPlanningProblem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Represents an instance of a path planning problem to be solved. Handles 
% discretization over continuous space as well as obstacle checking and 
% performance metric, as all these should remain the same regardless of how
% the problem is solved.

classdef PathPlanningProblem < matlab.mixin.Copyable
    
    properties (Access = public)
        %ObstacleModule - abstraction of obstacles. Must implement
        %CollisionFree(p1,p2)->{T,F}. Checks to see if the path from p1 to
        %p2 is collision free. May also implement a plot() method
        obstacleMod;  
    end
    
    
    properties (Access = private)
        %region = [x_max x_min y_max y_min]; Raw coordinates from problem
        region;
        
        %offset - [x_offset, y_offest]; Same as [x_min, y_min];
        offset;
        
        %gridRegion = [grid_x_max grid_x_min grid_y_max grid_y_min];
        %corrdinates divded by step size,may not result in integer
        % corrdinates. x and y min will alway be 0
        gridRegion;
        
        %resolution - # of tick marks between integer values in the
        %discretized space. Can also be thought of as 1/step size
        resolution;
        
        %source = [x_start, y_start]; Raw coordinates
        source;
        
        %secondary source = [x_start, y_start]; Used only when solving
        %mutli source problem (i.e. arrive at source 2 from source (1)
        %after passing through a connected point).
        source2;
        
        %goalRegion - GoalRegion object. The set of points to which the
        %robot is heading. 
        goalRegion;
                
        %costFunction - maps two points (n1, n2) to the cost between
        %them. n1, n2 are TreeNode objects.
        costFunction;
    end
    
    % Methods (public):
    % (constructor) - create new PathPlanningProblem instance
    % copyWithNewSource - create a copy of the PathPlanningProblem instance
    % pathCost - evaluates the cost of a path or path segment
    % nodeInDestGrid - checks if a node (class TreeNode) is in the goal
    %           region. Assumes the nodes position is in grid coordinates
    % pointInGoalRegion - checks if a point is in the goal region
    % collisionFree - finds portion of path that is collision free
    % gridPtInRegion - checks if a given point (in grid coordinates) is in
    %                   the workspace
    
    % getGridNeighbors - finds grid coordinates for all neighbors of a
    %                    given coordinate, including diagonals
    % toGridCoordinate - takes a raw poit and translates to nearest grid
    %                    coordinate.
    % toRawCoodinate - takes a grid coordinate and converts to raw
    %                   coordinate
    % plotProb - plots the obstacles, goal region, source, and destination
    % getRegion - returns the max and min values of the coordiante axis of
    %              the region as [x max, x min, y max, y min]
    % getGridRegion - returns max and min values of the coordinates axis
    %                   for the region after translating x min, y min to 0
    %                   and scaling so the step size in the grid equals 1
    % getStepSize - returns the step size (inverse of resolution)
    % getSource - returns the source/starting point
    % getSourceGrid - returns the source in grid coordinates
    % setSource2 - sets a secondary source. Used in e.g. double tree RRT
    % getSource2 - returns the secondary source if one has been provided.
    %               If no secondary source for this problem, returns []
    % getSourceGrid2 - returns the secondary source in grid coordinates 
    % getGoalRegion - returns the goal region object
    % getObstacleMod - returns the obstacle module object
    methods
        
        function this = PathPlanningProblem(region, resolution, source, destinations, obstacle_mod, cost_func, destination_resolution)
            this.resolution = resolution;
            
            %destination resolution not supplied
           if resolution == Inf
               if nargin ~=7
                   error('Must supply a destination resolution if overall problem is continuous'); 
               end
           else
               destination_resolution = resolution;
           end
              
            
            if length(region) ~= 4
               error('PathPlanningProblem.region must have length of 4 ( [x_max x_min y_max y_min])');
            else
                this.region = region;
                this.offset = [region(2), region(4)];
                this.gridRegion = reshape(this.toGridCoordinate(reshape(region, [2,2])),size(region));
            end
            
            if length(source) ~= 2
                error('PathPlanningProblem.source must have length of 2 ( [x_start, y_start])');
            else
                this.source = source;
            end
            
            this.goalRegion = GoalRegion(destination_resolution, destinations, this.region);
           
            this.obstacleMod = obstacle_mod;
            this.costFunction = cost_func;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % copyWithNewSource
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Creates an identical copy of the PathPlanningProblem but with a
        % different source
        % Input:
        % this - reference to the PathPlanningProblem object
        % new_src - the new source
        %
        % Output:
        % copy - the copy of the path planning problem
        function copy = copyWithNewSource(this, new_src)
           %shallow copy using inherited copy() method from matlab.mixin.Copyable
           copy = this.copy(); 
           copy.source = new_src;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % copyWithNewSource
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evalues the path bewteen n1 and n2 using the problem's cost
        % function.
        % Input:
        % this - reference to the PathPlanningProblem object
        % n1 - the starting node (TreeNode)
        % n2 - the ending node (TreeNode)
        % path  - waypoints between node 1 and 2
        % mode - optionally have cost function perform additional work for
        %           values of mode
        % Output:
        % cost - the cost associated with moving from n1 to n1
        function cost = pathCost(this, n1, n2, path, mode)
            %mode = 1 just calculate the distance
            %mode = 2 make any updates to node data (if necessary)
            if nargin == 4
                mode = 1;
            end

            cost = this.costFunction(n1, n2, path, mode);
        end 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % nodeInDestGrid
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cheks if node is in the goal region
        % Input:
        % this - reference to the PathPlanningProblem object
        % n1 - node (TreeNode)
        %
        % Output:
        % tf - true if the node is in the goal region. False otherwise
        function tf = nodeInDestGrid(this, n1)
            pos = n1.getPos();
            tf = this.pointInGoalRegion(pos);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % pointInGoalRegion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cheks if point is in the goal region
        % Input:
        % this - reference to the PathPlanningProblem object
        % grid_pt - point [x,y] in grid coordinates. If the problem has not
        %           been discretized, grid point will be a raw point
        %
        % Output:
        % tf - true if the point is in the goal region. False otherwise
        function tf = pointInGoalRegion(this, pt)
            tf = this.goalRegion.pointInGoalRegion(pt);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % collisionFree
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cheks if node is in the goal region
        % function.
        % Input:
        % this - reference to the PathPlanningProblem object
        % path - series of path waypoints
        %
        % Output:
        % viable_subpath - returns the portion of the path that is not
        %                   obstructed. For example if there is a collision
        %                   at path(i), will return path(1:i-1)
        function viable_subpath = collisionFree(this, path)
            terminate_early = 0;
            [lngth, ~] = size(path);
            
            for i = 1:lngth - 1
                %check if the point is in the goal region
                if this.pointInGoalRegion(path(i,:))
                    terminate_early = 1;
                    break
                end
                p1 = path(i,:);
                p2 = path(i+1,:);
                if ~ this.obstacleMod.collisionFree(p1, p2)
                    terminate_early = 1;
                   break; 
                end
            end
            if terminate_early
                viable_subpath = path(1:i,:);
            else
                viable_subpath = path;
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ptInRegion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cheks if point is in the workspace
        % Input:
        % this - reference to the PathPlanningProblem object
        % pt - point [x,y] 
        %
        % Output:
        % tf - true if the point is in the workspace. False otherwise
        function tf = ptInRegion(this, pt)
            tf = 0;
            x = pt(1);
            y = pt(2);
            if x <= this.region(1) && (x >= this.region(2)) && ...
                    y <= this.region(3) && (y >= this.region(4))
               tf = 1; 
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % getGridNeighbors
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Finds directly reachable neighbors of a point
        % Input:
        % this - reference to the PathPlanningProblem object
        % pt - point [x,y] in workspace coordinates. 
        %
        % Output:
        % neighbors - neighboring grid points (up to 8). Excludes neighbors
        %               that are not in workspace or for which the
        %               straightline path from grid_pt is obstructed
        function neighbors = getNeighbors(this, pt)
           grd_pt = this.toGridCoordinate(pt);
           neighbors = [];
           offsets = [[0,1];[1,0];[-1, 0];[0, -1];[1,1];[1, -1];[-1, -1];[-1, 1]];
           for i=1:length(offsets)
              neighbor = grd_pt + offsets(i,:); 
              if this.ptInRegion(this.toRawCoordinate(neighbor))
                  %now check for obstacles
                  sub_path = this.collisionFree(this.toRawCoordinate([grd_pt; neighbor]));
                  if sum(size(sub_path)) == 4  
                    neighbors = [neighbors; neighbor]; 
                  end
              end
           end
           if ~isempty(neighbors)
               neighbors = this.toRawCoordinate(neighbors);
           end
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % toGridCoordinate
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Converts a raw point (arbitrary offset, step size) to a
        % grid point ( (0,0) min, step size = 1). Rounds when necessary to
        % snap to grid.
        % Input:
        % this - reference to the PathPlanningProblem object
        % pt - point [x,y] in raw coordinates. 
        %
        % Output:
        % coord - grid coordinates of pt
        function coord =  toGridCoordinate(this,pt)
            
           coord = toGridFromRaw(this.region, this.resolution, pt);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % toRawCoordinate
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Converts a grid point ( (0,0) min, step size = 1) to a raw point
        % (arbitrary offset, step size = 1/res).
        % Input:
        % this - reference to the PathPlanningProblem object
        % grd_pt - point [x,y] in grid coordinates. 
        %
        % Output:
        % coord - the point in raw coordinates
        function coord = toRawCoordinate(this, grd_pt)
           coord = toRawFromGrid(this.region, this.resolution, grd_pt);
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotProb
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots the problem
        % Input:
        % this - reference to the PathPlanningProblem object
        %
        function plotProb(this)
           obs_scale = 1;
           obs_offset = [0,0];

           
           this.obstacleMod.plotObstacles(obs_scale, obs_offset);
           src = this.getSource();
           scatter(src(1), src(2), 90, 'filled', 'ys', 'MarkerEdgeColor', 'k');
           dests = this.goalRegion.goalPoints();
         
           scatter(dests(:,1), dests(:,2), 90, 'filled', 'yd', 'MarkerEdgeColor', 'k');
           
           xlim([this.region(2), this.region(1)]);
           ylim([this.region(4), this.region(3)]);   
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % isContinuous
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Indicates whether or not the problem has been posed as continuous
        % or discrete.
        % Input:
        % this - reference to the PathPlanningProblem object
        %
        % Output:
        % tf - 1 if problem is in continuous space, 0 if discretized
        function tf = isContinuous(this)
            tf = isinf(this.resolution);
        end
        
        function path = getSLPath(this, start, dest)
            start = this.toGridCoordinate(start);
            dest = this.toGridCoordinate(dest);
            if this.isContinuous()
                % may still need break up into waypoints
                path = [start; dest];
            else
                %approximate a straightline through the grid
                dist = start - dest;
                dist_abs = abs(dist);
                mabs = dist_abs(2)/dist_abs(1);
                min_path_length = sum(dist_abs) + 1;

                path = zeros([min_path_length, 2]);
                path(1,:) = start;
                last_path_index = min_path_length;
                for i = 2:min_path_length
                    %just straight up /down
                    if mabs == Inf
                       path(i,:) = path(i-1,:) + [0,-1*sign(dist(2))];
                    elseif mabs == 0
                        %just straight left or right
                        path(i,:) = path(i-1,:) + [-1*sign(dist(1)), 0];
                    else
                        %get slope from v_nearest to previous point.
                        prev_tot_delta = abs(path(1,:) - path(i-1, :));
                        delta = PathPlanningProblem.getGridSLPathDelta(mabs, dist, prev_tot_delta(1), prev_tot_delta(2));
                        path(i,:) = path(i-1,:) + delta;
                    end

                    if all(path(i,:) == dest)
                        last_path_index = i;
                        break;
                    end
                end

                path = path(1:last_path_index,:);

            end
            path = this.toRawCoordinate(path);
         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Getters and Setters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function region = getRegion(this)
            region = this.region;
        end
        
        function gridRegion = getGridRegion(this)
           gridRegion = this.gridRegion; 
        end
        
        function step = getStepSize(this)
           step = 1/this.resolution; 
        end
        
        function source = getSource(this)
           source = this.source; 
        end
        
        function grd_source = getSourceGrid(this)
           grd_source = this.toGridCoordinate(this.getSource()); 
        end
        
        function setSource2(this, src_2)
           this.source2 = src_2; 
        end
        
        function source_2 = getSource2(this)
           source_2 = this.source2; 
        end
        
        function grd_source_2 = getSourceGrid2(this)
           grd_source_2 = this.toGridCoordinate(this.getSource2()); 
        end
        
        function goal_region =  getGoalRegion(this)
           goal_region = this.goalRegion; 
        end
        
        function obstacle_mod = getObstacleMod(this)
           obstacle_mod = this.obstacleMod; 
        end
        
    end
    
    methods (Access = private, Static = true)
         function delta = getGridSLPathDelta(true_slope_abs, dist, cur_dx_abs, cur_dy_abs)
             next_slopes_abs = (cur_dy_abs + [1, 0, 1])./(cur_dx_abs + [0, 1, 1]);
             [~,dir] =min(abs(true_slope_abs - next_slopes_abs)); 
             delta = [0,0];
             if dir == 1 || dir == 3
               %need more rise
               delta =  [0, -1*sign(dist(2))];
             end
             if dir >= 2
                %needs to run more
               delta = delta + [-1*sign(dist(1)), 0];
             end
         end
    end
end

