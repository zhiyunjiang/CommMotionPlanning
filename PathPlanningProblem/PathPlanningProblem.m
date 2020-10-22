%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PathPlanningProblem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Represents an instance of a path planning problem to be solved. Handles 
% discretization over continuous space as well as obstacle checking and 
% performance metric, as all these should remain the same regardless of how
% the problem is solved.

classdef PathPlanningProblem < matlab.mixin.Copyable
    
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
        
        %ObstacleModule - abstraction of obstacles. Must implement
        %CollisionFree(p1,p2)->{T,F}. Checks to see if the path from p1 to
        %p2 is collision free. May also implement a plot() method
        obstacleMod;
        
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
                l_x = region(1) - region(2); l_y = region(3) - region(4);
                this.offset = [region(2), region(4)];
                this.gridRegion = [l_x*resolution, 0, l_y*resolution, 0];
            end
            
            if length(source) ~= 2
                error('PathPlanningProblem.source must have length of 2 ( [x_start, y_start])');
            else
                this.source = source;
            end
            
            this.goalRegion = GoalRegion(destination_resolution, destinations, this.offset, this.gridRegion);
           
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
        % grd_src - the new source in grid coordinates
        %
        % Output:
        % copy - the copy of the path planning problem
        function copy = copyWithNewSource(this, grd_src)
           %shallow copy using inherited copy() method from matlab.mixin.Copyable
           copy = this.copy(); 
           copy.source = this.toRawCoordinate(grd_src);
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
        function tf = pointInGoalRegion(this, grid_pt)
            if this.resolution == Inf
               %we're dealing with raw points, will need to map to grid first
               tf = this.goalRegion.rawPointInGoalRegion(grid_pt);
            else
                tf = this.goalRegion.gridPointInGoalRegion(grid_pt);
            end
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
                %path will be in terms of grid coordinates only. To
                %compare to objects, convert back to true scale

                p1 = this.toRawCoordinate(path(i,:));
                p2 = this.toRawCoordinate(path(i+1,:));
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
        % gridPtInRegion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cheks if point is in the workspace
        % Input:
        % this - reference to the PathPlanningProblem object
        % grd_pt - point [x,y] in grid coordinates. 
        %
        % Output:
        % tf - true if the point is in the workspace. False otherwise
        function tf = gridPtInRegion(this, grd_pt)
            tf = 0;
            x = grd_pt(1);
            y = grd_pt(2);
            if x <= this.gridRegion(1) && (x >= this.gridRegion(2)) && ...
                    y <= this.gridRegion(3) && (y >= this.gridRegion(4))
               tf = 1; 
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gridPtInRegion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Finds directly reachable neighbors of a grid point
        % Input:
        % this - reference to the PathPlanningProblem object
        % grd_pt - point [x,y] in grid coordinates. 
        %
        % Output:
        % neighbors - neighboring grid points (up to 8). Excludes neighbors
        %               that are not in workspace or for which the
        %               straightline path from grid_pt is obstructed
        function neighbors = getGridNeighbors(this, grd_pt)
           neighbors = [];
           offsets = [[0,1];[1,0];[-1, 0];[0, -1];[1,1];[1, -1];[-1, -1];[-1, 1]];
           for i=1:length(offsets)
              neighbor = grd_pt + offsets(i,:); 
              if this.gridPtInRegion(neighbor)
                  %now check for obstacles
                  sub_path = this.collisionFree([grd_pt; neighbor]);
                  if sum(size(sub_path)) == 4  
                    neighbors = [neighbors; neighbor]; 
                  end
              end
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
           %toGridCoordinate - takes raw point (x,y) and converts to point
           %on grid by (1) shifting by offset, (2) scaling by resolution,
           %and (3) rounding to nearest point
           %
           %Input
           % this - this PathPlanningProblem object
           % val - raw corrdinate (x,y)
           if this.resolution == Inf
               coord = pt - this.offset;
           else
               coord = round((pt - this.offset)/this.getStepSize());
           end
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
        % raw_pt - the point in raw coordinates
        function raw_pt = toRawCoordinate(this, grd_pt)
           if this.resolution == Inf
               raw_pt = grd_pt + this.offset;
           else
               raw_pt = (grd_pt*this.getStepSize()) + this.offset;
           end
        end
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotProb
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots the problem
        % Input:
        % this - reference to the PathPlanningProblem object
        % use_grid_scale - If true, plot with lower left being (0,0), step
        % size of 1. Otherwise, plot with original prolem scale and origin.
        %
        % Output:
        % raw_pt - the point in raw coordinates
        function plotProb(this, use_grid_scale)
           scale = 1;
           plot_offset = [0,0];
           obs_scale = this.resolution;
           obs_offset = this.offset;
           if ~ use_grid_scale
              scale = this.getStepSize(); 
              plot_offset = this.offset;
              obs_scale = 1;
              obs_offset = [0,0];
           end
           
           this.obstacleMod.plotObstacles(obs_scale, obs_offset);
           src = scale * this.getSourceGrid() + plot_offset;
           scatter(src(1), src(2), 90, 'filled', 'ys', 'MarkerEdgeColor', 'k');
           dests = scale*this.goalRegion.goalGridPoints() + plot_offset;
         
           scatter(dests(:,1), dests(:,2), 90, 'filled', 'yd', 'MarkerEdgeColor', 'k');
           
           xlim(scale*[this.gridRegion(2), this.gridRegion(1)] + plot_offset(1));
           ylim(scale*[this.gridRegion(4), this.gridRegion(3)] + plot_offset(2));   
        end
        
        % Getters and Setters
       
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
end

