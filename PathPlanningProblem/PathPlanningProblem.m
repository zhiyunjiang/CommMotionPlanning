classdef PathPlanningProblem < matlab.mixin.Copyable
    %PathPlanningProblem - Represents an instance of a path planning
    %problem to be solved. Handles discretization over continuous space as
    %well as obstacle checking and performance metric, as all these should
    %remain the same regardless of how the problem is solved.
    
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
        
        function copy = copyWithNewSource(this, grd_src)
           %shallow copy using inherited copy() method from matlab.mixin.Copyable
           copy = this.copy(); 
           copy.source = this.toRawCoordinate(grd_src);
        end
        
        function cost = pathCost(this, n1, n2, path, mode)
            %mode = 1 just calculate the distance
            %mode = 2 make any updates to node data (if necessary)
            if nargin == 4
                mode = 1;
            end
            cost = this.costFunction(n1, n2, path, mode);
        end 
        
        function tf = nodeInDestGrid(this, n1)
            pos = n1.getPos();
            tf = this.pointInGoalRegion(pos);
        end
        
        function tf = pointInGoalRegion(this, pt)
            if this.resolution == Inf
               %we're dealing with raw points, will need to map to grid first
               tf = this.goalRegion.rawPointInGoalRegion(pt);
            else
                tf = this.goalRegion.gridPointInGoalRegion(pt);
            end
        end
        
        function cost = Cost(this, n1, n2)
           cost = this.costFunction(n1, n2); 
        end
        
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
        
        function region = getRegion(this)
            region = this.region;
        end
        
        function gridRegion = getGridRegion(this)
           gridRegion = this.gridRegion; 
        end
        
        function tf = gridPtInRegion(this, pt)
            tf = 0;
            x = pt(1);
            y = pt(2);
            if x <= this.gridRegion(1) && (x >= this.gridRegion(2)) && ...
                    y <= this.gridRegion(3) && (y >= this.gridRegion(4))
               tf = 1; 
            end
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
        
        function neighbors = getGridNeighbors(this, grid_loc)
           neighbors = [];
           offsets = [[0,1];[1,0];[-1, 0];[0, -1];[1,1];[1, -1];[-1, -1];[-1, 1]];
           for i=1:length(offsets)
              neighbor = grid_loc + offsets(i,:); 
              if this.gridPtInRegion(neighbor)
                  %now check for obstacles
                  sub_path = this.collisionFree([grid_loc; neighbor]);
                  if sum(size(sub_path)) == 4  
                    neighbors = [neighbors; neighbor]; 
                  end
              end
           end
        end
       
       function coord =  toGridCoordinate(this,val)
           %toGridCoordinate - takes raw point (x,y) and converts to point
           %on grid by (1) shifting by offset, (2) scaling by resolution,
           %and (3) rounding to nearest point
           %
           %Input
           % this - this PathPlanningProblem object
           % val - raw corrdinate (x,y)
           if this.resolution == Inf
               coord = val - this.offset;
           else
               coord = round((val - this.offset)/this.getStepSize());
           end
       end
       
       function raw = toRawCoordinate(this, grd_val)
           if this.resolution == Inf
               raw = grd_val + this.offset;
           else
               raw = (grd_val*this.getStepSize()) + this.offset;
           end
       end
       
       function plotProb(this, use_grid_scale)
           scale = 1;
           offset = [0,0];
           obs_scale = this.resolution;
           obs_offset = this.offset;
           if ~ use_grid_scale
              scale = this.getStepSize(); 
              offset = this.offset;
              obs_scale = 1;
              obs_offset = [0,0];
           end
           
           this.obstacleMod.plotObstacles(obs_scale, obs_offset);
           src = scale * this.getSourceGrid() + offset;
           plot(src(1), src(2), 'gx');
           dests = scale*this.goalRegion.goalGridPoints() + offset;
         
           scatter(dests(:,1), dests(:,2), 'ko');
           
           xlim(scale*[this.gridRegion(2), this.gridRegion(1)] + offset(1));
           ylim(scale*[this.gridRegion(4), this.gridRegion(3)] + offset(2));   
       end
    end
end

