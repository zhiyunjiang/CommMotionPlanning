classdef PathPlanningProblem < handle
    %PathPlanningProblem - Represents an instance of a path planning
    %problem to be solved. Handles discretization over continuous space as
    %well as obstacle checking as performance metric, as all these should
    %remain the same regardless of how the problem is solved.
    
    properties (Access = private)
        %region = [x_max x_min y_max y_min];
        region;
        
        %gridRegion = [grid_x_max grid_x_min grid_y_max grid_y_min] 
        gridRegion;
        
        %resolution - # of tick marks between integer values in the
        %discretized space. Can also be thought of as 1/step size
        resolution;
        
        %source = [x_start, y_start]
        source;
        
        %secondary source
        source2;
        
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
                this.gridRegion = this.toGridCoordinate(region);
            end
            
            if length(source) ~= 2
                error('PathPlanningProblem.source must have length of 2 ( [x_start, y_start])');
            else
                this.source = source;
            end
            
            this.goalRegion = GoalRegion(destination_resolution, destinations);
           
            this.obstacleMod = obstacle_mod;
            this.costFunction = cost_func;
        end
        
        function cost = pathCost(this, n1, n2, path)
            cost = this.costFunction(n1, n2, path);
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
            was_obstructed = 0;
            [lngth, ~] = size(path);
            
            for i = 1:lngth - 1
                %path will be in terms of grid coordinates only. To
                %compare to objects, convert back to true scale
                scale = this.getStepSize();
                if scale == 0
                    %this is a continuous space problem, just use raw
                    %points
                    scale = 1;
                end
                p1 = scale*path(i,:);
                p2 = scale*path(i+1,:);
                if ~ this.obstacleMod.collisionFree(p1, p2)
                    was_obstructed = 1;
                   break; 
                end
            end
            if was_obstructed
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
        
        function regionDscrt = getRegionDscrt(this)
           regionDscrt = this.roundToResolution(this.region); 
        end
        
        function step = getStepSize(this)
           step = 1/this.resolution; 
        end
        
        function source = getSource(this)
           source = this.source; 
        end
        
        
        function dscrt_source = getSourceDscrt(this)
            dscrt_source = this.roundToResolution(this.getSource());
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
        
        
        function dscrt_source2 = getSourceDscrt2(this)
            dscrt_source2 = this.roundToResolution(this.getSource2());
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
           offsets = [[0,1];[1,0];[-1, 0];[0, -1]];
           for i=1:4
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
        
        %rounding so that 0,0 always on grid, that is, so that grid aligns 
        %with integer values
       function rounded = roundToResolution(this, point)
           if this.resolution == Inf
               rounded = point;
           else
               rounded = this.toGridCoordinate(point)/this.resolution;
           end
           
           
       end 
       
       function coord =  toGridCoordinate(this,val)
           if this.resolution == Inf
               coord = val;
           else
               step_size = 1/this.resolution;
           coord = round(val/step_size);
           end
       end
       
       
       function plotProb(this, use_grid_scale)
           scale = 1;
           obs_scale = 1/this.getStepSize(); 
           if ~ use_grid_scale
              scale = this.getStepSize(); 
              obs_scale = 1;
           end
           
           this.obstacleMod.plotObstacles(obs_scale);
           src = scale * this.getSourceGrid();
           plot(src(1), src(2), 'gx');
           dests = scale*this.goalRegion.goalGridPoints();
         
           scatter(dests(:,1), dests(:,2), 'ko');
           
           xlim(scale*[this.gridRegion(2), this.gridRegion(1)]);
           ylim(scale*[this.gridRegion(4), this.gridRegion(3)]);
           
       end
    end
end

