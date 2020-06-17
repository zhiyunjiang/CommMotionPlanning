classdef GoalRegion < handle
    %Grid - all coordinates are integer values
    %Rounded - all coordinates are multiples of the step size
    %(1/resolution)
    %Raw - the raw coordinates
    
    properties
        resolution;
        
        goalPoints;
        
        goalGridPoints;
        
        goalGridMap;
        
    end
    
    methods
        function this = GoalRegion(resolution, goal_points)
            if resolution == Inf
               error('Goal region must have finite resolution'); 
            end
            this.resolution = resolution;
            this.setGoals(goal_points);
        end
        
        function tf =  rawPointInGoalRegion(this, pt)
            grid_pt = this.toGridCoordinate(pt);
            tf = this.gridPointInGoalRegion(grid_pt);
        end
        
        function tf = gridPointInGoalRegion(this, pt)
            tf = 0;
            x = pt(1);

            if this.goalGridMap.isKey(x)
                y = pt(2);
                submap = this.goalGridMap(x);
                if submap.isKey(y)
                   tf = 1; 
                end
            end
        end
        
    end
    
    methods(Access = private)
       function setGoals(this, goal_points)
            
            this.goalGridPoints = this.toGridCoordinate(goal_points);
            
            this.goalGridMap = this.createDestinationsMap(); 
            
            kept_indices = this.pruneDestinations();
            
            this.goalPoints = goal_points(kept_indices,:);
            this.goalGridPoints = this.goalGridPoints(kept_indices,:);
            %consider keep full map so that if we land within a connected region, we
            %still recognize this as arriving.
            this.goalGridMap = this.createDestinationsMap(); 
        end
       
        function map = createDestinationsMap(this)
           map = containers.Map('KeyType','double','ValueType','any');
           
           [dest_count, ~] = size(this.goalGridPoints);

           for i=1:dest_count
              x = this.goalGridPoints(i,1);
              y = this.goalGridPoints(i,2);
              if map.isKey(x)
                 submap = map(x);
                 submap(y) = i;
              else
                  submap = containers.Map(y,i);
                  map([x]) = submap;
              end
           end
        end

       
        function kept_indices = pruneDestinations(this)
            
            kept_indices = [];
            %to avoid duplicates, actually want to iterate over the values
            %in the map
            x_keys = this.goalGridMap.keys;
            x_count = length(x_keys);
            offests = [[0, 1];[0, -1];[1, 0];[-1, 0]];
            
            for i=1:x_count
                x = cell2mat(x_keys(i));
                y_map = this.goalGridMap(x);
                y_keys = y_map.keys;
                y_count = length(y_keys);
                
                for j = 1:y_count
                    y = cell2mat(y_keys(j));
                    
                    dest = [x ,y]; 
                    keep = 0;
                    for k = 1:4
                        neighbor = dest + offests(k,:);
                        if  ~this.gridPointInGoalRegion(neighbor)
                            keep = 1;
                            break;
                        end
                    end

                    if keep
                        index_val = y_map(y);
                        kept_indices = [kept_indices, index_val];
                    end
                end
            end
            
        end
        
        % rounding so that 0,0 always on grid, that is, so that grid aligns 
        % with integer values
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
    end
end

