%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GoalRegion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Represents the destination of our path planning problem. The
% destinatation region is taken to be a set of points from the simulated
% communication channel, and will thus generally has the same resolution as
% the communication channel simulation. This may differ from the problem
% resolution.

classdef GoalRegion < handle
    
    % Properties (public):
    % resolution - the resolution used for enumerating the set of poitns in
    %               the goal region.
    % offset - the offset of points in terms of the problem instance
    %           resolution
    % problemGridRegion - [x_max, 0, y_max, 0] the problem instance region
    %                       with max dimensions in terms of problem
    %                       resolution
    % goalPoints - list of points (in problem resolution) within goal
    %               region, pruned so that only points on border are inluded.
    % goalGirdPoints - list of points (in  resolution) within goal
    %               region, pruned so that only points on border are inluded.
    % goalGridMap - map of maps. First map key is x corrdinate, second map
    %               keys are y values. Use for quickly looking up whether 
    %               or not a given point is in the goal region 
    properties
        resolution;
        offset;
        problemGridRegion;
        
        goalPoints;
        
        goalGridPoints;
        
        goalGridMap;
    end
    
    % Methods (public)
    % (Constructor) - create new GoalRegion object
    % rawPointInGoalRegion - check is a point not necessarily on the grid
    %                           is in the goal region
    % gridPointInGoalRegion - check if a point on the grid is in the goal
    %                           region
    methods
        function this = GoalRegion(resolution, goal_points, offset, grid_region)
            if resolution == Inf
               error('Goal region must have finite resolution'); 
            end
            this.resolution = resolution;
            this.offset = offset;
            this.problemGridRegion = grid_region;
            this.setGoals(goal_points);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % rawPointInGoalRegion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check is a point not necessarily on the grid is in the goal region
        % Input:
        % this - reference to the GoalRegion object
        % pt - the point [x,y] we're interested in
        %
        % Output:
        % tf - true if point is in goal region. False otherwise.
        function tf =  rawPointInGoalRegion(this, pt)
            grid_pt = this.toGridCoordinate(pt);
            tf = this.gridPointInGoalRegion(grid_pt);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gridPointInGoalRegion
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % check if a point on the grid is in the goal region
        % Input:
        % this - reference to the GoalRegion object
        % pt - the point [x,y] we're interested in
        %
        % Output:
        % tf - true if point is in goal region. False otherwise.
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
                 %submap(y) is a handle object - modification persists 
                 %outside this function call
                 submap(y) = i;
              else
                  submap = containers.Map(y,i);
                  map(x) = submap;
              end
           end
        end
       
        function kept_indices = pruneDestinations(this)
            
            kept_indices = [];
            %to avoid duplicates, actually want to iterate over the values
            %in the map
            x_keys = this.goalGridMap.keys;
            x_count = length(x_keys);
            offsets = [[0, 1];[0, -1];[1, 0];[-1, 0];...
                        [1,1];[-1,1];[1,-1];[-1,-1]];
            
            for i=1:x_count
                x = cell2mat(x_keys(i));
                y_map = this.goalGridMap(x);
                y_keys = y_map.keys;
                y_count = length(y_keys);
                
                for j = 1:y_count
                    y = cell2mat(y_keys(j));
                    
                    dest = [x ,y]; 
                    keep = 0;
                    for k = 1:length(offsets)
                        neighbor = dest + offsets(k,:);
                        if  (~this.gridPointInGoalRegion(neighbor) ...
                                && this.inProblemGridRegion(neighbor))
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
        
        function tf = inProblemGridRegion(this, pt)
            tf = 0;
            x = pt(1);
            y = pt(2);
            if x <= this.problemGridRegion(1) && (x >= this.problemGridRegion(2)) && ...
                    y <= this.problemGridRegion(3) && (y >= this.problemGridRegion(4))
               tf = 1; 
            end
        end
       
       function coord =  toGridCoordinate(this,val)
           if this.resolution == Inf
               coord = val - this.offset;
           else
                step_size = 1/this.resolution;
                coord = round((val - this.offset)/step_size);
           end
       end
    end
end

