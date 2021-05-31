classdef PointObstacle < IObstacle
    %POINTOBSTACLE
    
    properties
        map
    end
    
    % Methods(public)
    % (Constructor) - creates a new CircObs obstacle
    % obstructs - checks if the this CircObs obstructs the path between
    %                   two waypoints
    % plotObstacle - plots the obstacle
    methods (Access = public)
        function this = PointObstacle(point_table)
            this.map = point_table;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % obstructs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % checks if this set of points obstructs the path between two
        % waypoints.
        % Input:
        % this - reference to the PointObstacle object
        % a - first waypoint [x,y]
        % b - second waypoint [x,y]
        %
        % Output:
        % does_obstruct - true if obstacle blocks path. False otherwise.
        function does_obstruct = obstructs(this,a, b)
            does_obstruct = 0;
            if this.map.pointInTable(a) || this.map.pointInTable(b)
                does_obstruct = 1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotObstacle
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plots the obstacle
        % Input:
        % this - reference to the PointObstacle object
        % scale - scale on which plotting will occur relative to obstacle's
        %           scale
        % offset - how much to shift frame of reference. Offset in terms of
        %           scale.
        function plotObstacle(this, scale, offset)
            if nargin < 3
                offset = [0,0];
                if nargin < 2
                    scale = 1;
                end
            end
            
            %plot((scale*X)+offset(1),  (scale*Y)+offset(2))
            warning('Plotting of point obstacles not yet supported');
        end
    end
    
end

