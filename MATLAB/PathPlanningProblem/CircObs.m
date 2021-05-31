%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CircObs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements the IObstacle interface. Represents a circular obstacle.

classdef CircObs < IObstacle
    
    % Properties (public):
    % r - radius
    % c = circle center
    properties
        %radius
        r = 1;
        %center
        c = [0,0];
    end
    
    % Methods(public)
    % (Constructor) - creates a new CircObs obstacle
    % obstructs - checks if the this CircObs obstructs the path between
    %                   two waypoints
    % plotObstacle - plots the obstacle
    methods (Access = public)
        
        function obj = CircObs(radius, center)
           obj.r = radius;
           obj.c = center;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % obstructs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % checks if the this CircObs obstructs the path between two
        % waypoints.
        % Input:
        % this - reference to the CircObs object
        % a - first waypoint [x,y]
        % b - second waypoint [x,y]
        %
        % Output:
        % does_obstruct - true if obstacle blocks path. False otherwise.
        function does_obstruct = obstructs(this,a, b)
            does_obstruct = 1;
            %first, calculate the minimum distance to the line on which
            %line segment a-b lies.
            a_to_c = this.c - a;
            a_to_b = b - a;
            %find angle between ac and ab
            
            theta = real(acos(max(min(dot(a_to_c,a_to_b)/(norm(a_to_c)*norm(a_to_b)),1),-1)));
            
            %find mindist from center to line on which line segment lies
            if ~isreal(a_to_c) || ~isreal(theta)
                bad = 1;
            end
            min_dist = norm(a_to_c)*sin(mod(theta, pi));
            
            %if the min distance from the circle center to the line is
            %greater than the circle radius, we're set
            if min_dist > this.r
                does_obstruct = 0;
            else
                %oherwise, check if min occurs within the line segment. 
                theta_c = atan2(a_to_c(2), a_to_c(1));
                x_delta = abs(min_dist*cos(theta_c)) * sign(a(1)-this.c(1));
                x_min = x_delta + this.c(1);
                
                if (a(1)>x_min && b(1)>x_min) || (a(1)<x_min && b(1)<x_min)
                    %if not, then just check the end points
                    does_obstruct = ((norm(a_to_c) <= this.r) || (norm(this.c-b) <= this.r));
                end
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotObstacle
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plots the obstacle
        % Input:
        % this - reference to the CircObs object
        % scale - scale on which plotting will occur relative to obstacle's
        %           scale
        % offset - how much to shift frame of reference. Offset in terms of
        %           scale.
        function plotObstacle(this, scale, offset)
            th = 0:pi/50:2*pi;
            rad = scale*this.r;
            cent = scale*(this.c) - offset;
            X = rad*cos(th) + cent(1);
            Y = rad*sin(th) + cent(2);
            plot(X, Y)
            fill(X, Y, [0.25, 0.25 0.25]);
        end
        
    end
end

