classdef CircObs < IObstacle
    %CircObs A circular obstacle
    
    properties
        %radius
        r = 1;
        %center
        c = [0,0];
    end
    
    methods (Access = public)
        
        function obj = CircObs(radius, center)
           obj.r = radius;
           obj.c = center;
        end
        
        function does_obstruct = obstructs(obj,a, b)
            does_obstruct = 1;
            %first, calculate the minimum distance to the line on which
            %line segment a-b lies.
            a_to_c = obj.c - a;
            a_to_b = b - a;
            %find angle between ac and ab
            theta = acos(a_to_c*a_to_b'/(norm(a_to_c)*norm(a_to_b)));
            
            %find mindist from center to line on which line segment lies
            min_dist = norm(a_to_c)*sin(mod(theta, pi));
            
            %if the min distance from the circle center to the line is
            %greater than the circle radius, we're set
            if min_dist > obj.r
                does_obstruct = 0;
            else
                %oherwise, check if min occurs within the line segment. 
                theta_c = atan2(a_to_c(2), a_to_c(1));
                x_delta = abs(min_dist*cos(theta_c)) * sign(a(1)-obj.c(1));
                x_min = x_delta + obj.c(1);
                
                if (a(1)>x_min && b(1)>x_min) || (a(1)<x_min && b(1)<x_min)
                    %if not, then just check the end points
                    does_obstruct = ((norm(a_to_c) <= obj.r) || (norm(obj.c-b) <= obj.r));
                end
            end
        end
        
        function plotObstacle(o, scale, offset)
            th = 0:pi/50:2*pi;
            rad = scale*o.r;
            cent = scale*(o.c - offset);
            X = rad*cos(th) + cent(1);
            Y = rad*sin(th) + cent(2);
            plot(X, Y)
            fill(X, Y, 'g');
        end
        
    end
end

