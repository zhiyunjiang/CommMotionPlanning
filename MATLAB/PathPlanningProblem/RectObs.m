classdef RectObs < IObstacle
    %RectObs - A Rectangular obstacle
    %   Detailed explanation goes here
    
    properties
        x_min;
        x_max;
        y_min;
        y_max;
    end
    
    methods (Access = public)
        function obj = RectObs(region)
            %RectObs Construct a RectObs instance
            obj.x_max = region(1);
            obj.x_min = region(2);
            obj.y_max = region(3);
            obj.y_min = region(4);
            
            if (obj.x_max < obj.x_min) || (obj.y_max < obj.y_min)
               error('max values cannot be less than min values') 
            end
        end
        
        function does_obstruct = obstructs(o,a, b)
            %obstructs Checks whether the recangular object obstructs the
            % line segment a-b
            does_obstruct = 0;

            slope = (a(2) - b(2))/(a(1) - b(1));
            if (slope == Inf) || (slope == -Inf)
                %a vertical line
                does_obstruct = (a(1)>=o.x_min && a(1)<=o.x_max) && ...
                                 ( (o.y_min>=min(a(2), b(2)) && o.y_min<=max(a(2), b(2))) || ...
                                   (o.y_max>=min(a(2), b(2)) && o.y_max<=max(a(2), b(2))));   
            else
                intercept = [-slope, 1]*a';
                xinter_ymax = (o.y_max - intercept)/slope;

                
                xinter_ymin = (o.y_min - intercept)/slope;

                yinter_xmin = slope*o.x_min + intercept;
                yinter_xmax = slope*o.x_max + intercept;

                does_obstruct = o.validIntersection([xinter_ymax, o.y_max], a, b) || ...
                                o.validIntersection([xinter_ymin, o.y_min], a, b) || ...
                                o.validIntersection([o.x_min, yinter_xmin], a, b) || ...
                                o.validIntersection([o.x_max, yinter_xmax], a, b);
            end

            %no intersection with the sides, check if we're inside
            if ~does_obstruct
                does_obstruct = (o.containsPoint(a) && o.containsPoint(a));
            end

        end

        function plotObstacle(o, scale, offset)
            if nargin == 1
                scale = 1;
            end
            X = scale*([o.x_min, o.x_min, o.x_max, o.x_max] - offset(1));
            Y = scale*([o.y_min, o.y_max, o.y_max, o.y_min] - offset(2));
            plot(X, Y);
            fill(X, Y, [0.25, 0.25 0.25]);
            
        end
    end

    methods (Access = private)
        function tf = containsPoint(o, point)
            tf = ( (point(1) <= o.x_max) && (point(1) >= o.x_min) && ...
                   (point(2) <= o.y_max) && (point(2) >= o.y_min));
        end

        function tf = validIntersection(o, p, a, b)
            tf = o.containsPoint(p);
            tf = tf && (p(1)>=min(a(1), b(1)) && p(1)<=max(a(1), b(1)));
            tf = tf && (p(2)>=min(a(2), b(2)) && p(2)<=max(a(2), b(2)));
        end
    end
end

