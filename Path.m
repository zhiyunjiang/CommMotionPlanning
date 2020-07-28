classdef Path < handle
    %Path - A path
    
    properties (Access = public)
        points;
        
        pointCount;
        
        steps;
        
    end
    
    methods
        function this = Path(points)
            this.points = points;
            this.pointCount = Path.getPointCount(points);
            this.steps = Path.getSteps(points, this.pointCount);
        end
        
        function l = length(this, start_index, end_index)
            if nargin == 2 || end_index >= this.pointCount %no end_index given or overflows
                %just got total length
                end_index = this.pointCount - 1;
            end
            
            if nargin <= 1 || start_index <= 0 %no start_index given or underflows
                start_index = 1;
            end
            
            if start_index > end_index
                error('start_index must be less than or equal to end_index');
            end
            
            l = sum(this.steps(start_index:end_index));
        end
    end
    
    methods (Static)
        
        function point_count = getPointCount(points)
           point_count = length(points);
           %handle the case of a path consisting of single point
           if min(size(points)) == 1
                point_count = 1;
           end
        end
        
        function steps = getSteps(points, point_counts)
           steps = zeros([point_counts-1, 1]);
           for i = 1:point_counts-1
              steps(i) = norm(points(i,:) - points(i+1,:)); 
           end
        end
    end
end

