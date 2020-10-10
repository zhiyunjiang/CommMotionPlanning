classdef ObstacleMod < handle
    properties
        obstacles
    end
    
    methods
        function o = ObstacleMod(obstacles)
            o.obstacles = obstacles;
        end
        
        function tf = collisionFree(o, p1, p2)
            tf = 1;
            for i=1:length(o.obstacles)
               obstacle = o.obstacles(i);
               tf = ~obstacle.obstructs(p1, p2); 
               if ~tf
                   break
               end
            end
        end
        
        function plotObstacles(o, scale, offset)
           for i=1:length(o.obstacles)              
              o.obstacles(i).plotObstacle(scale, offset);
              hold on
           end
        end
    end
end

