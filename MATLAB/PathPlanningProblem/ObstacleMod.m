%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ObstacleMod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Containing class for the collection of obstacles the workspace of our
% path palnning problem

classdef ObstacleMod < handle
    % Properties (public):
    % obstacles - list of all obstacles
    properties
        obstacles
    end
    
    % Methods (public):
    % (constructor) - create new ObstacleMod object
    % collisionFree - checks if the straightline path form p1 to p2 is
    %                   obstacle free
    % plotObstacles - plots all obstacles
    methods
        function this = ObstacleMod(obstacles)
            this.obstacles = obstacles;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % collisionFree
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Checks if the straightline path from p1 to p2 is obstacle free
        %
        % Inputs:
        % this - reference to this ObstacleMod object
        % p1 - first waypoint
        % p2 - second waypoiny
        %
        % Output:
        % tf - true if the path is collision free, false otherwise
        function tf = collisionFree(this, p1, p2)
            tf = 1;
            for i=1:length(this.obstacles)
               obstacle = this.obstacles(i);
               tf = ~obstacle.obstructs(p1, p2); 
               if ~tf
                   break
               end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotObstacles
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots all obstacles
        %
        % Inputs:
        % this - reference to this ObstacleMod object
        % scale - scale to use when plotting obstacles
        % offset -offset to use when plotting obstacles. Offset in terms of
        %           problem resolution
        function plotObstacles(this, scale, offset)
           for i=1:length(this.obstacles)              
              this.obstacles(i).plotObstacle(scale, offset);
              hold on
           end
        end
    end
end

