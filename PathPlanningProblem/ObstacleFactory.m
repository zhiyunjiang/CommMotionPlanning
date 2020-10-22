%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ObstacleFactory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Collection of static methods for quickly generating different groups of
% obstacles. These all assum a 50 x 50 workspace. All methods return an
% ObsMod object.

classdef ObstacleFactory
    
    % Methods (static):
    % maze - create a zig zag type set of obstacles 
    % simple - three circular obstacles placed around vertical line at x=35
    % noObs - creates an ObstacleMod object with no obstacles
    % random - randomly places 6 circular objects with varying radii in
    %           workspace
    % cirlceGrid - places 7 circle obstacles in a grid, but with large
    %               jitter and random radii
    % custom - for ad hoc obstacle placement
    methods (Static)
        
        function obs_mod = maze()
            maze = [RectObs([7, 0, 35, 15]), RectObs([23, 16, 50, 27]), ...
            RectObs([38, 30, 30, 0]), RectObs([50, 42, 35, 25])];
            
            obs_mod = ObstacleMod(maze);
        end
        
        function obs_mod = simple()
            simple = [CircObs(3,[34,10]), CircObs(3,[37,44]), CircObs(3,[34,30])]; 
            obs_mod = ObstacleMod(simple);
        end
        
        function obs_mod = noObs()
            obs_mod = ObstacleMod([]);
        end
        
        function obs_mod = random(start, dest)
            valid = 0;
            while ~valid
               rand_obs = [CircObs(randi(3)+2, [randi(51)-1, randi(51)-1]), CircObs(randi(3)+2, [randi(51)-1, randi(51)-1]),...
                    CircObs(randi(3)+2, [randi(51)-1, randi(51)-1]), CircObs(randi(3)+2, [randi(51)-1, randi(51)-1]),...
                    CircObs(randi(3)+2, [randi(51)-1, randi(51)-1]), CircObs(randi(3)+2, [randi(51)-1, randi(51)-1])];

                obs_mod = ObstacleMod(rand_obs);
                if obs_mod.collisionFree(start, start) && obs_mod.collisionFree(dest, dest)
                    valid = 1;
                end
            end
        end
        
        function obs_mod = circleGrid(r_base, r_plus_max, pos_rng, start, dest)
            valid = 0;
            while ~valid
                circle_grid = [CircObs(randi(r_plus_max)+r_base, [10,40] + randi(pos_rng,[1,2])), CircObs(randi(r_plus_max)+r_base, [25, 10] + randi(pos_rng,[1,2])),...
                    CircObs(randi(r_plus_max)+r_base, [25,25] + randi(pos_rng,[1,2])), CircObs(randi(r_plus_max)+r_base, [25, 40] + randi(pos_rng,[1,2])),...
                    CircObs(randi(r_plus_max)+r_base, [40,10] + randi(pos_rng,[1,2])), CircObs(randi(r_plus_max)+r_base, [40, 25] + randi(pos_rng,[1,2])),...
                    CircObs(randi(r_plus_max)+r_base, [40,40] + randi(pos_rng,[1,2]))];

                obs_mod = ObstacleMod(circle_grid);
                if obs_mod.collisionFree(start, start) && obs_mod.collisionFree(dest, dest)
                    valid = 1;
                end
            end
        end
        
        function obs_mod = custom()
            obstacles = [CircObs(3.25, [10,32]), CircObs(3.75, [18,45]), ...
                CircObs(3.5, [37,41]), CircObs(3.5, [35,15]),...
                CircObs(3, [25,35]), CircObs(3.75, [13,10])];
            obs_mod = ObstacleMod(obstacles);
        end
    end
end

