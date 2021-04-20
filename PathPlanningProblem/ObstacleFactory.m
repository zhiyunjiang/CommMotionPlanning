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
        
        function obs_mod = circleGrid(r_base, r_plus_max, pos_rng, start, dest, scale)
            if nargin <=5
                scale = 1;
            end
            valid = 0;
            while ~valid
                circle_grid = [CircObs(r_plus_max*rand()+r_base, scale*[10,40] + pos_rng*(rand([1,2])-0.5) ), CircObs(r_plus_max*rand()+r_base, scale*[25, 10] + pos_rng*(rand([1,2])-0.5)),...
                    CircObs(r_plus_max*rand()+r_base, scale*[25,25] + pos_rng*(rand([1,2])-0.5)), CircObs(r_plus_max*rand()+r_base, scale*[25, 40] + pos_rng*(rand([1,2])-0.5)),...
                    CircObs(r_plus_max*rand()+r_base, scale*[40,10] + pos_rng*(rand([1,2])-0.5)), CircObs(r_plus_max*rand()+r_base, scale*[40, 25] + pos_rng*(rand([1,2])-0.5)),...
                    CircObs(r_plus_max*rand()+r_base, scale*[40,40] + pos_rng*(rand([1,2])-0.5))];

                obs_mod = ObstacleMod(circle_grid);
                if obs_mod.collisionFree(start, start) && obs_mod.collisionFree(dest, dest)
                    valid = 1;
                end
            end
        end
        
        function obs_mod = custom()
            obstacles = [CircObs(4.5, [15,40]), CircObs(3.25, [22,43]), ...
                CircObs(3.5, [36.75,44.25]), CircObs(4.5, [26,24]),...
                CircObs(3.5, [35.5,21]), CircObs(3.5, [28,11]),...
                CircObs(3.5, [38,7])];
            obs_mod = ObstacleMod(obstacles);
        end
        
        function obs_mod = fig2()
            obstacles = [CircObs(3.5, [12,45]), CircObs(4, [28,43]), ...
                CircObs(4, [29, 29]), CircObs(3.5, [38,38]),...
                CircObs(4, [43, 21]), CircObs(4, [21, 13]),...
                CircObs(3.5, [39,11])];
            obs_mod = ObstacleMod(obstacles);
        end
        
        function obs_mod = slide1()
            obstacles = [CircObs(2.8, [5, 19.5]), CircObs(2.4, [8,9]), ...
                CircObs(2.5, [13, 16.5]), CircObs(3, [21, 6]),...
                CircObs(1.8, [20, 20]), CircObs(1.4, [18, 12.5])];
            obs_mod = ObstacleMod(obstacles);
        end
    end
end

