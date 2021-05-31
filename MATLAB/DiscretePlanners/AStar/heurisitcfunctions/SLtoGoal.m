%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLtoGoal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Straight line euclidean distance to nearest goal point.
%
% Inputs:
% grid_pt - [x,y] position of the node from which we are calculating the
%           heuristic
% pppi - PathPlaningProblem instance
%
% Outputs:
% estimated2go - (under)estimated distance to go from grid_pt to goal.

function estimated2go = SLtoGoal(grid_pt,pppi)
    goals = pppi.getGoalRegion.goalGridPoints;
    
    estimated2go = Inf;
    %If we store goals in balanced box decomposition tree, we could look
    %this up more quickly
    [g_count, ~] = size(goals);
    for i=1:g_count
        dist = norm(grid_pt - goals(i,:));
        if dist < estimated2go
           estimated2go = dist; 
        end
    end
end

