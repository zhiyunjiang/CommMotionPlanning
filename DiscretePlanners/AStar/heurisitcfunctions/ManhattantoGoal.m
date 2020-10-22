%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ManhattantoGoal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manhattan distance to nearest goal point. Only admissible if we can only
% move up, down, left, right
%
% Inputs:
% grid_pt - [x,y] position of the node from which we are calculating the
%           heuristic
% pppi - PathPlaningProblem instance
%
% Outputs:
% estimated2go - (under)estimated distance to go from grid_pt to goal.

function estimated2go = ManhattantoGoal(grid_pt,pppi)
     goals = pppi.getDestGrid();
    
    estimated2go = Inf;
    %If we store goals in balanced box decomposition tree, we could look
    %tihs up more quickly
    [g_count, ~] = size(goals);
    for i=1:g_count
        dist = norm(grid_pt - goals(i,:), 1);
        if dist < estimated2go
           estimated2go = dist; 
        end
    end
end

