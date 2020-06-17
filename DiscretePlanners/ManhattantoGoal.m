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

