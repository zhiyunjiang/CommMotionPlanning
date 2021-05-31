%parent_node - parent tree node object
%path - [parent coords, ....., child coords]

function delta_dist = MinExpFPDInd(parent_node, path, prob_field, gamma_th)
    if length(path) == 1
        delta_dist = 0;
    else
        
        cost = parent_node.distToHere;
        for i = 2:length(path)
            step = norm(path(i-1,:) - path(i,:));
            p_conn = prob_field.posteriorPConn(path(i,:), gamma_th);
            cost = (1-p_conn)*step - p_conn*cost;
        end
        %handle numerical issues that are pushing cost to very small
        %negative numbers. Cost should never be able to fall below 0,
        %theoretically
        cost = max(0, cost);
        
        delta_dist = cost - parent_node.distToHere;
    end
end

