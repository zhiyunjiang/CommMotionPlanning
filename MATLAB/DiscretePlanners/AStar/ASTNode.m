%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASTNode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A node in the tree build out by A*. Extends the TreeNode class, which is
% also used in RRT.


classdef ASTNode < TreeNode & handle
    
    % Properties (public):
    % cost2go - the value of the admissible heuristic for this node, i.e.
    %           the underestimated cost from this node to the destination
    % visited - whether or not we have already expanded this node/pulled
    %           this node off our priority queue
    properties
        cost2go = 0;
        visited = 0;
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % totalDist
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Finds the total estimated cost of the unique path that from root
        % to destination via this node.
        % Input:
        % this - reference to the AStarSolver object
        function total = totalDist(this)
            total = this.distToHere + this.cost2go;
        end
    end
end

