classdef ASTNode < TreeNode & handle
    %ASTNode - A* treenode class
    
    properties
        cost2go = 0;
        visited = 0;
    end
    
    methods
        %just use TreeNode superclass constructor
        
        function total = totalDist(o)
            total = o.distToHere + o.cost2go;
        end
    end
end

