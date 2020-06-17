classdef BridgeNode < handle
    
    properties (Access = public)
        optimistic_heuristic;
        
        pos;
        
        leafFromSource;
        
        leafToDest;
    end
    
    methods
        function this = BridgeNode(pos, optimistic_heuristic)

            this.pos = pos;
            
            this.optimistic_heuristic = optimistic_heuristic;
        end
        
        function cost =  costToDest(this)
            if isempty(this.leafToDest)
                cost = Inf;
            else 
                cost = this.leafToDest.distToHere;
            end
        end
        
        function cost =  costFromSource(this)
            if isempty(this.leafFromSource)
                cost = Inf;
            else
                cost = this.leafFromSource.distToHere;
            end
        end
        
        function bsf = getBestPathCost(this)
            if isempty(this.leafToDest) || isempty(this.leafFromSource)
                bsf = Inf;
            else
                bsf = this.leafToDest.distToHere + this.leafFromSource.distToHere;
            end
        end
        
        function dist = distToHere(this)
           dist = this.getBestPathCost(); 
        end
        
        function path =  pathToRoot(this, do_print)
           path_to = this.leafFromSource.pathToRoot(do_print); 
           hold on
           path_from = this.leafToDest.pathToRoot(do_print);
           
           path = [path_to;flip(path_from)];
        end
        
    end
end

