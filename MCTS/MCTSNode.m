classdef MCTSNode < handle
    %MCTSNode Node in the tree built out by MCTS
    
    properties
        wrkspcpos;
        visited;
        parent;
        children = [];
        isRoot = 0;
        distToHere = Inf;
    end
    
    methods
        function this = MCTSNode(pos, parent)
            %MCTS - Contructor
            this.wrkspcpos = pos;
            if nargin == 2
                this.parent = parent;
                this.visited = parent.visited;
                this.addVisited(parent.wrkspcpos);
            else
                this.visited = containers.Map('KeyType','double','ValueType','any');
            end
        end
        
        function makeRoot(this)
           this.isRoot = 1; 
        end
        
        function pos = getPos(this)
            pos = this.wrkspcpos;
        end
        
        function addVisited(this, pos)
           if ~this.visited.isKey(pos(1))
               this.visited(pos(1)) = containers.Map('KeyType','double','ValueType','any');
           end
           
           xmap = this.visited(pos(1));
           
           xmap(pos(2)) = 1;
        end
        
        function tf = hasVisited(this, pos)
           tf = 0;
           if this.visited.isKey(pos(1)) && this.visited(pos(1)).isKey(pos(2))
                tf = 1;
           end
        end
        
        function path = pathToHere(this)
           node = this;
           path = node.getPos();
           while ~node.isRoot
              node = node.parent;
              path = [node.getPos(); path];
           end
        end

    end
end

