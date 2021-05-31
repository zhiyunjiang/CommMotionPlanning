classdef MCTSNode < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MCTSNode - A tree in the game tree built out by the MCTS algorithm
    
    properties
        %%%%%%%%%%%%%%%%%
        % Tree Variables
        %%%%%%%%%%%%%%%%%
        parent;
        children =[];
        
        %%%%%%%%%%%%%%%%%%
        % State Variables
        %%%%%%%%%%%%%%%%%%
        workspace_pos;
        visited_pos;
        
        %%%%%%%%%%%%%%%%%
        % MCTS Variables
        %%%%%%%%%%%%%%%%%
        num_visits = 0;
        r;
        expected_future_dist = 0;
        
        %%%%%%%%%%%%%%%%%%%%%%
        % Auxiliary Variables
        %%%%%%%%%%%%%%%%%%%%%%
        % moves we haven't strated exploring
        unexpanded_moves;
        % the J in Arjun's recursive characterization of the FPD integral
        J; 
        % probability of no connection up to here
        p_no_conn_to_here;
    end
    
    methods
        function this = MCTSNode(x, x_adj, parent)
            this.workspace_pos = x;
            this.unexpanded_moves = x_adj;
            
            ca = MCTSNode.setgetCa();
            if nargin == 3
               this.parent = parent;
               dist_to_parent = norm(x - parent.workspace_pos);
               parent.children = [parent.children(), this];
               
               % TODO - actually implement this function
               this.J = ca.JNext(parent.J, x, dist_to_parent);
               this.visited_pos = this.parent.visited_pos.deepCopy();
               this.visited_pos.add(parent.workspace_pos);
               this.filterUnexpandedMoves();
               
               this.p_no_conn_to_here = ca.IntegrateJ(this.J);
               this.r = dist_to_parent*this.parent.conditionalProbNoCon();
            else
                this.J = ca.J0(x);
                this.visited_pos = VisitedTable();
                this.p_no_conn_to_here = ca.IntegrateJ(this.J);
            end
            
        end
        
        function p = conditionalProbNoCon(this)
            denom = 1;
            if ~isempty(this.parent)
               denom =  this.parent.p_no_conn_to_here;
            end
           p = this.p_no_conn_to_here/denom;
        end
        
        function path = pathToHere(this)
            node = this;
            path=[];
            while ~isempty(node)
                path = [path; node.workspace_pos];
                node = node.parent;
            end
        end
        
        function d = expCost(this)
            d = this.r+ this.expected_future_dist;
        end
        
        function bst_child = getBestChild(this)
            bst_child = [];
            best_score = 0;
            for i=1:length(this.children)
               child = this.children(i);
               if best_score <  child.expCost()
                  best_score = child.expCost();
                  bst_child = child;
               end
            end
        end
        
        function tf = isFullyExpanded(this)
            tf = isempty(this.unexpanded_moves);
        end
        
        function addRewardFromOneSim(this, reward)
            if isempty(this.parent) % If we're the root
                scale = 1;
            else
                scale = this.parent.conditionalProbNoCon;
            end
                
            this.expected_future_dist = this.num_visits*this.expected_future_dist + scale*reward;
            this.num_visits = this.num_visits + 1;
            this.expected_future_dist = this.expected_future_dist/this.num_visits;
        end
        
        function [tf, won] =  sampleIsTerminal(this)
            if isempty(this.parent)
                %only the root has no parent
                tf = (rand() >= this.p_no_conn_to_here);
                
            else
                %need to calculate the probability of connectivity here
                %given no prior connectiivty, or rather, the probability of
                %no connectivity here given no prior connectivity)
                tf = (rand() <= (this.p_no_conn_to_here/this.parent.p_no_conn_to_here) );
                
            end
            won = tf;
            
            % handle the case that somehow we've backed ourselves into an
            % alley
            if isempty(this.children) && isempty(this.unexpanded_moves)
                tf = 1;
            end 
            
        end
        
        function new_child = expand(this, fpd_prob)
            if isempty(this.unexpanded_moves)
                new_child =[];
                warning('Tried to expand a node that was already fully expanded');
            else
                idx = randi(size(this.unexpanded_moves,1));
                new_x = this.unexpanded_moves(idx,:);
                %remove this from the unexpanded moves list
                this.unexpanded_moves(idx,:) = [];
                new_x_adjacent = fpd_prob.getNeighbors(new_x);
                new_child = MCTSNode(new_x, new_x_adjacent,this);
            end
        end
    end
    
    methods(Access = private)
        function filterUnexpandedMoves(this)
           to_remove = [];
           for i=1:length(this.unexpanded_moves)
              move = this.unexpanded_moves(i,:);
              [~, count] = this.visited_pos.pointInTable(move);
              % for now, we will disallow any cycles in the workspace
              if  count == 1 
                  %remove. this is a pointless cycle
                  to_remove = [to_remove,i];
              end
           end
           this.unexpanded_moves(to_remove,:) = [];
        end
        
        function cost = actionCost(action)
            if mod(action, 2) == 1%odd numbered action, moving diagonally
                cost = sqrt(2);
            else
                cost = 1;
            end
        end
    end
    
    methods(Access = public, Static)
        
        function out = setgetCa(ca)
            persistent Ca;
            if nargin
                Ca = ca;
            end
            out = Ca;
        end
    end
end

