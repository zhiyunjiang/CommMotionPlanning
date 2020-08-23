classdef PotentialGameSolver < handle
    %PotentialGameSolver - Solves Min-Exp-Cost-Simple_Path problem using
    %potential game formulation given in Arjun's Minimizing Expected Cost
    %to Success Paper
    
    properties (Constant = true)
        Actions = [[-1, 1]; [0, 1]; [1,1];...
                        [-1,0]; [1, 0];
                        [-1, -1];[0, -1];[1, -1]];
    end
    
    properties
        ca;
        pppi;
        mu;
        ups;
    end
    
    methods
        function this = PotentialGameSolver(pppi, ca)
            %PotentialGameSOlver - Construct an instance of
            %PotentialGameSolver
            %Inputs
            % pppi - PathPlanningProblem instance to be solved
            % ca - ChannelAnalyzer
            %Output
            % this - new PotentialGameSolver intance
            this.pppi = pppi;
            this.ca = ca;
        end
        
        function path = solve(this, max_iter)
            %solve -  Solves the path planning problem using log-linear
            %learning;
            grid_region = this.pppi.getGridRegion();
            dim_x = grid_region(1) + 1;
            dim_y = grid_region(3) + 1;
            %randomly with null policy
            this.initializeMuUps([dim_x, dim_y]);
            
            
            iter = 0;
            while iter <= max_iter
               %randomly select point
               
               x = randi(dim_x,1)-1; y = randi(dim_y,1)-1;
               point = [x,y];
               prev_action = this.mu(x+1, y+1);
               
               %get allowed actions
               actions = this.getActions(point);
               %find the costs associated with those action
               action_costs = this.getActionCosts(point, actions);
               %sample based on the action costs
               next_action = this.sampleNextAction(action_costs, 1/(iter+1) );
               
               %new action
               this.mu(x+1,y+1) = next_action;
               this.updateUps(point, prev_action, next_action);
               iter = iter + 1;
               
               if mod(iter, 500) == 0
                   fprintf('Completed %d of %d iterations\n',iter, max_iter);
               end
            end
            
            %now reconstruct the path from the given policy
            path = this.buildPath();
        end
    end
    
    methods (Access = private)
        
        function initializeMuUps(this, mu_dim)
            this.mu = zeros(mu_dim);
            this.ups = containers.Map('KeyType','double','ValueType','any');
            action_count = length(PotentialGameSolver.Actions);
            
            for i = 1:mu_dim(1)
                x = i - 1;
                for j = 1:mu_dim(2)
                   y = j-1;
                   pt = [x,y];
                   if min(size(this.pppi.collisionFree([pt;pt]))) == 2
                      action = randi(action_count,1);
                      this.mu(i, j) = action;
                      this.updateUps(pt, 0, action);
                   end
                end
            end
        end
        
        function path = buildPath(this)
            %getting source
            src = this.pppi.getSourceGrid();
            path = [src];
            next_action = this.mu(src(1)+1, src(2)+1);
            pt = src;
            path_length = 1;
            while next_action ~= 0
                path_length = path_length + 1;
                pt = pt + PotentialGameSolver.Actions(next_action,:);
                path(path_length) = pt;
                next_action = this.mu(pt(1)+1, pt(2)+1);
            end
        end
        
        function next_action = sampleNextAction(this, action_costs, t)
            exp_costs = exp((-1/t)*action_costs);
            if ~any(exp_costs)
                next_action = 0;
            else
                probs = exp_costs/sum(exp_costs);

                cdf = cumsum(probs);
                r = rand(1);
                next_action = find(cdf>r, 1, 'first');
            end
        end
        
        function allowed_actions = getActions(this, point)
            %getActions - finds the admissible actions from a specific
            %point in the configuration space.
            %Input
            % this - current PotentialGameSolver
            % point - the point in question
            %Output
            % allowed_actions - boolean vector indicating wheter or not
            % each action is admissible from the point in question.
            %Possible actions are up, upright, left, down left, down, down
            %right, right and up right - eigth actions total, numbered as
            %shown below:
            %   1   2   3
            %   4   c   5
            %   6   7   8
            allowed_actions = zeros([length(PotentialGameSolver.Actions),1]);
            %now check each action
            pppi = this.pppi;
           for i=1:8
              neighbor = point + PotentialGameSolver.Actions(i,:); 
              if pppi.gridPtInRegion(neighbor)
                  %now check for obstacles
                  sub_path = pppi.collisionFree([point; neighbor]);
                  if sum(size(sub_path)) == 4  && ~this.inUps(point, neighbor)
                    allowed_actions(i) = 1;
                  end
              end
           end   
        end
        
        function costs = getActionCosts(this, point, allowed_actions)
            %getActionCosts - Find the costs associated with each allowed
            %action, given the current policy mu
            %Input
            % point - current point, that is, the player whose action we're
            % interested in
            % allowed_actions = boolean vector of allowed actions
            %Output
            % costs - array of the costs associated with each action. If an
            % action is not allowed, the cost will be set of Inf.
           
            costs = Inf*ones([length(PotentialGameSolver.Actions), 1]);
            
            J0 = this.ca.J0(point);
            p_no_conn = this.ca.IntegrateJ(J0);
            
            for i = length(PotentialGameSolver.Actions)
              if allowed_actions(i)
                 next_point = point + PotentialGameSolver.Actions(i,:);
                 dist = norm(point - next_point);
                 cost = p_no_conn*dist;
                 J1 = this.ca.ItterativeJNextFromPoint(J0, point, dist);
                 costs(i) = cost + this.iterativeCost(next_point, J1);
              end
           end
               
        end
        
        function cost = iterativeCost(this, point, J)
            action = this.mu(point(1) + 1, point(2) + 1);
            if action == 0
                if this.pppi.pointInGoalRegion(point)
                    cost = 0;
                else
                    cost = Inf;
                end
            else
                next = point + PotentialGameSolver.Actions(action,:);
                [osc, J_next] = this.oneStepCost(J, point, next);
                cost = osc + this.iterativeCost(next, J_next);
            end
            
        end
        
        function [osc, J_next] = oneStepCost(this, J, point, next_point)
            p_no_conn = this.ca.IntegrateJ(J);
            dist = norm(point - next_point);
            osc = p_no_conn*dist;
            J_next = this.ca.ItterativeJNextFromPoint(J, point, dist);
        end
        
        function tf = inUps(this, start, needle)
           %inUps - Checks if the needle is in the haystock of points upstream
           %from the given starting node 
           %Inputs
           % ups - 2-d map x->y->[points whose current action is to move to (x,y)]
           % start - point to start from
           %needle - point we're looking for
           queue = [start];
           tf = 0;
           while ~isempty(queue)
               next = queue(1,:);
               if min(size(queue) == 1)
                   queue = [];
               else
                   queue = queue(2:end,:);
               end
               
               if this.ups.isKey(next(1))
                   submap = this.ups(next(1));
                   if submap.isKey(next(2))
                      direct_ups = submap(next(2)); 
                      diffs = direct_ups - needle;
                      if all(diffs) %check that all diffs are non-zero
                          queue = [queue;direct_ups];
                      else
                          %one of the diffs was non-zero
                         tf = 1;
                         break;
                      end
                   end
               end
           end
        end
        
        function updateUps(this, point, prev_action, current_action)
           % if it was 0, we weren't point to anything from this point, don't need to stop point to anything
            if prev_action ~= 0 
               prev_next_point = point + PotentialGameSolver.Actions(prev_action,:);
               this.removeFromUps(prev_next_point, point);
            end
           
            if current_action ~= 0
               curr_next_point = point + PotentialGameSolver.Actions(current_action,:);
               this.addToUps(curr_next_point, point); 
            end  
        end
        
        function addToUps(this, index_point, new_point)
            %first check if the x map exits
            x = index_point(1); y = index_point(2);
            if ~this.ups.isKey(x)
                this.ups(x) = containers.Map('KeyType','double','ValueType','any');
            end
            submap = this.ups(x);
            if ~ submap.isKey(y)
                submap(y) = [];
            end
            submap(y) = [submap(y);new_point];
        end
        
        function removeFromUps(this, index_point, new_point)
            submap = this.ups(index_point(1));
            point_list = submap(index_point(2));
            rm_index = find( sum(abs(point_list - new_point), 2)==0);
            point_list(rm_index) = [];
            submap(index_point(2)) = point_list;
        end
    end
end

