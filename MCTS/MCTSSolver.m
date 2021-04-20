classdef MCTSSolver < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MCTSSolver - class that solves a path planning problem using Monte
    %               Carlo Tree Search
    
    properties
        root; %the root MCTS node
        stop_criteria; %conditions on which we should stop
        reward_constant; %reward will be of the form (constant - length)/constant
        channel_analyzer;
        dp_rrt; % RRT Solver used to produce a default policy
        reward_offset;
    end
    
    methods (Access = public)
        function this = MCTSSolver(stop_criteria, channel_analyzer, rrt_solver)
            this.stop_criteria = stop_criteria;
            this.channel_analyzer = channel_analyzer;
            this.dp_rrt = rrt_solver;
            this.reward_offset = 50*50*5;%TODO - dynamically readjust based on size of region
            MCTSNode.setgetCa(this.channel_analyzer);% TODO - decouple the MCTsSolver and MCTSNode
        end
        
        function solve(this, pppi)
           iteration_num = 0;
           elapsed_time = 0;
           tic;

           %initialize the root
           this.initializeMCTree(pppi);
           
           %run iterations until stopping criteria met
           while ~this.stop_criteria.stop(elapsed_time, iteration_num, [])
              iteration_num = iteration_num + 1;
              %run the node selection policy
              selected_node = this.treePolicy(pppi);
              %play out the default policy
              reward = this.defaultPolicy(selected_node, pppi);
              % back-propagate the reward
              MCTSSolver.backProp(reward, selected_node);
              elapsed_time = toc;
              if mod(iteration_num, 20)== 0
                  fprintf("%d simulations complete\n",iteration_num);
              end
           end
        end
        
        function path = getOptiamlPolicy(this)
           
            node = this.root;
            path = [node.workspace_pos];
            while ~isempty(node.children)
                node = node.getBestChild();
                if isempty(node)
                    break;
                end
                try
                    path = [path; node.workspace_pos];
                catch
                    warning("this is a werid error that should never happen");
                end
            end
        end
    end
    
    methods (Access = private)
        function node = treePolicy(this, fpd_prob)
            %there are a number of ways to select the next node. We will
            %use the standard UCT (Upper confidence bound for tree search)
            
            parent = this.root;
            while parent.isFullyExpanded()%handle case where we've reached a terminal node
                best_index = -1;
                best_score = this.reward_offset;
                if isempty(parent.children)
                    break;
                end
                for i=1:length(parent.children)
                   child = parent.children(i);
                   %use UCT (Upper Confidence Bound (UCB) for Trees)
                   child_score = MCTSSolver.calculateUCB(parent, child);
                   if child_score < best_score
                      best_score = child_score;
                      best_index = i;
                   end
                end
                
                if best_index == -1
                   error("No child was better than default offset");
                end
                parent = parent.children(best_index);
            end
            
             if parent.isFullyExpanded()% there will be no children to expand, possibly because there are no more valid moves
                 node = parent;
             else
                %node wasn't fully expanded, so let's continue to expand it
                node = parent.expand(fpd_prob);
             end
            
        end
        
        function dist = defaultPolicy(this, node, pppi)
            % TODO - just use RRT for the default policy
            % don't do full RRT, just use the sampler to perform random
            % samples. This way we don't sample more than we need to, i.e.
            % Use RRT to find a path to the remote station
            dp_pppi = pppi.copyWithNewSource(node.workspace_pos);
            %dp_pppi.obstacleMod.obstacles(end+1) = PointObstacle(node.visited_pos);

            this.dp_rrt.solve(dp_pppi);

            bsf = this.dp_rrt.getBSF();
            rrt_path = bsf.pathToRoot(0);
            J = node.J;%need to find distance conditioned on no connection to here
            %now find expected distance along this path
            [approx_PMF, distances] = this.channel_analyzer.ApproxFPDPMF(rrt_path, 1, J);
            dist = approx_PMF'*distances;
        end
        
        function initializeMCTree(this,fpd_prob) 
           x0 = fpd_prob.getSource();
           x0_adjacent = fpd_prob.getNeighbors(x0);
           this.root = MCTSNode(x0, x0_adjacent);
        end
       
    end
    
    methods (Access = private, Static)
        
        
        function backProp(reward, node)
            node.addRewardFromOneSim(reward);
            while ~isempty(node.parent)
                reward = node.parent.conditionalProbNoCon()*reward + node.r;
                node = node.parent;
                node.addRewardFromOneSim(reward);
            end
        end
        
        function score = calculateUCB(parent, child)
            score = child.expCost ...
                       - 0.1*sqrt(log(parent.num_visits)/child.num_visits);
                   %subtract since we're acutally looking to minimize
        end
        
    end
end

