classdef MCTSSolver
    %MCTSSolver - Solves a path planning problem instance using Monte Carlo
    %Tree Search
    
    properties (Constant = true)
        Actions = [[-1, 1]; [0, 1]; [1,1];...
                        [-1,0]; [1, 0];
                        [-1, -1];[0, -1];[1, -1]];
    end
    
    methods (Static = true)
        
        function [root,node_2_expand] = solve(pppi, ca, max_depth, n_sims)
            %solve - solves the path planning problem instance using Monte
            %Carlo Tree Search
            %INPUT
            % pppi - path planning problem instance to be solved
            % max_depth - maximum number of nodes to add to tree
            % n_sims - number of simulations to perform
            %
            %OUTPUT
            % root - MCTS root node
            
            action_count = length(MCTSSolver.Actions);
            grid_region = pppi.getGridRegion();
            lower_bound = grid_region(1)*grid_region(3);
            %initialize tree (root node)
            root = MCTSNode(pppi.getSourceGrid());
            root.makeRoot();
            node_2_expand = root;
            iteration_count = 0;
            rrt_solver = MCTSSolver.getRRTSimulator();
            while iteration_count < max_depth && ~(pppi.nodeInDestGrid(node_2_expand))
                tic;
                iteration_count = iteration_count + 1;
                %already have node_2_expand selected
                admissible_actions = MCTSSolver.getAdmissibleActions(node_2_expand, pppi);
                %TODO - consider precomputing all available actions (or at
                %least collisions)

                %TODO - consider simulating each action in parallel
                bsf = 0;%best action so far
                bcsf = Inf;%best cost so far
                costs = zeros([action_count, 1]);
                counts = zeros([action_count, 1]);
                for i = 1:n_sims
                   action = MCTSSolver.getUCB1Next(admissible_actions, costs, counts, i-1, lower_bound);
                   counts(action) = counts(action) + 1;
                    
                   %run a single simulation, starting from 
                   next_pos = MCTSSolver.nextPoint(node_2_expand, action);
                   pppi_sim = pppi.copyWithNewSource(next_pos);
                   
                   costs(action) = costs(action) +  MCTSSolver.simulateFrom(pppi_sim, ca, rrt_solver);
                   
                   exp_cost = costs(action)/counts(action);

                   if exp_cost <= bcsf
                       bcsf = exp_cost;
                       bsf = action;
                   end
                   %since we're committing to an action at each step,
                   %we don't need to back prop
                end
                
                %now get the next point
                %TODO - handle case where bsf = 0;
                node_2_expand = MCTSSolver.getNextNode(node_2_expand, bsf);
                fprintf('Completed iteration %d. Duration %d\n', iteration_count, toc);
            end
            
        end
        
        function next_action = getUCB1Next(allowed, costs, plays, total_plays, lower_bound)
            i_allowed = find(allowed==1, 8);
            allowed_costs = costs(i_allowed);
            allowed_plays = plays(i_allowed);
            
            if ~all(allowed_plays)
                ind = find(allowed_plays==0,1);
            else
                %use UCB1 to decide which to use
                explore = sqrt(2*log(total_plays)/allowed_plays);
                %normalize to keep in [0, 1]
                exploit = 1 - (allowed_costs./allowed_plays)/lower_bound;
                [~,ind] = max(explore'+exploit);
                if length(ind) > 1
                   ind = randi(length(ind)); 
                end
            end
            
            next_action = i_allowed(ind);
        end
        
        function next_node = getNextNode(node_2_expand, action)
            next_pos = MCTSSolver.nextPoint(node_2_expand, action);
            next_node = MCTSNode(next_pos, node_2_expand);
        end
        
        function next_pos = nextPoint(node, action)
            next_pos = node.getPos() + MCTSSolver.Actions(action,:);
        end
        
        function admissible_actions = getAdmissibleActions(node, pppi)
            %getAdmissibleActions - finds the admissible actions from a specific
            %point in the configuration space.
            %Input
            % node - current node
            % pppi - the path planning problem instance
            %Output
            % admissible_actions - boolean vector indicating wheter or not
            % each action is admissible from the point in question.
            %Possible actions are up, upright, left, down left, down, down
            %right, right and up right - eigth actions total, numbered as
            %shown below:
            %   1   2   3
            %   4   c   5
            %   6   7   8
           admissible_actions = zeros([length(MCTSSolver.Actions),1]);
           %now check each action
           point = node.getPos();
           for i=1:8
              neighbor = point + MCTSSolver.Actions(i,:); 
              if pppi.gridPtInRegion(neighbor)
                  %now check for obstacles
                  sub_path = pppi.collisionFree([point; neighbor]);
                  if sum(size(sub_path)) == 4  && ~node.hasVisited(neighbor)
                    admissible_actions(i) = 1;
                  end
              end
           end   
        end
        
        function cost = simulateFrom(pppi_sim, ca, rrt_solver)
            rrt_solver.solve(pppi_sim);
            bsf = rrt_solver.getBSF();
            path = bsf.pathToRoot(0);
            cost = ca.ExpectedFPD(path,4);
        end
        
        function rrt_solver = getRRTSimulator()
            %solve using RRT
            %TODO - just use the same solver, don't waste time
            %re-intializing
            %stop on first solution
            solution_not_required = 0;
            max_iterations = Inf;
            %max run time in seconds
            max_run_time = 0;
            stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);
            %1 - uniform random (RRT), 2 - deterministic sequence (RDT) 
            %3 - informed set sampling, 4 - continuous uniform random
            type = 1;
            dest_freq = 20;
            sequence = [];%sampler.sequence;%use previous run's sequence
            sampler = Sampler(type, dest_freq, sequence);

            do_rewire = 0;%just need RRT, not RRT*
            steer_rad = 5;

            rrt_solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);
        end
    end
end

