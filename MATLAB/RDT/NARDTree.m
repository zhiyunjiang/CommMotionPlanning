%non-additive RD tree
classdef NARDTree < handle
    
    properties
        %doRewire - if true, will rewire, that is, will run RRT* rather
        %than RRT
        doRewire;
        
        %steerRad - maximum distance the tree can grow at a time, that is
        %maximum distance between p_nearest and p_new
        steerRad;
        
        %BSF = TreeNode. Path from BSF to root traces the best path found
        %so far
        BSF;
        
        %treeNodes = [root, n1, n2, ....], the list of nodes in the tree,
        %in order in which they are added.
        treeNodes;
        
        %see "Sampling-based Algorithms for Optimal Motion Planning",
        %s. Karaman, E. Frazzoli, 2011 
        gamma;
        
        %
        solCosts;
        solDists;
        
        %neighborBoxes - data structure presented in "Minimising 
        %computational complexity of the RRT algorithm a practical
        %approach" M Sventrup et. al, 2011
        %neighborBoxes;
        
        %
        isContinuous;
        
        theta;
    end
    
    methods
        function this = RDTree(do_rewire, steer_rad)
           this.doRewire = do_rewire;
           this.steerRad = steer_rad;
            
        end
        
        %setup asspects of the tree before beginning the algorithm
        function initialize(this, gamma, root_pos, is_continuous, theta)
            if nargin == 4
                theta = 1;
            end
            
            this.theta = theta;
            this.gamma = gamma;
            
            root = TreeNode(root_pos);
            root.isRoot = 1;
            root.distToHere = 0;
            this.treeNodes = root;
            
            %create a placeholder BSF with infinite cost
            this.BSF = TreeNode([-Inf,-Inf]);
            this.BSF.distToHere = Inf;
            
            if nargin == 4
                this.isContinuous = is_continuous;
            else
                this.isContinuous = 0;
            end
            
        end
        
        function executeIteration(this, x_rand, pppi)
            
            [success, new_node] = this.tryAddNode(x_rand, pppi);
            
            %check if the new_node has given us a new BSF
            if success && pppi.nodeInDestGrid(new_node) && this.BSF.distToHere > new_node.distToHere 
                this.BSF = new_node;
            end
        end
        
        function [success, n_new] = tryAddNode(this, x_rand, pppi)
            success = 0;
            %initialize dummy new node
            n_new = TreeNode([Inf, Inf]);
            [nearest, min_dist] = this.nearest(x_rand);
            %check to see if we already have this node added to the tree
            %should be able to just check if zero again, but let's validate
            if  min_dist ~= 0 
                path = this.steer(x_rand, nearest, pppi);
                %check if the path is collision free, trimming if necessary
                viable_path = pppi.collisionFree(path);
                
                if min(size(viable_path)) > 1 % more than one point in path
                    success = 1;
                    n_new = this.addNode(nearest, viable_path, this.doRewire, pppi);
                end
            else
                %handle the case where we resample an already sampled
                %point. Need to rewire. Look into literature
                this.rewire(nearest, pppi)
            end
        end
         
        function recordBSFCost(this, time)
            if this.BSF.distToHere ~= Inf
                this.solCosts = [this.solCosts; time, this.BSF.distToHere];
                path = this.BSF.pathToRoot(0);
                this.solDists = [this.solDists; time, GridDist(path)];
            end
        end
    end
    
    methods (Access = private)

         %addNode
         function new_node = addNode(this, nearest, path, do_rewire, pppi)
            new_node = TreeNode(path(end,:));
            %mode = 2 - update node data (if distance metric requires it)
            cost = pppi.pathCost(new_node, nearest, path, 2);
            new_node.setParent(nearest, cost, path); 
            if do_rewire
                this.rewire(new_node, pppi);
            end
            this.treeNodes = [this.treeNodes, new_node];
         end
        
        %rewire
        % TODO - rewire also includes finding shortest path (not just
        % closest TO the new node)
        % TODO - implement with balanced box decomposition
        %can be implemented in O(log n) time using Balanced Box Decomposition
        %currently is O(n) (brute force)
        function rewire(this, new_node, pppi)
            %See Karaman & Frazzoli, 2011
            cardV = length(this.treeNodes);
            radius = min(this.steerRad, this.gamma*sqrt(log(cardV)/cardV));
            neighbors = this.near(new_node.wrkspcpos, radius);
            
            
            cost_btwn = -1*ones(size(neighbors));
            x_current = new_node.wrkspcpos;
            if new_node.isRoot
                min_cost = 0;
                min_cost_btwn = 0;
            else
                best_neighbor = new_node.parent;
                min_cost = new_node.distToHere;
                min_cost_btwn = min_cost - best_neighbor.distToHere;
                best_path = new_node.pathFromParent;
            end
            
            %first, find the least cost path to current
            for i = 1:length(neighbors)
                neighbor = neighbors(i);
                if ~new_node.isRoot && neighbor == new_node.parent
                    continue;
                end
                %want path from from  neighbor to current
                path = pppi.getSLPath(neighbor.wrkspcpos, x_current);
                
                %check to make sure there's nothing obstructing
                viable_path = pppi.collisionFree(path);
                if norm(size(viable_path) - size(path)) == 0
                    %full path works!
                    cost_btwn(i) = this.theta*pppi.pathCost(new_node, neighbor, path, 1);
                    cost_to_new_via_neighbor = cost_btwn(i) + neighbor.distToHere;
                    if cost_to_new_via_neighbor < min_cost              
                           best_path = path;
                           best_neighbor = neighbor;
                           min_cost = cost_to_new_via_neighbor;
                           min_cost_btwn = cost_btwn(i);
                    end
                end
            end
            
            if ~new_node.isRoot
                new_node.setParent(best_neighbor, min_cost_btwn, best_path); 
            end
           
            %now check if there are any nodes that we can reach with less
            %cost from our new_node
            for j = 1:length(neighbors)
               neighbor = neighbors(j);
               cost = cost_btwn(j);
               if cost >= 0

                    if (neighbor.distToHere > (new_node.distToHere + cost)) 
                       %we've already seen that there's nothing
                       %obstructing from the first time we looped through

                       %path from new to neighbor
                       path = pppi.getSLPath(x_current, neighbor.wrkspcpos);
                       %now add the new parent
                       neighbor.setParent(new_node, cost, path);
                    end
               end
            end
        end
          
        function path = steer(this, x_rand, v_nearest, pppi)

            start = v_nearest.wrkspcpos;
            diff = x_rand - start;
            %don't extend beyond the newly sampled point
            if norm(diff) > this.steerRad
            
                %find the direction
                phi = atan2(diff(2), diff(1));

                new_diff = [this.steerRad*cos(phi), this.steerRad*sin(phi)];
                new = start + new_diff;
                %snap to grid
                new = pppi.toRawCoordinate(pppi.toGridCoordinate(new));
            else
               new =  x_rand;
            end
            
           path = pppi.getSLPath(start, new);
            
        end
        
        function neighbors = near(this, pt, radius)
            neighbors = [];
             for i = 1:length(this.treeNodes)
               node = this.treeNodes(i);
               dist = norm(node.getPos() - pt);
               if dist < radius && dist ~= 0%ignore if it's the same point
                  neighbors = [neighbors, node];
               end
            end
        end

        % TODO - implement with balanced box decomposition (space O(n), 
        % time O(log(n))) or boxed (space = grid_x_dim * grid_y_dim, time O(1))
        %currently is O(n) (brute force)
        
        % Aslo, not sure we need to do this AND the rewire routine. Prollay
        % just need the rewire. Will investigate in a hot min tho.
        function [nearest, min_dist] = nearest(this, x_rand)
            root = this.treeNodes(1);
            min_dist = norm(root.getPos() - x_rand);
            nearest = root;
            for i = 2:length(this.treeNodes)
               dist = norm(this.treeNodes(i).getPos() - x_rand);
               if dist < min_dist
                  min_dist = dist;
                  nearest = this.treeNodes(i);
               end
            end
        end
         
    end
    
end

