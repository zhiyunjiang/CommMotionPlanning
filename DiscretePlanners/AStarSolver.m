classdef AStarSolver < handle
    %AStarSolver Summary of this class goes here
    
    properties
        heuristicFcn;
        isDijkstras;
        pppi;
        BST;
        root;
    end
    
    methods (Access = public)
        function this = AStarSolver(heuristic_fcn)
            if nargin == 1
               this.heuristicFcn = heuristic_fcn;
               this.isDijkstras = 0;
            else
                this.heuristicFcn = @(n, dest) AStarSolver.DfltHeuristic(n, dest);
                this.isDijkstras = 1;
            end
        end
        
        function solve(this, pppi)
            this.pppi = pppi;
            existing_nodes = containers.Map('KeyType','double','ValueType','any');
            
            %implement priority queue using RedBlack tree so that
            %insertion, searching, takes log(n) time where n is queue size.
            
            p_queue = RedBlackTree();
            start_pt = pppi.getSourceGrid();
            this.initializeTree(start_pt, existing_nodes);
            %get starting point, create a node, pop it onto the queue
            p_queue.Insert(this.root.totalDist, this.root);
            
            found = 0;
            
            while ~p_queue.IsEmpty() && ~found
               %pop off the shortest from the tree/queue
               current = AStarSolver.Dequeue(p_queue);
               if current.visited
               %Will be using consistent heuristic, therefore once we've
               %visited once, we're set
                   continue;
               end
               
               current.visited = 1;
               
               if pppi.pointInGoalRegion(current.getPos())
                  found = 1;
                  break;
               end
               
               %siblings, per the grid
               neighbor_pts = pppi.getGridNeighbors(current.getPos());
               sz = size(neighbor_pts);
               neighbor_count = sz(1);
               for i = 1:neighbor_count
                   n_grid_pt = neighbor_pts(i,:);
                   [neighbor_node, is_new] = this.getNodeAtGridPt(existing_nodes, n_grid_pt);
                   
                   if is_new
                       %set the new parent, which will automatically update
                       %distance to there
                       neighbor_node.setParent(current, 1);
                        
                        %now push onto the priority queue
                        p_queue.Insert(neighbor_node.totalDist, neighbor_node);
                        
                   elseif ~neighbor_node.visited && neighbor_node.distToHere > current.distToHere + 1
                       
                       neighbor_node.parent.removeChild(neighbor_node);
                       
                       neighbor_node.setParent(current, 1);
                       
                       %now push onto the priority queue
                        p_queue.Insert(neighbor_node.totalDist, neighbor_node);

                   end
               end

            end
            
            this.BST = current;
        end
        
    end

    methods (Access = private)
        
        function root = initializeTree(this, start_pt, existing_nodes)
           [root, ~] = this.getNodeAtGridPt(existing_nodes, start_pt); 
           root.isRoot = 1;
           root.distToHere = 0;
           this.root = root;
        end
        
        function [node, is_new] = getNodeAtGridPt(this, existing_nodes, grid_pt)
            x = grid_pt(1);
            y = grid_pt(2);
            is_new = 1;
            if existing_nodes.isKey(x)
                ymap = existing_nodes(x);
                
                if ymap.isKey(y)
                   node = ymap(y);
                   is_new = 0;
                else
                    node = ASTNode(grid_pt); 
                    node.cost2go = this.heuristicFcn(grid_pt, this.pppi);
                    ymap(y) = node;
                end
                
            else
                
                node = ASTNode(grid_pt); 
                node.cost2go = this.heuristicFcn(grid_pt, this.pppi);
                ymap = containers.Map(y,node);
                existing_nodes(x) = ymap;
                
            end 
        end
        

    end
    
    methods (Access = private, Static = true)
        
        %define the default heuristic function
        function cost2go = DfltHeuristic(p1, dest)
           cost2go = 0;
        end
        
        function next_AS_node = Dequeue(p_queue)
            queue_element = p_queue.Minimum();
            next_AS_node = queue_element.value;
            p_queue.Delete(queue_element);
        end

    end
end

