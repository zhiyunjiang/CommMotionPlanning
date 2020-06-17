classdef MidSetSolver < handle
    
    properties (Constant)
       SOURCE = 0;
       DEST = 1;
    end
    
    
    properties (Access = private)
        %sampler - sampling module to use
        sampler;
        
        %stopCriteria - function handle. Function determines whether or not
        %we should stop. Arguments are elapsed time, iteration number, and cost
        %of best solution so far
        stopCriteria;
        
        %Best Bridge so Far
        BBSF
        
        sourceRDT;
        
        destRDT;
        
        currentTreeID;
        
        bridgeNodes;
        
        histBBSFCost;
        
    end
    
    methods (Access = public)
        function this = MidSetSolver(sampler, stop_criteria, do_rewire, steer_rad)
           this.sampler = sampler;
           sampler.setHalfCount(1);
           
           this.stopCriteria = stop_criteria;
           
           this.sourceRDT = RDTree(do_rewire, steer_rad);
           this.destRDT = RDTree(do_rewire, steer_rad);
        end
        
        function solve(this, pppi, theta)
            if nargin == 2
                theta = 1;
            end
            
            tic;
            this.initializeTrees(pppi, theta);
            this.sampler.setParams(pppi);
            this.initializeBridgeNodes(pppi)
            
            iteration_count = 0;
            elapsed_time = toc;
            this.currentTreeID = MidSetSolver.SOURCE;
            current_tree = this.sourceRDT;
            
            series_delta = 0.2;
            time_for_next_recording = series_delta;
            
            while ~this.stopCriteria.stop(elapsed_time , iteration_count, this.BBSF)
                current_BSF = current_tree.BSF;
                x_rand = this.sampler.sample(iteration_count, pppi, current_BSF.distToHere());
                
                [success, new_node] = current_tree.tryAddNode(x_rand, pppi);
                
                %check if the new_node has given us a new BSF
                if success && pppi.nodeInDestGrid(new_node)
                    if current_BSF.distToHere() > new_node.distToHere 
                        current_tree.BSF = new_node;
                    end
                    
                    this.updateBridgeNodes(new_node);
                end
                
                
                iteration_count = iteration_count + 1;
                
                elapsed_time = toc;
                
                 if elapsed_time > time_for_next_recording
                    time_for_next_recording = time_for_next_recording + series_delta;
                   this.recordBBSFHist(elapsed_time); 
                 end
                 
                current_tree = this.getOtherTree();
                
            end
            
        end
        
        function series = getHistBBSFCost(this)
            series = this.histBBSFCost;
        end
        
        function recordBBSFHist(this, time)
            if this.BBSF.distToHere() ~= Inf
                this.histBBSFCost = [this.histBBSFCost; time, this.BBSF.distToHere()]; 
            end
        end
        
        function tail = getBBSF(this)
            tail = this.BBSF;
        end
    end
    
    methods (Access = private)
        
        function other = getOtherTree(this)
           if this.currentTreeID == MidSetSolver.SOURCE
               other = this.destRDT;
               this.currentTreeID = MidSetSolver.DEST;
           else
               other = this.sourceRDT;
               this.currentTreeID = MidSetSolver.SOURCE;
           end
        end
        
        function initializeTrees(this, pppi, theta)
           %see Karaman & Frazzoli, 2011, page 29. Here, d=2 
           gridRegion = pppi.getGridRegion();
           length = (gridRegion(1) - gridRegion(2));
           width = (gridRegion(3) - gridRegion(4));
           grid_area = length*width;%total area will be >= area of X_free
           gamma =  sqrt(3*(grid_area/2*pi));
           
           is_continuous = ( pppi.getStepSize() == 0 );
            
            src_root_pos = pppi.getSourceGrid();
            this.sourceRDT.initialize(gamma, src_root_pos, is_continuous);
            
            dest_root_pos = pppi.getSourceGrid2();
            this.destRDT.initialize(gamma, dest_root_pos, is_continuous, theta);
        end
        
        
        function improved_connection =  updateBridgeNodes(this, new_node)
            improved_connection = 0;
            pos =  new_node.getPos();
            sub_map = this.bridgeNodes(pos(1));
            bridge_node = sub_map(pos(2));
            if this.currentTreeID == MidSetSolver.SOURCE
               bsf = bridge_node.costFromSource();
               if new_node.distToHere() < bsf
                  bridge_node.leafFromSource = new_node;
                  improved_connection = 1;
               end
            else
                bsf = bridge_node.costToDest();
               if new_node.distToHere() < bsf
                  bridge_node.leafToDest = new_node;
                  improved_connection = 1;
               end
            end
             
            if this.bothTreesHaveSolution() && improved_connection
               %sample more frequently from destinations
               current_freq = this.sampler.getDestFreq(); 
               new_freq = round((3/4)*current_freq) + 1;
               this.sampler.setDestFreq(new_freq);
            end
            
            if bridge_node.getBestPathCost() < this.BBSF.getBestPathCost()
               this.BBSF = bridge_node; 
            end
        end
        
        function tf = bothTreesHaveSolution(this)
           src_BSF = this.sourceRDT.BSF;
           dest_BSF = this.destRDT.BSF;
           tf = (src_BSF.distToHere < Inf) &&  (dest_BSF.distToHere < Inf);
        end
        
        
        
        function initializeBridgeNodes(this, pppi)
           mid_pts = pppi.getGoalRegion().goalGridPoints(); 
           source = pppi.getSourceGrid();
           dest = pppi.getSourceGrid2();
           
           map = containers.Map('KeyType','double','ValueType','any');
           
           [mp_count, ~] = size(mid_pts);

           for i=1:mp_count
              mid_pt = mid_pts(i,:);
              x = mid_pt(1); y = mid_pt(2);
              best_case_cost = norm(mid_pt - source) + norm(mid_pt - dest);
              if map.isKey(x)
                 submap = map(x);
                 %containers.Map are handle objcets, so stop complaining
                 %MATLAB
                 submap(y) = BridgeNode(mid_pt, best_case_cost);
              else
                  submap = containers.Map(y,BridgeNode(mid_pt, best_case_cost));
                  map([x]) = submap;
              end
           end
           
           this.bridgeNodes = map;
           
           %empty bridge node, will have total cost of infinity
           this.BBSF = BridgeNode([0,0],0);
           
        end
    end
end


