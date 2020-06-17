classdef RDTSolver < handle
    %RDTSolver - can solve instances of PathPlanningProblem using Randomly
    %Exploring Random Trees, including RRT and RRT*. Flexible enough to
    %handle some variants 
    
    properties (Access = private)
        %sampler - sampling module to use
        sampler;
        
        %stopCriteria - function handle. Function determines whether or not
        %we should stop. Arguments are elapsed time, iteration number, and cost
        %of best solution so far
        stopCriteria;
        
        %rdTree - The tree being built out
        rdTree;
        
    end
    
    methods (Access = public)
        function this = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad)
           this.sampler = sampler;
           this.stopCriteria = stop_criteria;
            
            this.rdTree = RDTree(do_rewire, steer_rad);
            
        end
        
        function solve(this, pppi)
            tic;
            this.initializeTree(pppi);
            this.sampler.setParams(pppi);
            
            iteration_count = 0;
            elapsed_time = toc;
            series_delta = 0.2;
            time_for_next_recording = series_delta;
            
            while ~this.stopCriteria.stop(elapsed_time , iteration_count, this.getBSF())
                x_rand = this.sampler.sample(iteration_count, pppi, this.minCostSoFar());
                
                this.rdTree.executeIteration(x_rand, pppi);
                
                iteration_count = iteration_count + 1;
                elapsed_time = toc;
                if elapsed_time > time_for_next_recording
                    time_for_next_recording = time_for_next_recording + series_delta;
                   this.rdTree.recordBSFCost(elapsed_time); 
                end
            end
            
        end
        
        function tail = getBSF(this)
            tail = this.rdTree.BSF;
        end
        
        function series = getBSFTimeSeries(this)
            series = this.rdTree.solCosts;
        end
        
        function cost = minCostSoFar(this)
            cost = this.getBSF().distToHere;
        end
    end
    
    methods (Access = private)

        function initializeTree(this, pppi)
           root_pos = pppi.getSourceGrid();
           %see Karaman & Frazzoli, 2011, page 29. Here, d=2 
           gridRegion = pppi.getGridRegion();
           length = (gridRegion(1) - gridRegion(2));
           width = (gridRegion(3) - gridRegion(4));
           grid_area = length*width;%total area will be >= area of X_free
           gamma =  sqrt(3*(grid_area/2*pi));
            
           is_continuous = (pppi.getStepSize() == 0);
            this.rdTree.initialize(gamma, root_pos, is_continuous);
        end
    end
end

