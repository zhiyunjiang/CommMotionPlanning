classdef StopCriteria < handle
    
    properties
        solutionNotRequired;
        maxIterations;
        maxRunTime;
    end
    
    methods
        function o = StopCriteria(solution_not_required, max_iterations, max_run_time)
            o.solutionNotRequired = solution_not_required;
            o.maxIterations = max_iterations;
            o.maxRunTime = max_run_time;
        end
        
        function do_stop = stop(o, elapsed_time, iteration_count, BSF)
            do_stop = 0;
            if (BSF.distToHere < Inf) || (o.solutionNotRequired)
               if iteration_count > o.maxIterations || elapsed_time > o.maxRunTime
                  do_stop = 1; 
               end
            end
        end
    end
end

