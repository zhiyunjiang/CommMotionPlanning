%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getRRTSolver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates RRTSolver objects to be used in simulations
% Output:
% rrt_solver - the rrt_solver object used to solve a path planning problem
%               instance


function rrt_solver = getRRTSolver()

    %stopping criteria
    solution_not_required = 0; max_iterations = Inf; mins = 3; max_run_time = mins*60;
    stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);
    %sampler setup
    %1 - uniform random (RRT), 2 - deterministic sequence (RDT) 
    %3 - informed set sampling, 4 - continuous uniform random
    type = 1; dest_freq = 100; sequence = [];%sampler.sequence;%use previous run's sequence
    sampler = Sampler(type, dest_freq, sequence);

    %0 - RDT, 1 - RDT*
    do_rewire = 1;steer_rad = 2;
    rrt_solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);
end