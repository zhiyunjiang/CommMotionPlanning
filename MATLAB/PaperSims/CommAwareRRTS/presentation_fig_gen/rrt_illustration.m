%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4
% For pre-defined workspace, obstacles, channel params, etc... runs
% multiple sims to find average costs. This is for variable transmission
% power, all tasks (uploading, broadcasting, and relaying).
% 
% Fig 4(left) - Upload
% Fig 4(middle) - Broadcast
% Fig 4(right) - Relay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fetch predifned params
[n_PL, K_PL, sh_decorr, sigma_SH, PSD_at_f_c, mp_decorr, lambda, K_ric,...
    corr_mp, sigma_mp, region, N_sin, res, n_samples, n_sims, gamma_TH] = getConstParams(); 

res = 10;
%Workspace params
obs_mod = ObstacleFactory.maze();
source = [1, 43]; goal = [49, 17];

pppi = PathPlanningProblem(region, res, source, goal, obs_mod, getCostFxnPredicted(0));


%stopping criteria
solution_not_required = 0; max_iterations = 5000; mins = 20; max_run_time = mins*60;
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);
%sampler setup
%1 - uniform random (RRT), 2 - deterministic sequence (RDT) 
%3 - informed set sampling, 4 - continuous uniform random
type = 1; dest_freq = 100; sequence = [];%sampler.sequence;%use previous run's sequence
sampler = Sampler(type, dest_freq, sequence);

%0 - RDT, 1 - RDT*
do_rewire = 1;steer_rad = 3;
rrt_solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);


rrt_solver.solve(pppi);
bsf = rrt_solver.getBSF();
rrt_path = bsf.pathToRoot(0);
