%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing out the old MCTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setup the stopping criteria:
solution_not_required = 1;
max_iterations = 50000;
max_run_time = 60*(20);
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);

%setup channel and channel Analyzer
[n_PL, K_PL, sh_decorr, sigma_SH, PSD_at_f_c, mp_decorr, lambda, K_ric,...
    corr_mp, sigma_mp, region, N_sin, res, n_samples, n_sims, gamma_TH] = getConstParams(); 
region = min(region, 25);
q_b = [2,7];
cp = ChannelParams(q_b, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);

cc = CommChannel(cp, N_sin, region, res);
%setup Channel Analyzer
no_mp = 0;
ca = ChannelAnalyzer(cc, gamma_TH, no_mp);

%sampler setup
%1 - uniform random (RRT), 2 - deterministic sequence (RDT) 
%3 - informed set sampling, 4 - continuous uniform random
type = 1; dest_freq = 10; sequence = [];
sampler = Sampler(type, dest_freq, sequence);


solver = MCTSSolver(stop_criteria, ca, getDPRRTSolver());
%setup path planning problem
obs_mod = ObstacleFactory.simple();
pppi = PathPlanningProblem(region, res, [25,25], [2,7], obs_mod, getCostFxnPredicted(0));

solver.solve(pppi);
%%
path = solver.getOptiamlPolicy();
%%
cc.plot( 1, 1)
hold on
plot(path(:,1), path(:,2))

function rrt_solver = getDPRRTSolver()

    %stopping criteria
    solution_not_required = 0; max_iterations = Inf; mins = 0; max_run_time = mins*60;
    stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);
    %sampler setup
    %1 - uniform random (RRT), 2 - deterministic sequence (RDT) 
    %3 - informed set sampling, 4 - continuous uniform random
    type = 1; dest_freq = 10;
    sampler = Sampler(type, dest_freq, []);

    %0 - RDT, 1 - RDT*
    do_rewire = 0;steer_rad = 6;
    rrt_solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);
end