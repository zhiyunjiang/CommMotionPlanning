%%
x_max = 6.2;
x_min = -4.7;
y_max = 5.5;
y_min = -5.1;
region = [x_max x_min y_max y_min];
resolution = 10;

source = [5.25, -4.75];
dests = [[-3.76, 5.1]];

%obstacles
r1 = RectObs([2,-1, 0, -1.25]);
r2 = RectObs([4, 2.5, -2, -3]);
obstacles = [r1, r2];
obs_mod = ObstacleMod(obstacles);
%cost function
cost_fxn = @(n1, n2, path, mode) GridDist(path);

problem_instance = PathPlanningProblem(region, resolution, source, dests, obs_mod, cost_fxn);

%Now setup the solver
solution_not_required = 0;
max_iterations = -1;
max_run_time = 10;
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);

%1 - uniform random (RRT), 2 - deterministic sequence (RDT)
type = 1;
dest_freq = 10;
sequence = [];%sampler.sequence;%use previous run's sequence
sampler = Sampler(type, dest_freq, sequence);

%1 - RDT*, 0 - RDT
do_rewire = 1;
steer_rad = 1;
solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);

%%
solver.solve(problem_instance);
%%
hold off
clf
problem_instance.plotProb();
%%

best_path_tail = solver.getBSF();

hold on
best_path_tail.pathToRoot(1);



%%
%try a problem in continuous space
resolution2 = Inf;
goal_region_resolution = 12;
problem_instance = PathPlanningProblem(region, resolution2, source, dests, obs_mod, cost_fxn, goal_region_resolution);

%Now setup the solver
solution_not_required = 0;
max_iterations = -1;
max_run_time = 10;
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);

%1 - uniform random (RRT), 2 - deterministic sequence (RDT)
type = 4;
dest_freq = 10;
sequence = [];%sampler.sequence;%use previous run's sequence
sampler = Sampler(type, dest_freq, sequence);

%1 - RDT*, 0 - RDT
do_rewire = 1;
steer_rad = 0.5;
solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);

%%
solver.solve(problem_instance);
%%
hold off
clf
problem_instance.plotProb();
%%

best_path_tail = solver.getBSF();

hold on
best_path_tail.pathToRoot(1);