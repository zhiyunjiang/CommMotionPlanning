%%
%RDT take 2
x_max = 6;
x_min = -4;
y_max = 5.5;
y_min = -5;
region = [x_max x_min y_max y_min];
resolution = 10;

source = [6, -5];
dests = [[-4, 5]];

%obstacles
r1 = RectObs([2,-1, 0, -1.25]);
r2 = RectObs([4, 2.5, -2, -3]);

%obstacles = [RectObs([-1, -1.5, 3, -5]), RectObs([0, -0.5, 5.5, -2]), ...
            % RectObs([1, 0.5, 3, -5]), RectObs([2, 1.5, 5.5, -2])];
obstacles = [RectObs([-1.5, -2, 3, -5]), RectObs([0, -0.5, 5.5, -2])];
obs_mod = ObstacleMod(obstacles);
%cost function
cost_fxn = @(n1, n2, pppi) ManhattanDistance(n1, n2, pppi);

problem_instance = PathPlanningProblem(region, resolution, source, dests, obs_mod, cost_fxn);

%%
%Now setup the solver
solution_not_required = 0;
max_iterations = -1;
max_run_time = -1;
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);

%1 - uniform random (RRT), 2 - deterministic sequence (RDT)
type = 1;
dest_freq = 10;
sequence = [];%sampler.sequence;%use previous run's sequence
sampler = Sampler(type, dest_freq, sequence);

%1 - RDT*, 0 - RDT
do_rewire = 1;
steer_rad = 6;
solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);

%%
solver.solve(problem_instance);
%%
hold off
clf
problem_instance.plotProb(1);
%%

best_path_tail = solver.getBSF();

hold on
best_path_tail.pathToRoot(1);



