%%
% test DP solver
x_max = 6;
x_min = -4;
y_max = 5.5;
y_min = -5;
region = [x_max x_min y_max y_min];
resolution = 3;

source = [6, -5];
dests = [[-4, 5]];

%obstacles
r1 = RectObs([2,-1, 0, -1.25]);
r2 = RectObs([4, 2.5, -2, -3]);

obstacles = [RectObs([-1, -1.5, 3, -5]), RectObs([0, -0.5, 5.5, -2]), ...
             RectObs([1, 0.5, 3, -5]), RectObs([2, 1.5, 5.5, -2])];
obs_mod = ObstacleMod(obstacles);
%cost function
cost_fxn = @(n1, n2) ManhatanDistance(n1, n2);

problem_instance = PathPlanningProblem(region, resolution, source, dests, obs_mod, cost_fxn);

%%
dp = DPSolver();
dp.solve(problem_instance);