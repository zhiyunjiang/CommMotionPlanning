%Deterministic Comm Path Planning
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_max = 6;
x_min = -4;
y_max = 5.5;
y_min = -5;
region = [x_max x_min y_max y_min];

source = [-4, 5];
dest = [6, -5];

maze = [RectObs([-1, -1.5, 3, -5.5]), RectObs([0, -0.5, 6, -2]), ...
            RectObs([1, 0.5, 3, -5.5]), RectObs([2, 1.5, 6, -2])];
simple = [CircObs(1.5, [0,4]), RectObs([4, 2, 0, -2 ]), CircObs(1, [2,2])];

dense = RectObs.empty;
dx = 0.25;
dy = 0.25;
for i=1:floor(x_max - x_min)
    
    for j = 1:floor(y_max - y_min)
        y_jitter = 2*dy*rand(1); 
        x_jitter = 2*dx*rand(1); 
        obs = RectObs([i+x_min+dx + x_jitter, i+x_min + x_jitter, j+y_min+dy + y_jitter, j+y_min+y_jitter]);
        dense((i-1)*floor(y_max - y_min) + j) = obs;
    end 
end

obs_mod = ObstacleMod(dense);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Connectivity Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position of the base station (remote station or transmitter)
%q_b = [3 5.5];
q_b = dest;
% Path loss parameters
% The resulting path loss component in dB is 
% gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
% K_PL is the path loss constant in dB and n_PL is the path loss exponent.
n_PL = 3;  
K_PL = -12.89;
%Shadowing decorrelation distance
beta = 2;    

%Multipath decorrelation, used to set grid resolution
ss_decorr = 0.05;         
res = 2/ss_decorr;

%shadowing power
alpha = 10;
%standard deviation of shadowing and power conn
sigma_SH = sqrt(alpha);

N_sin = 5000;                    
PSD_at_f_c = 30;
lambda = 0.125;
K_ric = 10;     

% corr_mp = 1 -> correlated multipath 
% corr_mp = 0 -> uncorrelated multipath 
corr_mp = 0;

[gamma_TOT_dB, gamma_PL_SH_dB, gamma_PL_dB, gamma_SH_dB, gamma_MP_LIN, g_x, g_y] = channel_simulator(region, ...
    q_b, K_PL, n_PL, alpha, beta, N_sin, PSD_at_f_c, lambda, K_ric, res, corr_mp);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Set of Connected Points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Channel power connectiivty threshold (dB)
gamma_TH = -40;
sz = size(gamma_TOT_dB);
connected = [];
for i = 1:sz(1)
    y_loc = (i-1)/res + y_min;
    for j = 1:sz(2)
        if gamma_TOT_dB(i, j) >= gamma_TH
            %let problem instance handle pruning
            x_loc = (j-1)/res + x_min;
            connected = [connected; [x_loc, y_loc]];
        end
    end
end

%scatter(connected(:,1), connected(:,2)); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path Plan to Connected Spots Using Manhattan distance cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cost function
cost_fxn = @(n1, n2, path) ManhattanDistance(n1, n2, path);
path_res = 10;
step_size = 1/path_res;
problem_instance = PathPlanningProblem(region, path_res, source, connected, obs_mod, cost_fxn);
%%
figure (1)
clf
hold on
problem_instance.plotProb(0)
%%
cost_fxn_continuous = @(n1, n2, path) norm(n1.getPos() - n2.getPos());
problem_instance_continuous = PathPlanningProblem(region, Inf, source, connected, obs_mod, cost_fxn_continuous, res);

%%
%solve using RRT*

%stop on first solution, stop after 2 seconds, unless a solution hasn't
%been found
solution_not_required = 0;
max_iterations = Inf;
%max run time in seconds
max_run_time = 20;
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);

%1 - uniform random (RRT), 2 - deterministic sequence (RDT) 
%3 - informed set sampling, 4 - continuous uniform random
type = 1;
dest_freq = 20;
sequence = [];%sampler.sequence;%use previous run's sequence
sampler = Sampler(type, dest_freq, sequence);


do_rewire = 1;
steer_rad = 5;

rrt_solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);
tic
rrt_solver.solve(problem_instance);
toc
%%
figure(1)
hold on
bsf = rrt_solver.getBSF();
bsf.pathToRoot(1, step_size);

figure(2)
clf
series = rrt_solver.getBSFTimeSeries();
plot(series(:,1), step_size*series(:,2));
%%
%Solve with Dijkstras
dijkstraSolver = AStarSolver();
tic
dijkstraSolver.solve(problem_instance);
toc
%%
bst = dijkstraSolver.BST;
figure(1)
hold on
bst.pathToRoot(1, step_size);
%%
sl_heuristic = @(pt, inst) SLtoGoal(pt, inst);
mh_heuristic = @(pt, inst) ManhattantoGoal(pt, inst);
astar_solver = AStarSolver(sl_heuristic);
tic
astar_solver.solve(problem_instance);
toc
%%
bst = astar_solver.BST;
figure(1)
hold on
bst.pathToRoot(1, step_size);
legend('RRT*', 'Dijkstra''s', 'A*');
%%
%The conintuous spce problem
solution_not_required = 0;
max_iterations = Inf;
%max run time in seconds
max_run_time = 30;
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);

%1 - uniform random (RRT), 2 - deterministic sequence (RDT) 
%3 - informed set sampling, 4 - continuous uniform random
type = 4;
dest_freq = 20;
sequence = [];
sampler = Sampler(type, dest_freq, sequence);


do_rewire = 1;
steer_rad = 0.5;

cont_rrt_solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);
tic
cont_rrt_solver.solve(problem_instance_continuous);
toc
%%
figure(1)
hold on
bsf = cont_rrt_solver.getBSF();
bsf.pathToRoot(1);

figure(2)
hold on
series = cont_rrt_solver.getBSFTimeSeries();
plot(series(:,1), series(:,2));
