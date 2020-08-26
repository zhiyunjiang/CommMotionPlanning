%Example of how to use MCTSSolver
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Connectivity Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position of the base station (remote station or transmitter)
%q_b = [3 5.5];
q_b = [-25 -25];
% Path loss parameters
% The resulting path loss component in dB is 
% gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
% K_PL is the path loss constant in dB and n_PL is the path loss exponent.
n_PL = 4.2;  
K_PL = 27;

%Shadowing Parameters
%Shadowing decorrelation distance
sh_decorr = 2;    
%shadowing power
alpha = 8.41;
%standard deviation of shadowing and power conn
sigma_SH = sqrt(alpha);
PSD_at_f_c = 30;

%Multipath Parameters
%Multipath decorrelation, used to set grid resolution
%not relevant if modeling mp as uncorrelated
mp_decorr = 0.05;         

lambda = 0.125;
K_ric = 1.59;     

% corr_mp = 1 -> correlated multipath 
% corr_mp = 0 -> uncorrelated multipath 
corr_mp = 0;

cp = ChannelParams(q_b, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Comm Channel Simulation & Generation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_max = 25;
x_min = 0;
y_max = 25;
y_min = 0;
region = [x_max x_min y_max y_min];
%simulation parameter
N_sin = 5000;

%we're not going to be looking at MP component for a quick minute, so
%setting the resolution to something more meaningful
%res = 2/mp_decorr;
res = 1;
cc = CommChannel(cp, N_sin, region, res);
cc.plot(1)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the Channel Analyzer, Preview Prior Probability of Connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%minimum power required for connectivity
gamma_TH = -44;
no_mp = 0;
ca = ChannelAnalyzer(cc, gamma_TH, no_mp);
figure(2)
clf
ca.plot_priors()

%%
% %mininum probability of channel power being greater than or equal to
% %gamma_TH require to count spot as a connected spot
p_TH = 0.95;

gamma_diff = gamma_TH - cc.getGammaPLdB();% what will need to be accounted for with shadowing.
%find min value
min_diff = norminv(1 - p_TH,0,sigma_SH);
sz = size(gamma_diff);
connected = [];
for i = 1:sz(1)
    y_loc = (i-1)/res + y_min;
    for j = 1:sz(2)
        if gamma_diff(i, j) <= min_diff
            x_loc = (j-1)/res + x_min;
            connected = [connected; [x_loc, y_loc]];
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path Plan to Connected Spots Using Expected distance cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup obstacles
source = [20, 20];

maze = [RectObs([-1, -1.5, 3, -5.5]), RectObs([0, -0.5, 6, -2]), ...
            RectObs([1, 0.5, 3, -5.5]), RectObs([2, 1.5, 6, -2])];
simple = [CircObs(1.5, [5,5]), CircObs(1.5, [5, 10]),...
    CircObs(1.5, [5,15]), CircObs(1.5, [5, 20]),...
    CircObs(1.5, [10,7.5]), CircObs(1.5, [10, 12.5]),...
    CircObs(1.5, [10,17.5]), CircObs(1.5, [10, 22.5]),...
    CircObs(1.5, [15, 5]), CircObs(1.5, [15, 10]),...
    CircObs(1.5, [15,15]), CircObs(1.5, [15, 20])];

obs_mod = ObstacleMod(simple);
%cost function - this will only be used in RRT simulations
cost_fxn = @(n1, n2, path, mode) GridDist(n1, n2, path, mode);
path_res = res;
problem_instance = PathPlanningProblem(region, path_res, source, connected, obs_mod, cost_fxn);
%%
figure (3)
clf
hold on
problem_instance.plotProb(1)
%%
%solve using MCTS
sims_per_stage = 50;
MCTS_start = cputime;
[root, leaf] = MCTSSolver.solve(problem_instance, ca,  25*25*res^2, sims_per_stage);
MCTS_time = cputime - MCTS_start
%%
node = leaf;
path = node.getPos();
while ~node.isRoot
   node = node.parent;
   path = [node.getPos();path];
end
hold on
plot(path(:,1), path(:,2),'o-');
exp_dist = ca.ExpectedFPD(path,4)
%%
%solve using RRT

%stop on first solution, unless a solution hasn't
%been found
solution_not_required = 0;
max_iterations = Inf;
%max run time in seconds
max_run_time = 0;
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);

%1 - uniform random (RRT), 2 - deterministic sequence (RDT) 
%3 - informed set sampling, 4 - continuous uniform random
type = 1;
dest_freq = 20;
sequence = [];%sampler.sequence;%use previous run's sequence
sampler = Sampler(type, dest_freq, sequence);


do_rewire = 0;
steer_rad = 5;

rrt_solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);


tic
run_time = 0;
best_cost = Inf;
best_path = [];
cost_sum = 0;
rrt_counts = 0;
while MCTS_time > run_time
   rrt_counts = rrt_counts + 1;
   run_time = run_time + toc;
   rrt_solver.solve(problem_instance);
   bsf = rrt_solver.getBSF();
   rrt_path = bsf.pathToRoot(0);
   cost = ca.ExpectedFPD(rrt_path,4);
   cost_sum = cost_sum + cost;
   if cost < best_cost
       best_path = rrt_path;
       best_cost = cost;
   end
end
%%
hold on
plot(best_path(:,1), best_path(:,2), '+-');
best_cost 
avg_cost = cost_sum/rrt_counts

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some alternate methods for solving the path planning problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Solve with Dijkstras
% dijkstraSolver = AStarSolver();
% tic
% dijkstraSolver.solve(problem_instance);
% toc
% %%
% bst = dijkstraSolver.BST;
% figure(3)
% hold on
% bst.pathToRoot(1);
% %%
% sl_heuristic = @(pt, inst) SLtoGoal(pt, inst);
% mh_heuristic = @(pt, inst) ManhattantoGoal(pt, inst);
% astar_solver = AStarSolver(sl_heuristic);
% tic
% astar_solver.solve(problem_instance);
% toc
% %%
% bst = astar_solver.BST;
% figure(1)
% hold on
% bst.pathToRoot(1);
% legend('RRT*', 'Dijkstra''s', 'A*');
