%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Connectivity Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position of the base station (remote station or transmitter)
%q_b = [3 5.5];
q_b = [-7 -7];
% Path loss parameters
% The resulting path loss component in dB is 
% gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
% K_PL is the path loss constant in dB and n_PL is the path loss exponent.
n_PL = 6;  
K_PL = 0;

%Shadowing Parameters
%Shadowing decorrelation distance
sh_decorr = 1;    
%shadowing power
alpha = 20;
%standard deviation of shadowing and power conn
sigma_SH = sqrt(alpha);
PSD_at_f_c = 30;

%Multipath Parameters
%Multipath decorrelation, used to set grid resolution
mp_decorr = 0.05;         

lambda = 0.125;
K_ric = 10;     

% corr_mp = 1 -> correlated multipath 
% corr_mp = 0 -> uncorrelated multipath 
corr_mp = 0;

cp = ChannelParams(q_b, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Comm Channel Simulation & Generation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_max = 5;
x_min = -5;
y_max = 5;
y_min = -5;
region = [x_max x_min y_max y_min];
%simulation parameter
N_sin = 5000;

%we're not going to be looking at MP component for a quick minute, so
%setting the resolution to something more meaningful
%res = 2/mp_decorr;
res = 16;
cc = CommChannel(cp, N_sin, region, res);
%cc.simulateShadowing();
cc.plot()
%%
%minimum power required for connectivity
gamma_TH = -40;
no_mp = 1;
ca = ChannelAnalyzer(cc, gamma_TH, no_mp);
figure(2)
clf
ca.plot_priors()

%%
%mininum probability of channel power being greater than or equal to
%gamma_TH require to count spot as a connected spot
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
source = [0, 0];

maze = [RectObs([-1, -1.5, 3, -5.5]), RectObs([0, -0.5, 6, -2]), ...
            RectObs([1, 0.5, 3, -5.5]), RectObs([2, 1.5, 6, -2])];
simple = [CircObs(0.25, [-1,-1]), CircObs(0.25, [-1, -2]),...
    CircObs(0.25, [-1,-3]), CircObs(0.25, [-1, -4]),...
    CircObs(0.25, [-3,-1.5]), CircObs(0.25, [-3, -2.5]),...
    CircObs(0.25, [-3,-3.5]), CircObs(0.25, [-3, -4.5]),...
    CircObs(0.25, [-2,-1.5]), CircObs(0.25, [-2, -2.5]),...
    CircObs(0.25, [-2,-3.5]), CircObs(0.25, [-2, -0.5])];

obs_mod = ObstacleMod(simple);
%cost function
cost_fxn = @(n1, n2, path) ExpectedDistance(n1, n2, ca, 2);
path_res = res;
problem_instance = PathPlanningProblem(region, path_res, source, connected, obs_mod, cost_fxn);

figure (3)
clf
hold on
problem_instance.plotProb(1)
%%
%solve using RRT*

%stop on first solution, stop after 20 seconds, unless a solution hasn't
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
steer_rad = 1;

rrt_solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);
tic
rrt_solver.solve(problem_instance);
toc
%%
figure(3)
hold on
bsf = rrt_solver.getBSF();
path = bsf.pathToRoot(1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify that metric is working as expected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
root = [50,50];
qB = cc.getGridBase();
dsrc = norm(root - qB);
thetaos = atan2(-qB(2),-qB(1));
%calculation only works when root is to the right of qB
thetasrc = 0;
step = sqrt(2)/10;
num_steps = 49*10;

pathx = cumsum(ones([50,1]));
pathy = cumsum(ones([50,1]));

path = [flip(pathx), flip(pathy)];
%calculate the expected distance and CMF (CDF)
% tic
% exp_dist_M = ca.ExpectedFPD(path,2)
% toc
% % tic
% % exp_dist_B = ca.ExpectedFPD(path,1)
% % toc
% tic
% exp_dist_A = ca.ExpectedFPD(path,3)
% toc
%[approx_PMF_B, distances] = ca.ApproxFPDPMF(path,1);
% [approx_PMF_M, distances] = ca.ApproxFPDPMF(path,2);
% fprintf('Found PMF M');
% [approx_PMF_A, ~] = ca.ApproxFPDPMF(path,3);
%fprintf('Found PMF A');
[approx_PMF_P, pd] = ca.ApproxFPDPMF2(step, num_steps, dsrc, thetasrc, thetaos, root);
fprintf('Found PMF P');
%approx_CDF_B = cumsum(approx_PMF_B);
% approx_CDF_M = cumsum(approx_PMF_M);
% approx_CDF_A = cumsum(approx_PMF_A);
approx_CDF_P = cumsum(approx_PMF_P);
%%
step_size = 1/res;
figure(5)
clf;
%plot(cumsum(distances)*step_size, approx_CDF_B);

hold on
%plot(cumsum(distances)*step_size, approx_CDF_M);
%plot(cumsum(distances)*step_size, approx_CDF_A);
plot(cumsum(pd)*step_size, approx_CDF_P);

%%
%verify by simulating several channels and calculating expected distance
%to connectivity

n_sims = 200;
sims_run = 0;
all_dists = zeros([n_sims, 1]);

while sims_run < n_sims
   %generate a new shadowing component
   cc.simulateShadowing();
   gamma_TOT = cc.getGammaTOTdB();
   %check if root is not connected (event we're conditioned on)
   root_gamma = cc.getGammaTOTdBAtPoint(path(1,:));
   if root_gamma > gamma_TH
       %if not, we're not interested, simulate again.
       continue;
   end
   
   current_dist = 0;
   for j=1:length(path)-1
      gamma_here = cc.getGammaTOTdBAtPoint(path(j,:));
      if gamma_here < gamma_TH
        current_dist = current_dist + norm(path(j,:) - path(j+1,:));
      else
          break;
      end
   end
   sims_run = sims_run + 1;
   all_dists(sims_run) = current_dist;
end
histogram(all_dists*step_size, 'Normalization', 'cdf', 'BinEdges',(cumsum(distances))*step_size);
%legend('Hist','MATLAB', 'A', 'P', 'Hist');
avg_dist = sum(all_dists)/n_sims

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
