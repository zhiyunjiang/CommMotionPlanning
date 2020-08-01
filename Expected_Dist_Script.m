%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Connectivity Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position of the base station (remote station or transmitter)
%q_b = [3 5.5];
q_b = [-450 0];
% Path loss parameters
% The resulting path loss component in dB is 
% gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
% K_PL is the path loss constant in dB and n_PL is the path loss exponent.
n_PL = 4.2;  
K_PL = 27;

%Shadowing Parameters
%Shadowing decorrelation distance
sh_decorr = 12.92;    
%shadowing power
alpha = 8.41;
%standard deviation of shadowing and power conn
sigma_SH = sqrt(alpha);
PSD_at_f_c = 30;

%Multipath Parameters
%Multipath decorrelation, used to set grid resolution
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
x_max = 100;
x_min = 0;
y_max = 100;
y_min = 0;
region = [x_max x_min y_max y_min];
%simulation parameter
N_sin = 5000;

%we're not going to be looking at MP component for a quick minute, so
%setting the resolution to something more meaningful
%res = 2/mp_decorr;
res = 10;
cc = CommChannel(cp, N_sin, region, res);
%cc.simulateShadowing();
cc.plot()
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the Channel Analyzer, Preview Prior Probability of Connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%minimum power required for connectivity
gamma_TH = -80;
no_mp = 0;
ca = ChannelAnalyzer(cc, gamma_TH, no_mp);
% figure(2)
% clf
% ca.plot_priors()

%%
% %mininum probability of channel power being greater than or equal to
% %gamma_TH require to count spot as a connected spot
% p_TH = 0.95;
% 
% gamma_diff = gamma_TH - cc.getGammaPLdB();% what will need to be accounted for with shadowing.
% %find min value
% min_diff = norminv(1 - p_TH,0,sigma_SH);
% sz = size(gamma_diff);
% connected = [];
% for i = 1:sz(1)
%     y_loc = (i-1)/res + y_min;
%     for j = 1:sz(2)
%         if gamma_diff(i, j) <= min_diff
%             x_loc = (j-1)/res + x_min;
%             connected = [connected; [x_loc, y_loc]];
%         end
%     end
% end
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Path Plan to Connected Spots Using Expected distance cost function
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %setup obstacles
% source = [0, 0];
% 
% maze = [RectObs([-1, -1.5, 3, -5.5]), RectObs([0, -0.5, 6, -2]), ...
%             RectObs([1, 0.5, 3, -5.5]), RectObs([2, 1.5, 6, -2])];
% simple = [CircObs(0.25, [-1,-1]), CircObs(0.25, [-1, -2]),...
%     CircObs(0.25, [-1,-3]), CircObs(0.25, [-1, -4]),...
%     CircObs(0.25, [-3,-1.5]), CircObs(0.25, [-3, -2.5]),...
%     CircObs(0.25, [-3,-3.5]), CircObs(0.25, [-3, -4.5]),...
%     CircObs(0.25, [-2,-1.5]), CircObs(0.25, [-2, -2.5]),...
%     CircObs(0.25, [-2,-3.5]), CircObs(0.25, [-2, -0.5])];
% 
% obs_mod = ObstacleMod(simple);
% %cost function
% cost_fxn = @(n1, n2, path) ExpectedDistance(n1, n2, ca, 2);
% path_res = res;
% problem_instance = PathPlanningProblem(region, path_res, source, connected, obs_mod, cost_fxn);
% 
% figure (3)
% clf
% hold on
% problem_instance.plotProb(1)
% %%
% %solve using RRT*
% 
% %stop on first solution, stop after 20 seconds, unless a solution hasn't
% %been found
% solution_not_required = 0;
% max_iterations = Inf;
% %max run time in seconds
% max_run_time = 20;
% stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);
% 
% %1 - uniform random (RRT), 2 - deterministic sequence (RDT) 
% %3 - informed set sampling, 4 - continuous uniform random
% type = 1;
% dest_freq = 20;
% sequence = [];%sampler.sequence;%use previous run's sequence
% sampler = Sampler(type, dest_freq, sequence);
% 
% 
% do_rewire = 1;
% steer_rad = 1;
% 
% rrt_solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);
% tic
% rrt_solver.solve(problem_instance);
% toc
% %%
% figure(3)
% hold on
% bsf = rrt_solver.getBSF();
% path = bsf.pathToRoot(1);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify that metric is working as expected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%root = [0, 0];
% qB = cc.getGridBase();
% dsrc = norm(root - qB);
% thetaos = atan2(-qB(2),-qB(1));
% %calculation only works when root is to the right of qB
% thetaphi = 5*pi/4;
% step = 0.1; num_steps = 1000;

pathy = zeros([100*res, 1]);
pathx = cumsum(ones([100*res, 1]));

path = [flip(pathx), pathy];

% fscale = 10;
% fpathx = cumsuones([70*fscale,1]))/fscale;
% fpathy = cumsum(ones([70*fscale,1]))/fscale;
% 
% fpath = [fpathx, fpathy];
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
%[approx_PMF_M, distances] = ca.ApproxFPDPMF(path,2);
%fprintf('Found PMF M\n');
% [approx_PMF_Pc, pdc] = ca.ApproxFPDPMF2(step, num_steps, dsrc, thetasrc, thetaos, root);
% fprintf('Found PMF P Coarse\n');
% [approx_PMF_Pf, pdf] = ca.ApproxFPDPMF2(step/10, num_steps*10, dsrc, thetasrc, thetaos, root);
% fprintf('Found PMF P Fine\n');
% [approx_PMF_Pef, pdef] = ca.ApproxFPDPMF2(step/20, num_steps*20, dsrc, thetasrc, thetaos, root);
% fprintf('Found PMF P Extra Fine\n');
% [approx_PMF_A, distances] = ca.ApproxFPDPMF(path,3);
% fprintf('Found PMF A\n');
% [approx_PMF_Af, fdistances] = ca.ApproxFPDPMF(fpath,3);
% fprintf('Found PMF A Fine\n');
[approx_PMF_4, distances] = ca.ApproxFPDPMF(path,4);
fprintf('Found PMF 4\n');
%%
%approx_CDF_M = cumsum(approx_PMF_M);
% approx_CDF_Pc = cumsum(approx_PMF_Pc);
% approx_CDF_Pf = cumsum(approx_PMF_Pf);
% approx_CDF_Pef = cumsum(approx_PMF_Pef);
% approx_CDF_A = cumsum(approx_PMF_A);
% approx_CDF_Af = cumsum(approx_PMF_Af);
approx_CDF_4 = cumsum(approx_PMF_4);

step_size = 1/res;
figure(5)
clf;
%plot(cumsum(distances)*step_size, approx_CDF_B);

hold on
%plot(cumsum(distances)*step_size, approx_CDF_M);
% plot(cumsum(distances)*step_size, approx_CDF_A);
% plot(cumsum(fdistances)*step_size, approx_CDF_Af);
% plot(cumsum(pdc)*step_size, approx_CDF_Pc);
% plot(cumsum(pdf)*step_size, approx_CDF_Pf);
% plot(cumsum(pdef)*step_size, approx_CDF_Pef);
plot(cumsum(distances)*step_size, approx_CDF_4);

%verify by simulating several channels and calculating expected distance
%to connectivity

all_dists = ca.simulateFPD(path, 1000, 1);
histogram(all_dists*step_size, 'Normalization', 'cdf', 'BinEdges',(cumsum(distances))*step_size);
legend('No MP', 'sims');

ylabel('Probability')
xlabel('Distance (m)')
title('FPD CDF for Straight Line Path')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testing with non-straight path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stair Step Path
% root = [32,32];
% deltax = ones([62, 1]);deltax(1:2:end) = 0;
% deltay = ones([62, 1]);deltay(2:2:end) = 0;
% pathx = 32 - cumsum(deltax);
% pathy = 32 - cumsum(deltay);
% path = [root;pathx, pathy];
% 
% fdeltax = ones([62*10, 1])/10;fdeltax(1:2:end) = 0;
% fdeltay = ones([62*10, 1])/10;fdeltay(2:2:end) = 0;
% fpathx = 32 - cumsum(fdeltax);
% fpathy = 32 - cumsum(fdeltay);
% fpath = [root; fpathx, fpathy];

root = [150, 1];
deltax = ones([149, 1]);%deltax(1:2:end) = 0;
deltay = zeros([149, 1]);deltay(2:3:end) = 1;
pathx = 150 - cumsum(deltax);
pathy = 1 + cumsum(deltay);
path = [root;pathx, pathy];

fdeltax = ones([149*10, 1])/10;
fdeltay = zeros([149*10, 1])/10;fdeltay(2:3:end) = 1;
fpathx = 150 - cumsum(fdeltax);
fpathy = 1 + cumsum(fdeltay);
fpath = [root; fpathx, fpathy];

% [approx_PMF_M, distances] = ca.ApproxFPDPMF(path,2);
% fprintf('Found PMF M\n');
% approx_CDF_M = cumsum(approx_PMF_M);

[approx_PMF_A, distances] = ca.ApproxFPDPMF(path,3);
fprintf('Found PMF A\n');
approx_CDF_A = cumsum(approx_PMF_A);

[approx_PMF_Af, fdistances] = ca.ApproxFPDPMF(fpath,3);
fprintf('Found PMF A Fine\n');
approx_CDF_Af = cumsum(approx_PMF_Af);
%%
step_size = 1/res;
figure(6)
clf;
hold on
%plot(cumsum(distances)*step_size, approx_CDF_M);
plot(cumsum(distances)*step_size, approx_CDF_A);
plot(cumsum(fdistances)*step_size, approx_CDF_Af);
%verify by simulating several channels and calculating expected distance
%to connectivity

all_dists = ca.simulateFPD(path, 200);
histogram(all_dists*step_size, 'Normalization', 'cdf', 'BinEdges',(cumsum(distances))*step_size);
legend( 'Point-based', 'Point-based fine', 'Sims');
ylabel('Probability')
xlabel('Distance (m)')
title('FPD CDF for Stair-Stepping Path')
%avg_dist = sum(all_dists)/n_sims
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
