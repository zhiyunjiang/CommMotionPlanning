%Empirical Evidence of Neat Little Bound
%Solving Path Planning Problems for variable transmit power case
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given channel params and a handful of readings, find
% "communication rich" path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Connectivity Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position of the base station (remote station or transmitter)
%q_b = [3 5.5];
q_b = [2 7];
% Path loss parameters
% The resulting path loss component in dB is 
% gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
% K_PL is the path loss constant in dB and n_PL is the path loss exponent.\
% n_PL = 3.86;
n_PL = 2;
% K_PL = -41.34;
K_PL = -40;


%Shadowing Parameters
%Shadowing decorrelation distance
sh_decorr = 3.09;
%shadowing power
alpha = 10.24;
%standard deviation of shadowing and power conn
sigma_SH = sqrt(alpha);
PSD_at_f_c = 30;

%Multipath Parameters
%Multipath decorrelation, used to set grid resolution
%not relevant if modeling mp as uncorrelated
mp_decorr = 0.05;         

lambda = 0.125;
%K_ric = 14.56;     
K_ric = 50;

% corr_mp = 1 -> correlated multipath 
% corr_mp = 0 -> uncorrelated multipath 
corr_mp = 0;

%
sigma_mp = 10;%=0 -> no mp when modeling MP fading as zero-mean, lognormal

cp = ChannelParams(q_b, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Comm Channel Simulation & Generation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_max = 50;
x_min = 0;
y_max = 50;
y_min = 0;
region = [x_max x_min y_max y_min];
%simulation parameter
N_sin = 5000;

%we're not going to be looking at MP component for a quick minute, so
%setting the resolution to something more meaningful
%res = 2/mp_decorr;
res = 5;
cc = CommChannel(cp, N_sin, region, res);
cc.simulateSH(); cc.simulateMP();
%%


%req_power = cc.getReqTXPowermW(receiver_noise, BER, R) + cc_poi.getReqTXPowermW(receiver_noise, BER, R);
 
%cc.plotField(10*log10(req_power), 'Required TX Power for Relay Communication', 'Power (dBm)');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path Plan to Connected Spots Using Expected distance cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup obstacles
source = [5, 45];
goal = [49, 10];
%%
maze = [RectObs([7, 0, 35, 15]), RectObs([23, 16, 50, 27]), ...
            RectObs([38, 30, 30, 0]), RectObs([50, 42, 35, 25])];
radpmax = 2;
rbase = 2.5;
posrng = [-5,5];
%pruned from grid: CircObs(randi(radpmax)+rbase, [10,10] + randi(posrng,[1,2])), CircObs(randi(radpmax)+rbase, [10, 25] + randi(posrng,[1,2])),...
circle_grid = [CircObs(randi(radpmax)+rbase, [10,40] + randi(posrng,[1,2])), CircObs(randi(radpmax)+rbase, [25, 10] + randi(posrng,[1,2])),...
    CircObs(randi(radpmax)+rbase, [25,25] + randi(posrng,[1,2])), CircObs(randi(radpmax)+rbase, [25, 40] + randi(posrng,[1,2])),...
    CircObs(randi(radpmax)+rbase, [40,10] + randi(posrng,[1,2])), CircObs(randi(radpmax)+rbase, [40, 25] + randi(posrng,[1,2])),...
    CircObs(randi(radpmax)+rbase, [40,40] + randi(posrng,[1,2]))];
simple = [CircObs(3,[34,10]), CircObs(3,[37,44]), CircObs(3,[34,30])]; 
no_obs = [];

rand_obs = [CircObs(randi(3)+2, [randi(51)-1, randi(51)-1]), CircObs(randi(3)+2, [randi(51)-1, randi(51)-1]),...
            CircObs(randi(3)+2, [randi(51)-1, randi(51)-1]), CircObs(randi(3)+2, [randi(51)-1, randi(51)-1]),...
            CircObs(randi(3)+2, [randi(51)-1, randi(51)-1]), CircObs(randi(3)+2, [randi(51)-1, randi(51)-1])];
obs_mod = ObstacleMod(no_obs);
%%

cost_fxn = @(n1, n2, path, mode) GridDist(n1,n2, path, mode);
path_res = res;
problem_instance = PathPlanningProblem(region, path_res, source, goal, obs_mod, cost_fxn);

hold on
problem_instance.plotProb(0)
hold on
scatter(q_b(1), q_b(2), 'v', 'filled');
%%
%(1) find shortest path using Dijkstra's
astar_solver = AStarSolver();
tic
astar_solver.solve(problem_instance);
toc
bst = astar_solver.BST;
astar_path = bst.pathToRoot(0);
d_min = PathDistance(astar_path);
%%
%solve using RRT
%stop on first solution, unless a solution hasn't
%been found
num_splits = 5;
theoretical_bounds = zeros([2,num_splits]);
empirical_rel = zeros([2,num_splits]);
for i = 1:num_splits
    
    K_pl_new = K_PL + (i-1)*10;
    cc.setKPl(K_pl_new);
    ca = ChannelAnalyzer(cc, -90, 0);
    req_tx_power = ca.getReqTXPower(receiver_noise, R, K);
    min_power = min(req_tx_power,[],'all')/1000; %to convert to Joules
    max_power = max(req_tx_power,[],'all')/1000;% to convert to Joules
    
    power_ratio = (max_power - min_power)/(motion_power + min_power);
    theoretical_bounds(:,i) = [power_ratio, power_ratio*d_min];
    
    
    receiver_noise = 1e-10; R = 6; BER = 1e-6; K = -1.5/log(5*BER);
    k_1 = 2.5; k_2 = 0.29; v_const = (1)*cc.res; %1 m/s
    motion_power = (k_1 + (k_2/v_const));
    cost_fxn = @(n1, n2, path, mode) LITotalEnergy(n1, n2, path, cc, receiver_noise, R, K, motion_power);
    path_res = res;
    problem_instance = PathPlanningProblem(region, path_res, source, goal, obs_mod, cost_fxn);
    
    astar_solver = AStarSolver();
    tic
    astar_solver.solve(problem_instance);
    toc
    bst = astar_solver.BST;
    astar_path = bst.pathToRoot(0);
    dist = PathDistance(astar_path);
    empirical_rel(:,i) = [power_ratio, dist - d_min]; 
end
%%
plot(empirical_rel(1,:), empirical_rel(2,:), 'r');
hold on
plot(theoretical_bounds(1,:), theoretical_bounds(2,:), 'g');