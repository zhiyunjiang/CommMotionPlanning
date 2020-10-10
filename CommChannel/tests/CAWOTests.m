%CAWithObservationsTests
%%
%Setup
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
sh_decorr = 6;    
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

%
sigma_mp = 0;%=0 -> no mp when modeling MP fading as zero-mean, lognormal

cp = ChannelParams(q_b, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);


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
% Setup the Channel Analyzer With Observations, View Posterior Probability of Connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%minimum power required for connectivity
gamma_TH = -42;
obs_pos = [10,15; 20,24; 1,1; 10,20; 10, 5; 20, 15]*res;
obs_vals = [-38; -38; -45; -36; -46; -39];
cawo = CAWithObservations(cp, cc, obs_pos, obs_vals, gamma_TH);
figure(2)
clf
cawo.plotMeans()

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path Plan to Connected Spots Using Expected distance cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup obstacles
source = [13, 25];

maze = [RectObs([-1, -1.5, 3, -5.5]), RectObs([0, -0.5, 6, -2]), ...
            RectObs([1, 0.5, 3, -5.5]), RectObs([2, 1.5, 6, -2])];
circle_grid = [CircObs(1.5, [5,5]), CircObs(1.5, [5, 10]),...
    CircObs(1.5, [5,15]), CircObs(1.5, [5, 20]),...
    CircObs(1.5, [10,7.5]), CircObs(1.5, [10, 12.5]),...
    CircObs(1.5, [10,17.5]), CircObs(1.5, [10, 22.5]),...
    CircObs(1.5, [15, 5]), CircObs(1.5, [15, 10]),...
    CircObs(1.5, [15,15]), CircObs(1.5, [15, 20])];
simple = [CircObs(2,[10,10]), CircObs(2,[20,5]), CircObs(2,[5,20])]; 

goal = [13,0];

obs_mod = ObstacleMod(simple);
cost_fxn = @(n1, n2, path, mode) LIPNoConn(n1, n2, path, cawo);
%receiver_noise = -100; R = 8; BER = 0.000001; K = -1.5/log(5*BER);
%cost_fxn = @(n1, n2, path, mode) LICommEnergy(n1, n2, path, cawo, receiver_noise, R, K);
path_res = res;
problem_instance = PathPlanningProblem(region, path_res, source, goal, obs_mod, cost_fxn);

%%
hold on
problem_instance.plotProb(1)
%%

%solve using RRT

%stop on first solution, unless a solution hasn't
%been found
solution_not_required = 0;
max_iterations = Inf;
%max run time in seconds
max_run_time = 60;
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);

%1 - uniform random (RRT), 2 - deterministic sequence (RDT) 
%3 - informed set sampling, 4 - continuous uniform random
type = 1;
dest_freq = 20;
sequence = [];%sampler.sequence;%use previous run's sequence
sampler = Sampler(type, dest_freq, sequence);


do_rewire = 1;
steer_rad = sqrt(2);

rrt_solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);

rrt_solver.solve(problem_instance);
bsf = rrt_solver.getBSF();
rrt_path = bsf.pathToRoot(0);
%%
hold on
plot(rrt_path(:,1), rrt_path(:,2), '+-');