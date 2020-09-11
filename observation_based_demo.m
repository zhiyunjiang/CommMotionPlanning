%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Given channel params and a handful of readings, find
% "communication rich" path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Connectivity Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position of the base station (remote station or transmitter)
% Position of the base station (remote station or transmitter)
%q_b = [3 5.5];
q_b = [2 7];
% Path loss parameters
% The resulting path loss component in dB is 
% gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
% K_PL is the path loss constant in dB and n_PL is the path loss exponent.
n_PL = 3.86;  
K_PL = -41.34;

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
K_ric = 14.56;     

% corr_mp = 1 -> correlated multipath 
% corr_mp = 0 -> uncorrelated multipath 
corr_mp = 0;

%
sigma_mp = 1;%=0 -> no mp when modeling MP fading as zero-mean, lognormal

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
res = 2;
cc = CommChannel(cp, N_sin, region, res);
%simualte the full channel so that we can accurately sample
USE_LL = 2;
cc.simulateSH(); cc.simulateMP(USE_LL);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the Channel Analyzer With Observations, View Posterior Probability of Connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%minimum power required for connectivity
n_samples = 20;
[obs_pos, obs_vals] = cc.randSample(n_samples);
%%
gamma_TH = -105;
cawo = CAWithObservations(cp, cc, obs_pos, obs_vals, gamma_TH);
%%
figure(1)
clf
% cawo.plotMeans2D();
% figure(2)
% cawo.plotPosteriors2D();
% figure(3)
%receiver_noise = 1e-10; R = 8; BER = 1e-6; K = -1.5/log(5*BER);
%cawo.plotRequiredTXPower2D(receiver_noise, BER, R)
% figure(4)
p_min = 0.8;
cawo.plotConnected2D(p_min);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path Plan to Connected Spots Using Expected distance cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup obstacles
source = [45, 49];

maze = [RectObs([-1, -1.5, 3, -5.5]), RectObs([0, -0.5, 6, -2]), ...
            RectObs([1, 0.5, 3, -5.5]), RectObs([2, 1.5, 6, -2])];
circle_grid = [CircObs(1.5, [5,5]), CircObs(1.5, [5, 10]),...
    CircObs(1.5, [5,15]), CircObs(1.5, [5, 20]),...
    CircObs(1.5, [10,7.5]), CircObs(1.5, [10, 12.5]),...
    CircObs(1.5, [10,17.5]), CircObs(1.5, [10, 22.5]),...
    CircObs(1.5, [15, 5]), CircObs(1.5, [15, 10]),...
    CircObs(1.5, [15,15]), CircObs(1.5, [15, 20])];
simple = [CircObs(3,[41,15]), CircObs(3,[37,30]), CircObs(3,[44,43])]; 
no_obs = [];

goal = [45, 2];

obs_mod = ObstacleMod(simple);
thetas = [0,0.2, 0.4, 0.6, 0.8, 1]; 
for i = 1:length(thetas)
    theta = thetas(i);
cost_fxn = @(n1, n2, path, mode) MinPNoConn(n1, n2, path, cawo, p_min, theta);
%receiver noise in mW
%cost_fxn = @(n1, n2, path, mode) LIExpectedTotalEnergy(n1, n2, path, cawo, receiver_noise, R, K);
%path_res = res;
problem_instance = PathPlanningProblem(region, res, source, goal, obs_mod, cost_fxn);
if i == 1
hold on
problem_instance.plotProb(1)
end
%%%

%solve using RRT
%stop on first solution, unless a solution hasn't
%been found
solution_not_required = 0;
max_iterations = Inf;
%max run time in seconds
mins = 10;
max_run_time = mins*60;
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);

%1 - uniform random (RRT), 2 - deterministic sequence (RDT) 
%3 - informed set sampling, 4 - continuous uniform random
type = 1;
dest_freq = 100;
sequence = [];%sampler.sequence;%use previous run's sequence
sampler = Sampler(type, dest_freq, sequence);


do_rewire = 1;
steer_rad = 5;

rrt_solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);

rrt_solver.solve(problem_instance);
bsf = rrt_solver.getBSF();
rrt_path = bsf.pathToRoot(0);
%%%
hold on
plot_handles(i) = plot(rrt_path(:,1), rrt_path(:,2), 'LineWidth', 2, ...
            'DisplayName', strcat('RRT* Path, \theta = ', sprintf('%0.2f', theta)));

end
%%
legend(plot_handles);
%%
%get series of BSF and plot over time
series_data = rrt_solver.getBSFTimeSeries();
figure(2)
clf
plot(series_data(:,1)/60, series_data(:,2))
title('Minimum Energy Path');
xlabel('Time (min)')
xlim([0,mins])
ylabel('Energy (J)')