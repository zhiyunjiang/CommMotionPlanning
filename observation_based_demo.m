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
%n_PL = 3.86;
n_PL = 1.8;
K_PL = -41.34;

%Shadowing Parameters
%Shadowing decorrelation distance
sh_decorr = 3.09;  
%sh_decorr = 5;  
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

% q_poi = [45, 49];
% cp_poi = ChannelParams(q_poi, n_PL, K_PL - 5, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);
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
%simualte the full channel so that we can accurately sample
USE_LL = 2;
cc.simulateSH(); cc.simulateMP(USE_LL);

% cc_poi = CommChannel(cp_poi, N_sin, region, res);
% cc_poi.simulateSH(); cc_poi.simulateMP(USE_LL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the Channel Analyzer With Observations, View Posterior Probability of Connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%minimum power required for connectivity
n_samples = 500;
[obs_pos, obs_vals] = cc.randSample(n_samples);
% [obs_pos_poi, obs_vals_poi] = cc_poi.randSample(n_samples);
%%
gamma_TH = -95;
cawo = CAWithObservations(cp, cc, obs_pos, obs_vals, gamma_TH);
% cawo_poi = CAWithObservations(cp_poi, cc_poi, obs_pos_poi, obs_vals_poi, gamma_TH);
%%
%1 - min/OR
%2 - max/AND
%3 - sum
scenario = 3;
receiver_noise = 1e-10; R = 6; BER = 1e-6;
qos = QoSParams(BER, R, receiver_noise);
p_th = 0.65;

%single base station
connection_field = cawo.getConnectionField(p_th);
%%
cc.plotField(connection_field, 'Connected Regions')
%cawo.plotPosteriors2D();

%multi - base station, energy based
% if scenario == 1
%     req_power = min(cc.getReqTXPowermW(qos), cc_poi.getReqTXPowermW(qos));
% elseif scenario == 2
%     req_power = max(cc.getReqTXPowermW(qos), cc_poi.getReqTXPowermW(qos));
% elseif scenario == 3
%     req_power = cc.getReqTXPowermW(qos) + cc_poi.getReqTXPowermW(qos);
% end
% cc.plotField(10*log10(req_power), 'Required TX Power for Communication', 'Power (dB)');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path Plan to Connected Spots Using Expected distance cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup obstacles
source = [49, 1];
goal = [20, 49];

obs_mod = ObstacleFactory.circleGrid(2.5, 2, [-5,5] ,source, goal);

%%
%set motion parameters
k1 = 5; k2 = 0.3; v_const=1;%m/s
mp = MotionParams(k1 , k2, v_const);

theta = 0.01;
cost_fxn = @(n1, n2, path, mode) MinPNoConn(n1, n2, path, cawo, p_th, theta);
% cost_fxn = @(n1, n2, path, mode) MinPNoConnPOI(n1, n2, path, cawo, cawo_poi, p_min, theta, scenario);
% R = 6;
%cost_fxn = @(n1, n2, path, mode) LIExpTotalEnergyWithPOI(n1, n2, path, cawo, receiver_noise, R, K, cawo_poi);
problem_instance = PathPlanningProblem(region, res, source, goal, obs_mod, cost_fxn);
hold on
problem_instance.plotProb(0)
scatter(q_b(1), q_b(2), 'v', 'filled');
%%
%solve using RRT
%stop on first solution, unless a solution hasn't
%been found
solution_not_required = 0;
max_iterations = Inf;
%max run time in seconds
mins = 5;
max_run_time = mins*60;
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);

%1 - uniform random (RRT), 2 - deterministic sequence (RDT) 
%3 - informed set sampling, 4 - continuous uniform random
type = 1;
dest_freq = 100;
sequence = [];%sampler.sequence;%use previous run's sequence
sampler = Sampler(type, dest_freq, sequence);


do_rewire = 1;
steer_rad = 20;

rrt_solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);

rrt_solver.solve(problem_instance);
bsf = rrt_solver.getBSF();
rrt_path = bsf.pathToRoot(0);
%%
hold on
% plot_handles(i) = plot(rrt_path(:,1)/res, rrt_path(:,2)/res, 'LineWidth', 1, ...
%             'DisplayName', strcat('RRT* Path, R = ', sprintf('%0.2f', R)));
% legend(plot_handles);
plot_handle = plot(rrt_path(:,1)/res, rrt_path(:,2)/res, 'r', 'LineWidth', 1);
legend(plot_handle, 'RRT* path');
%%
%get series of BSF and plot over time
series_data = rrt_solver.getBSFTimeSeries();
figure(2)
clf
plot(series_data(:,1)/60, series_data(:,2)/res)
title('RRT* Best Path So Far');
xlabel('Time (min)')
xlim([0,mins])
ylabel('Distance (m)')