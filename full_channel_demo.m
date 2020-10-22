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
q_b = [25 1];
% Path loss parameters
% The resulting path loss component in dB is 
% gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
% K_PL is the path loss constant in dB and n_PL is the path loss exponent.\
% n_PL = 3.86;
n_PL = 4;
% K_PL = -41.34;
K_PL = -65;


%Shadowing Parameters
%Shadowing decorrelation distance
sh_decorr = 3.09;
%shadowing power
% alpha = 10.24;
alpha = 8.5;
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

%make the transmitter on the mobile robot weaker than that of the base station
% q_poi = [49,38];
% cp_poi = ChannelParams(q_poi, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);
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

% cc_poi = CommChannel(cp_poi, N_sin, region, res);
% cc_poi.simulateSH(); cc_poi.simulateMP();
%%
%1 - min/OR
%2 - max/AND
%3 - sum
scenario = 3;
receiver_noise = 1e-10; R = 6; BER = 1e-6;
qos = QoSParams(BER, R, receiver_noise);

cc.plotRequiredTXPower(qos)
% figure();
%gamma_th = -126;   
% conn_field = cc_poi.getConnectionField(gamma_th) + 2*cc.getConnectionField(gamma_th);
% cc.plotField(conn_field,'Regions of Connectivity');
% cc.plotDoubleConnectField(conn_field, 'Connectivity Regions')
% figure(2)

% if scenario == 1
%     req_power = min(cc.getReqTXPowermW(qos), cc_poi.getReqTXPowermW(qos));
% elseif scenario == 2
%     req_power = max(cc.getReqTXPowermW(qos), cc_poi.getReqTXPowermW(qos));
% elseif scenario == 3
%     req_power = cc.getReqTXPowermW(qos) + cc_poi.getReqTXPowermW(qos);
% end
% cc.plotField(10*log10(req_power), 'Required TX Power for Communication', 'Power (dB)');
%%
%setup obstacles
source = [1, 25];
goal = [49, 40];

%generat obstacles
obs_mod = ObstacleFactory.circleGrid(2.5, 2, [-5,5] ,source, goal);
%%
%set motion parameters
k1 = 5; k2 = 0.3; v_const=1;%m/s
mp = MotionParams(k1 , k2, v_const);
%%
% theta = 0.01;
% cost_fxn = @(n1, n2, path, mode) MinNoConn(n1, n2, path, cc, gamma_th, theta);
% use_or = 0;
% cost_fxn = @(n1, n2, path, mode) MinNoConnWithPOI(n1, n2, path, cc, cc_poi, gamma_th, theta, use_or);
cost_fxn = @(n1, n2, path, mode) LITotalEnergy(path, cc, qos, mp);

% cost_fxn = @(n1, n2, path, mode) LITotalEnergyWithPOI(path, cc, cc_poi, qos, mp, scenario);
path_res = res;
problem_instance = PathPlanningProblem(region, path_res, source, goal, obs_mod, cost_fxn);

hold on
problem_instance.plotProb(0)
hold on
scatter(q_b(1), q_b(2), 'v', 'filled');
% scatter(q_poi(1), q_poi(2), '^', 'filled');
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

hold on
plot_handle = plot(rrt_path(:,1)/res, rrt_path(:,2)/res,  'k-.', 'LineWidth', 1, 'DisplayName', 'RRT*, Max Power');
%legend([OR_plothandle, plot_handle]);
%%
%compare to discrete planners (A*, Dijkstra's)
%how long does it take without any Heuristic (Dijkstra's)
% astar_solver = AStarSolver();
% tic
% astar_solver.solve(problem_instance);
% toc
% bst = astar_solver.BST;
% astar_path = bst.pathToRoot(0);
% 
% hold on
% plot(astar_path(:,1), astar_path(:,2),'g.', 'LineWidth', 2);
%%
%get series of BSF and plot over time
series_data = rrt_solver.getBSFTimeSeries();
figure(2)
clf
plot(series_data(:,1)/60, series_data(:,2))
title('RRT* Minimum Energy Path');
% title('RRT* Minimum Disconnected Distance Path');
xlabel('Time (min)')
xlim([0,mins])
ylabel('Energy (J)')
% ylabel('Disconnected Distance')
