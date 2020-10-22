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
q_b = [7 14];
% Path loss parameters
% The resulting path loss component in dB is 
% gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
% K_PL is the path loss constant in dB and n_PL is the path loss exponent.
%n_PL = 3.86;
n_PL = 4;
% K_PL = -41.34;
K_PL = -65;
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

q_poi = [45, 38];
cp_poi = ChannelParams(q_poi, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);
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

cc_poi = CommChannel(cp_poi, N_sin, region, res);
cc_poi.simulateSH(); cc_poi.simulateMP(USE_LL);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the Channel Analyzer With Observations, View Posterior Probability of Connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%minimum power required for connectivity
n_samples = 500;
[obs_pos, obs_vals] = cc.randSample(n_samples);
[obs_pos_poi, obs_vals_poi] = cc_poi.randSample(n_samples);
%%
%actually calculate this
gamma_TH = -95;
cawo = CAWithObservations(cp, cc, obs_pos, obs_vals, gamma_TH);
cawo_poi = CAWithObservations(cp_poi, cc_poi, obs_pos_poi, obs_vals_poi, gamma_TH);
%%
%Objective
%1 - min disconnected distance
%2 - min energy
objective = 1;
% Scenario
%1 - min/OR
%2 - max/AND
%3 - sum
scenario = 1;
%Number of Base Stations
num_bs = 2;

%disconnected distance related parameters
eps = 0.01;
p_th = 0.7;

%min-energy related parameters
receiver_noise = 1e-10; R = 6; BER = 1e-6;
qos = QoSParams(BER, R, receiver_noise);

k1 = 7; k2 = 0.3; v_const=1;%m/s
mp = MotionParams(k1 , k2, v_const);
    
%single base station
connection_field = cawo.getConnectionField(p_th);
%%
%cc.plotField(connection_field, 'Connected Regions')
%cawo.plotPosteriors2D();
cawo.plotMeans2D()
%%
figure();
req_power = cawo.getEReqTXPowerW(qos);
cc.plotField(10*log10(req_power'), 'Required TX Power for Communication', 'Power (dB)');
%%
%multi-base station, min-disconnected distance
% conn_field = cawo.getConnectionField(p_th) + 2*cawo_poi.getConnectionField(p_th);
cc.plotDoubleConnectField(conn_field, 'Connected Regions')
%multi - base station, energy based
%%
if scenario == 1
    req_power = min(cawo.getEReqTXPowerW(qos), cawo_poi.getEReqTXPowerW(qos));
elseif scenario == 2
    req_power = max(cawo.getEReqTXPowerW(qos), cawo_poi.getEReqTXPowerW(qos));
elseif scenario == 3
    req_power = cawo.getEReqTXPowerW(qos) + cawo_poi.getEReqTXPowerW(qos);
end

cc.plotField(10*log10(req_power'), 'Required TX Power for Communication', 'Required Transmit Power (dB)');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path Plan to Connected Spots Using Expected distance cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup obstacles
source = [2, 45];
goal = [49, 15];
%%
%obs_mod = ObstacleFactory.circleGrid(2.5, 2, [-5,5] ,source, goal);
%obs_mod = ObstacleFactory.custom();

%%
cost_fxn = getCostFxn(objective, scenario, num_bs, cawo, cawo_poi, mp, qos, eps, p_th);
problem_instance = PathPlanningProblem(region, res, source, goal, obs_mod, cost_fxn);
problem_instance_sp = PathPlanningProblem(region, res, source, goal, obs_mod, getCostFxn(0));
%%
plotComAwareRRT(problem_instance, q_b, q_poi)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Setup RRT* Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%stopping criteria
solution_not_required = 0; max_iterations = Inf; mins = 3; max_run_time = mins*60;
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);
%sampler setup
%1 - uniform random (RRT), 2 - deterministic sequence (RDT) 
%3 - informed set sampling, 4 - continuous uniform random
type = 1; dest_freq = 100; sequence = [];%sampler.sequence;%use previous run's sequence
sampler = Sampler(type, dest_freq, sequence);

%0 - RDT, 1 - RDT*
do_rewire = 1;steer_rad = 20;
rrt_solver = RDTSolver(sampler, stop_criteria, do_rewire, steer_rad);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run several times to find average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_sims = 10;
[avg_final_cost, avg_first_time, avg_first_cost] = runMultiSim(n_sims, rrt_solver, problem_instance,...
                                                objective, eps, res);
                                            
[sp_avg_final_cost, sp_avg_first_time, sp_avg_first_cost] = runMultiSim(n_sims, rrt_solver, problem_instance_sp,...
                                                0, eps, res);
diffs = [sp_avg_final_cost - avg_final_cost, sp_avg_first_time - avg_first_time, sp_avg_first_cost-avg_first_cost];
%%
hold on
% plot_handles(i) = plot(rrt_path(:,1)/res, rrt_path(:,2)/res, 'LineWidth', 1, ...
%             'DisplayName', strcat('RRT* Path, R = ', sprintf('%0.2f', R)));
% legend(plot_handles);
plot_handle_1 = plot(rrt_path(:,1)/res, rrt_path(:,2)/res, 'k:',...
    'LineWidth', 2.5);
title('')
%%
legend([plot_handle_2, handle1]);
%%
%get series of BSF and plot over time
[costs, dists] = rrt_solver.getBSFTimeSeries();
disc_dist = (costs(:,2) - eps*dists(:,2));
figure();
%plot(costs(:,1), disc_dist/res)
plot(costs(:,1), costs(:,2)/1000);
%%
xlabel('Time (s)')
xlim([0,mins*60])
ylabel('Energy (kJ)')
%%

function rrt_path = runOneSim(rrt_solver, pppi)
    rrt_solver.solve(pppi);
    bsf = rrt_solver.getBSF();
    rrt_path = bsf.pathToRoot(0);
end

function cost =  calcPresentedCost(cost, dist, objective, eps, res)
    
    if objective == 1 %minimum disconnected distance, in meters
        cost = (cost - dist*eps)/res;
    elseif objective == 2
        cost = cost/1000;%convert joules to kJ, a much more tractable unit
    end
end

function [avg_final_cost, avg_first_time, avg_first_cost] = runMultiSim(n_sims, rrt_solver,...
                                                                pppi, objective, eps, res)
    avg_final_cost = 0; avg_first_time = 0; avg_first_cost = 0;
    for i = 1:n_sims
        runOneSim(rrt_solver, pppi);
        [costs, dists] = rrt_solver.getBSFTimeSeries();
        avg_final_cost = avg_final_cost + calcPresentedCost(costs(end,2), dists(end,2), objective, eps, res);
        avg_first_time = costs(1,1);
        avg_first_cost = costs(1,2);
    end
    avg_final_cost = avg_final_cost/n_sims;
    avg_first_time = avg_first_time/n_sims;
    avg_first_cost = avg_first_cost/n_sims;
end

function cost_fxn = getCostFxn(objective, scenario, num_bs, cawo, cawo_poi, mp, qos, eps, p_th)

    if objective == 0%min distance cost function
        cost_fxn = @(n1, n2, path, mode) GridDist(path);
    elseif objective == 1%min disconnected distance
        if num_bs == 1
            cost_fxn = @(n1, n2, path, mode) MinPNoConn(n1, n2, path, cawo, p_th, eps);
        elseif num_bs == 2
            cost_fxn = @(n1, n2, path, mode) MinPNoConnWithPOI(path, cawo, cawo_poi, p_th, eps, scenario);
        else
            error('Current implementation of disconnected distance cost functions only handles up to 2 remote stations');
        end
    elseif objective == 2
        if num_bs == 1
            cost_fxn = @(n1, n2, path, mode) LIExpectedTotalEnergy(path, cawo, qos, mp);
        elseif num_bs == 2
            cost_fxn = @(n1, n2, path, mode) LIExpTotalEnergyWithPOI(path, cawo, cawo_poi, qos, mp, scenario);
        else
            error('Current implementation of energy cost functions only handles up to 2 remote stations');
        end
    else
        error('Invalid objective id');
    end
end