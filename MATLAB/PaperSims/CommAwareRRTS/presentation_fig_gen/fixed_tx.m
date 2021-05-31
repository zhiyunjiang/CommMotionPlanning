%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slide 1
% For pre-defined workspace, obstacles, channel params, etc... runs
% multiple sims to find average costs. This is for fixed transmission
% power, all tasks (uploading, broadcasting, and relaying).
% 
% Fig 2 - min disconnected distance (OR & AND, single plot)
% Fig 3 - Convergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Objective
%1 - min disconnected distance %2 - min energy
objective = 1;
%Number of Base Stations
num_bs = 2;

%Setup Comm Channel Params
[n_PL, K_PL, sh_decorr, sigma_SH, PSD_at_f_c, mp_decorr, lambda, K_ric,...
    corr_mp, sigma_mp, ~, N_sin, ~, n_samples, n_sims, ~] = getConstParams(); 

region = [25,0, 25, 0];
res = 10;
gamma_TH = -90;

q_b = [1,1];
cp = ChannelParams(q_b, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);
 
q_poi = [24, 24];
cp_poi = ChannelParams(q_poi, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);

%Multipath model: 1 - rician, 2 - log normal
MP_Mdl = 2;%use log normal model for simulating multipath

% Energy Related Params
%QoS Params (fig 4)
receiver_noise = 1e-10; R = 2; BER = 1e-6;
qos = QoSParams(BER, R, receiver_noise);
% Calculate Required TX power based on QoS requirements
p_th = 0.85;
fixed_tx_power_W = qos.reqTXPower(gamma_TH);

%Motion Params
k1 = 3.5; k2 = 0.3; v_const=1;%m/s
mp = MotionParams(k1 , k2, v_const);


%Workspace/path params
source = [2,12]; goal = [23, 18];
%use same obstacles for accurate comparison
obs_mod = ObstacleFactory.slide1();


%Formatting params
line_width = 2;

%Problem Setup

% generate the full channel to sample
cc1 = CommChannel(cp, N_sin, region, res);
cc1.simulateSH(); cc1.simulateMP(MP_Mdl);

cc2 = CommChannel(cp_poi, N_sin, region, res);
cc2.simulateSH(); cc2.simulateMP(MP_Mdl);


%Sample
[obs_pos, obs_vals] = cc1.randSample(n_samples);
[obs_pos_poi, obs_vals_poi] = cc2.randSample(n_samples);


%Create Predicted Channel
pc1 = PredictedChannel(cp, cc1, obs_pos, obs_vals);
pc2 = PredictedChannel(cp_poi, cc2, obs_pos_poi, obs_vals_poi);


deltas = [0.01, 0.5, 1, 2];%, 100];
%%

%UPLOADING
rrt_OR_true_cost = zeros([length(deltas),1]);
rrt_OR_length = zeros([length(deltas),1]);
rrt_OR_disconnected_length = zeros([length(deltas),1]);
cost_fxn_dl = getCostFxnPredicted(objective, scenario, num_bs, pc1, pc2, mp, qos, 0, p_th, gamma_TH);
scenario = 1;
figure(1)
plotPredictedField(objective, scenario, num_bs, pc1, pc2, qos, p_th, gamma_TH);


for i=1:length(deltas)
   delta = deltas(i);
    
    cost_fxn = getCostFxnPredicted(objective, scenario, num_bs, pc1, pc2, mp, qos, delta, p_th, gamma_TH);
    problem_instance_OR = PathPlanningProblem(region, res, source, goal, obs_mod, cost_fxn);
    if i == 1
        rs_handle = plotComAwareRRT(problem_instance_OR, q_b, q_poi);
    end
    
    rrt_solver_OR = getRRTSolver();
    rrt_path_OR = runOneSim(rrt_solver_OR, problem_instance_OR);

    hold on
    path_handles = plotPaths(rrt_path_OR, 'RRT* Path (Upload)');

    fprintf('Completed for eps = %f\n',delta);
    
    rrt_OR_true_cost(i) = cost_fxn([],[],rrt_path_OR,1);

    rrt_OR_length(i) = GridDist(rrt_path_OR);
    
    rrt_OR_disconnected_length(i) = cost_fxn_dl([],[],rrt_path_OR,1);
end
%%
%BROADCASTING/RELAY
scenario = 2;%broadcasting - AND
figure(2)

p_th = 0.9;
gamma_TH = -95;
source = [2,16]; goal = [23, 2];

rrt_AND_true_cost = zeros([length(deltas),1]);
rrt_AND_length = zeros([length(deltas),1]);
rrt_AND_disconnected_length = zeros([length(deltas),1]);
cost_fxn_dl = getCostFxnPredicted(objective, scenario, num_bs, pc1, pc2, mp, qos, 0, p_th, gamma_TH);


plotPredictedField(objective, scenario, num_bs, pc1, pc2, qos, p_th, gamma_TH);

for i=1:length(deltas)
    delta = deltas(i);
    cost_fxn = getCostFxnPredicted(objective, scenario, num_bs, pc1, pc2, mp, qos, delta, p_th, gamma_TH);
    problem_instance_AND = PathPlanningProblem(region, res, source, goal, obs_mod, cost_fxn);
     if i == 1
        rs_handle = plotComAwareRRT(problem_instance_AND, q_b, q_poi);
    end
    rrt_solver_AND = getRRTSolver();
    rrt_path_AND = runOneSim(rrt_solver_AND, problem_instance_AND);
    hold on
    path_handles = plotPaths(rrt_path_AND, sprintf('\\delta = %f', delta) );
    
    rrt_AND_true_cost(i) = cost_fxn([],[],rrt_path_AND,1);

    rrt_AND_length(i) = GridDist(rrt_path_AND);
    
    rrt_AND_disconnected_length(i) = cost_fxn_dl([],[],rrt_path_AND,1);
    
end
