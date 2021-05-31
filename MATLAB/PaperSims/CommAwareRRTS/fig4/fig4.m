%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4
% For pre-defined workspace, obstacles, channel params, etc... runs
% multiple sims to find average costs. This is for variable transmission
% power, all tasks (uploading, broadcasting, and relaying).
% 
% Fig 4(left) - Upload
% Fig 4(middle) - Broadcast
% Fig 4(right) - Relay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Objective
objective = 2;
%Number of Base Stations
num_bs = 2;
%fetch predifned params
[n_PL, K_PL, sh_decorr, sigma_SH, PSD_at_f_c, mp_decorr, lambda, K_ric,...
    corr_mp, sigma_mp, region, N_sin, res, n_samples, n_sims, gamma_TH] = getConstParams(); 

%Setup Comm Channel Params objects
q_b = [2 2];
cp = ChannelParams(q_b, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);

q_poi = [43, 48];
cp_poi = ChannelParams(q_poi, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);

%Multipath model: 1 - rician, 2 - log normal
MP_Mdl = 2;%use log normal model for simulating multipath

p_th = 0.7;%doesn't actually matter here

%QoS related params
receiver_noise = 1e-10;  BER = 1e-6;

%k1 = 7.4; k2 = 0.29; v_const=1;%m/s
k1 = 7.4; k2 = 0.29; v_const=1;%m/s
mp = MotionParams(k1 , k2, v_const);
%motion weight
%delta = 1/100;works well, but is it reasonable?
delta = 1/10;

%Workspace params
obs_mod = ObstacleFactory.custom();
source = [1, 30]; goal = [49, 25];

% Setup arrays to store costs and related statistics
rrt_MIN_true_cost = zeros([n_sims, 1]); sp_MIN_true_cost = zeros([n_sims,1]);
rrt_MIN_length = zeros([n_sims, 1]); rrt_MIN_avg_com_power = zeros([n_sims, 1]);
sp_MIN_avg_com_power = zeros([n_sims, 1]);

rrt_MAX_true_cost = zeros([n_sims, 1]); sp_MAX_true_cost = zeros([n_sims,1]);
rrt_MAX_length = zeros([n_sims, 1]); rrt_MAX_avg_com_power = zeros([n_sims, 1]);
sp_MAX_avg_com_power = zeros([n_sims, 1]);

rrt_SUM_true_cost = zeros([n_sims, 1]); sp_SUM_true_cost = zeros([n_sims,1]);
rrt_SUM_length = zeros([n_sims, 1]); rrt_SUM_avg_com_power = zeros([n_sims, 1]);
sp_SUM_avg_com_power = zeros([n_sims, 1]);

for i = 1:n_sims
%find the straightline path - only need to do once per obstacle
if i == 1
    problem_instance_SP = PathPlanningProblem(region, res, source, goal, obs_mod, getCostFxnPredicted(0));
    rrt_solver_SP = getRRTSolver();
    rrt_path_SP_f1 = runOneSim(rrt_solver_SP, problem_instance_SP);
    sp_f1_length = GridDist(rrt_path_SP_f1);
end
       
%Simulate a new channel
cc1_fig_4 = CommChannel(cp, N_sin, region, res);
cc1_fig_4.simulateSH(); cc1_fig_4.simulateMP(MP_Mdl);

cc2_fig_4 = CommChannel(cp_poi, N_sin, region, res);
cc2_fig_4.simulateSH(); cc2_fig_4.simulateMP(MP_Mdl);

% Sample the channel
[obs_pos, obs_vals] = cc1_fig_4.randSample(n_samples);
[obs_pos_poi, obs_vals_poi] = cc2_fig_4.randSample(n_samples);

% Predict the channel
pc1 = PredictedChannel(cp, cc1_fig_4, obs_pos, obs_vals);
pc2 = PredictedChannel(cp_poi, cc2_fig_4, obs_pos_poi, obs_vals_poi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4 (left) - min energy - Upload
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scenario = 1;

%Energy Related Params
R1 = 8;
qos1 = QoSParams(BER, R1, receiver_noise);


plotPredictedField(objective, scenario, num_bs, pc1, pc2, qos1, p_th, gamma_TH);
caxis([-40, 10])
cost_fxn = getCostFxnPredicted(objective, scenario, num_bs, pc1, pc2, mp, qos1, delta, p_th, gamma_TH);
problem_instance_MIN = PathPlanningProblem(region, res, source, goal, obs_mod, cost_fxn);
rs_handle = plotComAwareRRT(problem_instance_MIN, q_b, q_poi);

rrt_solver_MIN = getRRTSolver();
rrt_path_MIN = runOneSim(rrt_solver_MIN, problem_instance_MIN);
true_costfxn_MIN = getCostFxnTrue(objective, scenario, num_bs, cc1_fig_4, cc2_fig_4, mp, qos1, delta, gamma_TH);

path_handles = plotPaths(rrt_path_MIN, 'RRT* Path (Upload)');
legend([path_handles, rs_handle]);
caxis([-20,0])
saveas(gcf, sprintf('fig4a_run_%d', i));
close(gcf);
fprintf('Iteration %d: created fig4a\n',i);

rrt_MIN_true_cost(i) = true_costfxn_MIN([],[],rrt_path_MIN,1);
sp_MIN_true_cost(i) = true_costfxn_MIN([],[],rrt_path_SP_f1,1);

rrt_MIN_length(i) = GridDist(rrt_path_MIN);
rrt_MIN_avg_com_power(i) = (rrt_MIN_true_cost(i) - delta*mp.motionEnergy(rrt_MIN_length(i)))/(rrt_MIN_length(i));
sp_MIN_avg_com_power(i) = (sp_MIN_true_cost(i) - delta*mp.motionEnergy(sp_f1_length))/(sp_f1_length);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4 (mid) - min energy - Broadcast
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scenario = 2;

%Energy Related Params
R2 = 6;
qos2 = QoSParams(BER, R2, receiver_noise);

plotPredictedField(objective, scenario, num_bs, pc1, pc2, qos2, p_th, gamma_TH);

cost_fxn = getCostFxnPredicted(objective, scenario, num_bs, pc1, pc2, mp, qos2, delta, p_th, gamma_TH);
problem_instance_MAX = PathPlanningProblem(region, res, source, goal, obs_mod, cost_fxn);
caxis([0,20])

rs_handle = plotComAwareRRT(problem_instance_MAX, q_b, q_poi);
rrt_solver_MAX = getRRTSolver();
rrt_path_MAX = runOneSim(rrt_solver_MAX, problem_instance_MAX);
true_costfxn_MAX = getCostFxnTrue(objective, scenario, num_bs, cc1_fig_4, cc2_fig_4, mp, qos2, delta, gamma_TH);

path_handles = plotPaths(rrt_path_MAX, 'RRT* Path (Broadcast)');
legend([path_handles, rs_handle]);
caxis([-20,0])
saveas(gcf, sprintf('fig4b_run_%d', i));
close(gcf);
fprintf('Iteration %d: created fig4b\n',i);

rrt_MAX_true_cost(i) = true_costfxn_MAX([],[],rrt_path_MAX,1);
sp_MAX_true_cost(i) = true_costfxn_MAX([],[],rrt_path_SP_f1,1);

rrt_MAX_length(i) = GridDist(rrt_path_MAX);
rrt_MAX_avg_com_power(i) = (rrt_MAX_true_cost(i) - delta*mp.motionEnergy(rrt_MAX_length(i)))/(rrt_MAX_length(i));
sp_MAX_avg_com_power(i) = (sp_MAX_true_cost(i) - delta*mp.motionEnergy(sp_f1_length))/(sp_f1_length);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 4 (right) - min energy - Relay
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scenario = 3;

plotPredictedField(objective, scenario, num_bs, pc1, pc2, qos2, p_th, gamma_TH);
caxis([0,20])

cost_fxn = getCostFxnPredicted(objective, scenario, num_bs, pc1, pc2, mp, qos2, delta, p_th, gamma_TH);
problem_instance_SUM = PathPlanningProblem(region, res, source, goal, obs_mod, cost_fxn);

rs_handle = plotComAwareRRT(problem_instance_SUM, q_b, q_poi);
rrt_solver_SUM = getRRTSolver();
rrt_path_SUM = runOneSim(rrt_solver_SUM, problem_instance_SUM);
true_costfxn_SUM = getCostFxnTrue(objective, scenario, num_bs, cc1_fig_4, cc2_fig_4, mp, qos2, delta, gamma_TH);

path_handles = plotPaths(rrt_path_SUM, 'RRT* Path (Relay)');
legend([path_handles, rs_handle])
caxis([-20,0])
saveas(gcf, sprintf('fig4c_run_%d', i));
close(gcf);
fprintf('Iteration %d: created fig4c\n',i);

rrt_SUM_true_cost(i) = true_costfxn_SUM([],[],rrt_path_SUM,1);
sp_SUM_true_cost(i) = true_costfxn_SUM([],[],rrt_path_SP_f1,1);

rrt_SUM_length(i) = GridDist(rrt_path_SUM);
rrt_SUM_avg_com_power(i) = (rrt_SUM_true_cost(i) - delta*mp.motionEnergy(rrt_SUM_length(i)))/(rrt_SUM_length(i));
sp_SUM_avg_com_power(i) = (sp_SUM_true_cost(i) - delta*mp.motionEnergy(sp_f1_length))/(sp_f1_length);

end
save(strcat('fig4_ws_',date));