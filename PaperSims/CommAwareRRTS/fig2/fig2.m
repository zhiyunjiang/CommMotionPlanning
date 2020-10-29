%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fig 2
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
    corr_mp, sigma_mp, region, N_sin, res, n_samples, n_sims, gamma_TH] = getConstParams(); 

q_b = [2,7];
cp = ChannelParams(q_b, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);
 
q_poi = [45, 49];
cp_poi = ChannelParams(q_poi, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);

%Multipath model: 1 - rician, 2 - log normal
MP_Mdl = 2;%use log normal model for simulating multipath


% Energy Related Params
%QoS Params (fig 4)
receiver_noise = 1e-10; R = 2; BER = 1e-6;
qos = QoSParams(BER, R, receiver_noise);
% Calculate Required TX power based on QoS requirements
p_th = 0.70;
fixed_tx_power_W = qos.reqTXPower(gamma_TH);

%Motion Params
k1 = 3.5; k2 = 0.3; v_const=1;%m/s
mp = MotionParams(k1 , k2, v_const);
%motion weight
eps = 0.01;


%Workspace/path params
source = [1, 35]; goal = [40, 49];
%use same obstacles for accurate comparison
obs_mod = ObstacleFactory.fig2();


%Formatting params
line_width = 2;

%setup arrays to store costs and related values
rrt_OR_true_cost = zeros([n_sims,1]); sp_OR_true_cost = zeros([n_sims,1]);
rrt_OR_length = zeros([n_sims,1]); rrt_OR_disc_dist = zeros([n_sims,1]);
sp_OR_disc_dist = zeros([n_sims,1]); rrt_OR_per_disc = zeros([n_sims,1]);
sp_OR_per_disc = zeros([n_sims,1]);

rrt_AND_true_cost = zeros([n_sims,1]); sp_AND_true_cost = zeros([n_sims,1]);
rrt_AND_length = zeros([n_sims,1]); rrt_AND_disc_dist = zeros([n_sims,1]);
sp_AND_disc_dist = zeros([n_sims,1]); rrt_AND_per_disc = zeros([n_sims,1]);
sp_AND_per_disc = zeros([n_sims,1]);

for i = 1:n_sims

% generate the full channel to sample
cc1 = CommChannel(cp, N_sin, region, res);
cc1.simulateSH(); cc1.simulateMP(MP_Mdl);

cc2 = CommChannel(cp_poi, N_sin, region, res);
cc2.simulateSH(); cc2.simulateMP(MP_Mdl);


%Sample
[obs_pos, obs_vals] = cc1.randSample(n_samples);
[obs_pos_poi, obs_vals_poi] = cc2.randSample(n_samples);


%Create Predicted Channel
cawo = CAWithObservations(cp, cc1, obs_pos, obs_vals);
cawo_poi = CAWithObservations(cp_poi, cc2, obs_pos_poi, obs_vals_poi);


%UPLOADING
scenario = 1;

plotPredictedField(objective, scenario, num_bs, cawo, cawo_poi, qos, p_th, gamma_TH);

cost_fxn = getCostFxnPredicted(objective, scenario, num_bs, cawo, cawo_poi, mp, qos, eps, p_th, gamma_TH);
problem_instance_OR = PathPlanningProblem(region, res, source, goal, obs_mod, cost_fxn);

rs_handle = plotComAwareRRT(problem_instance_OR, q_b, q_poi);
rrt_solver_OR = getRRTSolver();
rrt_path_OR = runOneSim(rrt_solver_OR, problem_instance_OR);
true_costfxn_OR = getCostFxnTrue(objective, scenario, num_bs, cc1, cc2, mp, qos, eps, gamma_TH);

%BROADCASTING/RELAY
scenario = 2;%broadcasting - AND
cost_fxn = getCostFxnPredicted(objective, scenario, num_bs, cawo, cawo_poi, mp, qos, eps, p_th, gamma_TH);
problem_instance_AND = PathPlanningProblem(region, res, source, goal, obs_mod, cost_fxn);
rrt_solver_AND = getRRTSolver();
rrt_path_AND = runOneSim(rrt_solver_AND, problem_instance_AND);
true_costfxn_AND = getCostFxnTrue(objective, scenario, num_bs, cc1, cc2, mp, qos, eps, gamma_TH);

path_handles = plotPaths(rrt_path_OR/res, 'RRT* Path (Upload)', rrt_path_AND/res, 'RRT* Path (Broadcast)');

legend([path_handles, rs_handle])
saveas(gcf, sprintf('fig2_run_%d', i));
close(gcf);
fprintf('Iteration %d: created fig2\n',i);
%should only ever need to run this once per obstacle/start&stop configuration
if i == 1
    problem_instance_SP = PathPlanningProblem(region, res, source, goal, obs_mod, getCostFxnPredicted(0));
    rrt_solver_SP = getRRTSolver();
    rrt_path_SP_f2 = runOneSim(rrt_solver_SP, problem_instance_SP);
    sp_f2_length = GridDist(rrt_path_SP_f2);
end

%Compare costs to baseline
% Upload Costs (OR)
rrt_OR_true_cost(i) = true_costfxn_OR([],[],rrt_path_OR,1);
sp_OR_true_cost(i) = true_costfxn_OR([],[],rrt_path_SP_f2,1);

rrt_OR_length(i) = GridDist(rrt_path_OR);
rrt_OR_disc_dist(i) = calcDiscDistance(rrt_OR_true_cost(i), rrt_OR_length(i), eps);
sp_OR_disc_dist(i) = calcDiscDistance(sp_OR_true_cost(i), sp_f2_length, eps);
rrt_OR_per_disc(i) = rrt_OR_disc_dist(i)/rrt_OR_length(i);
sp_OR_per_disc(i) = sp_OR_disc_dist(i)/sp_f2_length;

%Broadcast/Relay (AND)
rrt_AND_true_cost(i) = true_costfxn_AND([],[],rrt_path_AND,1);
sp_AND_true_cost(i) = true_costfxn_AND([],[],rrt_path_SP_f2,1);

rrt_AND_length(i) = GridDist(rrt_path_AND);
rrt_AND_disc_dist(i) = calcDiscDistance(rrt_AND_true_cost(i), rrt_AND_length(i), eps);
sp_AND_disc_dist(i) = calcDiscDistance(sp_AND_true_cost(i), sp_f2_length, eps);
rrt_AND_per_disc(i) = rrt_AND_disc_dist(i)/rrt_AND_length(i);
sp_AND_per_disc(i) = sp_AND_disc_dist(i)/sp_f2_length;


%Save off time series in plot
%get series of BSF and plot over time
[costs, dists] = rrt_solver_AND.getBSFTimeSeries();
disc_dist = (costs(:,2) - eps*dists(:,2));
plot(costs(:,1), disc_dist/res, 'LineWidth', line_width)
xlabel('Time (s)')
xlim([0,3*60])
ylabel('Cost (m)')
saveas(gcf, sprintf('fig3a_run_%d', i));
close(gcf);

fprintf('Iteration %d: created fig3a\n',i);
[costs, dists] = rrt_solver_OR.getBSFTimeSeries();
disc_dist = (costs(:,2) - eps*dists(:,2));
plot(costs(:,1), disc_dist/res, 'LineWidth', line_width)
xlabel('Time (s)')
xlim([0,3*60])
ylabel('Cost (m)')
saveas(gcf, sprintf('fig3b_run_%d', i));
close(gcf);
fprintf('Iteration %d: created fig3b\n',i);

end
%%
%compute averages and compare
fprintf('UPLOAD STATS\n')
rrt_OR_true_cost_mean = mean(rrt_OR_true_cost);
sp_OR_true_cost_mean = mean(sp_OR_true_cost);
OR_true_diff = sp_OR_true_cost_mean - rrt_OR_true_cost_mean;
fprintf('RRT True: %.2f, SL True: %.2f, Raw Improvement: %.2f, %% Improvement: %.2f\n',...
    rrt_OR_true_cost_mean, sp_OR_true_cost_mean, OR_true_diff , OR_true_diff/sp_OR_true_cost_mean);

rrt_OR_length_mean = mean(rrt_OR_length);
rrt_OR_disc_dist_mean = mean(rrt_OR_disc_dist);
sp_OR_disc_dist_mean = mean(sp_OR_disc_dist);
rrt_OR_per_disc_mean = mean(rrt_OR_per_disc);
sp_OR_per_disc_mean = mean(sp_OR_per_disc);



fprintf('\nBROADCAST/RELAY STATS\n')
rrt_AND_true_cost_mean = mean(rrt_AND_true_cost(i));
sp_AND_true_cost_mean = mean(sp_AND_true_cost(i));
AND_true_diff = sp_AND_true_cost_mean - rrt_AND_true_cost_mean;
fprintf('RRT True: %.2f, SL True: %.2f, Raw Improvement: %.2f, %% Improvement: %.2f\n',...
    rrt_AND_true_cost_mean, sp_AND_true_cost_mean, AND_true_diff , AND_true_diff/sp_AND_true_cost_mean);

rrt_AND_length_mean = mean(rrt_AND_length);
rrt_AND_disc_dist_mean = mean(rrt_AND_disc_dist);
sp_AND_disc_dist_mean = mean(sp_AND_disc_dist);
rrt_AND_per_disc_mean = mean(rrt_AND_per_disc);
sp_AND_per_disc_mean = mean(sp_AND_per_disc);
%%
save(strcat('fig2_ws_',date));
%%
function disc_distance = calcDiscDistance(full_cost, dist, eps)
    disc_distance = (full_cost - dist*eps);
    if abs(disc_distance) < 0.25
        disc_distance = 0;
    end
end