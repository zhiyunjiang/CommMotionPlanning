% Simple Tour Adjustment

%Set up a simple workspace consisting of three points of interest in
%addition to the starting position.
%For now, no obstalces

%Setup Comm Channel Params
%Using parameters chosen for Comm Aware RRT* paper
[n_PL, K_PL, sh_decorr, sigma_SH, PSD_at_f_c, mp_decorr, lambda, K_ric,...
    corr_mp, sigma_mp, region, N_sin, res, n_samples, n_sims, gamma_TH] = getConstParams(); 

%Setup start location, other figures
q_start = [1, 1];
q1 = [47,3];
cp1 = ChannelParams(q1, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);
 
q2 = [45, 49];
cp2 = ChannelParams(q2, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);

q3 = [4, 46];
cp3 = ChannelParams(q3, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);

% generate the full channel to sample
%cc1 = CommChannel(cp1, N_sin, region, res);
%cc1.simulateSH(); cc1.simulateMP(MP_Mdl);

cc2 = CommChannel(cp2, N_sin, region, res);
cc2.simulateSH(); cc2.simulateMP(2);

% cc3 = CommChannel(cp3, N_sin, region, res);
% cc3.simulateSH(); cc3.simulateMP(MP_Mdl);


%Sample
% [obs_pos1, obs_vals1] = cc1.randSample(n_samples);
[obs_pos2, obs_vals2] = cc2.randSample(n_samples);
% [obs_pos3, obs_vals3] = cc3.randSample(n_samples);


%Create Predicted Channel
% pc1 = PredictedChannel(cp1, cc1, obs_pos1, obs_vals1);
pc2 = PredictedChannel(cp2, cc2, obs_pos2, obs_vals2);
% pc3 = PredictedChannel(cp3, cc3, obs_pos3, obs_vals3);

%setup the problem
obs_mod = ObstacleFactory.noObs();
%for now, assume we only make one-step increments
cost_fxn = @(n1, n2, path, mode) MinExpFPDInd(n2, path, pc2, gamma_TH);
%Run RRT for each separately, so we have nice trees for each
%src will be the start location.
%goal could really be anything, just make it far away
%for now, just plot probabilities
pc2.plotPosteriors2D(gamma_TH);

problem_instance = PathPlanningProblem(region, res, q2, [0,0], obs_mod, cost_fxn);

rrt_solver = getRRTSolver();
rrt_solver.solve(problem_instance);
root = rrt_solver.getRoot();
hold on
root.plotTree('k');

