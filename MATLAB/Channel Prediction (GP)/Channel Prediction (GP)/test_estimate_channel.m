% test estimation code
clear;
clc;

save_results = 1;

load('channel.mat');

channel_true = gamma_TOT_dB;
num_loc = length(channel_true(:));
perc_obs = 0.01;
num_obs = ceil(perc_obs*num_loc);
idx_o = randsample(1:num_loc,num_obs);
idx_u = setdiff(1:num_loc,idx_o);

x_o = g_x(idx_o);
y_o = g_y(idx_o);
z_o = channel_true(idx_o);
x_a = g_x(:);
y_a = g_y(:);


%% prediction
tic

% predict in batches rather than the entire map to avoid jamming the memory
batch = 2000;
curr = 1;
pre_mean = [];
pre_var = [];
for i_batch = 1:ceil(length(x_a)/batch)-1
    x_a_curr = x_a(curr:curr+batch-1);
    y_a_curr = y_a(curr:curr+batch-1);
    [pre_mean_curr, ~, pre_var_curr] = estimate_channel(x_o, y_o, z_o,... % obs. var.
                                                        x_a_curr, y_a_curr,...
                                                        alpha, beta, eta,...
                                                        K_PL, n_PL, q_b);
    pre_mean = [pre_mean;pre_mean_curr];
    pre_var = [pre_var;pre_var_curr];
    
    curr = curr+batch;
end

% last batch
x_a_curr = x_a(curr:end);
y_a_curr = x_a(curr:end);
[pre_mean_curr, ~, pre_var_curr] = estimate_channel(x_o, y_o, z_o,... % obs. var.
    x_a_curr, y_a_curr,...
    alpha, beta, eta,...
    K_PL, n_PL, q_b);
pre_mean = [pre_mean;pre_mean_curr];
pre_var = [pre_var;pre_var_curr];

% the predicted values are in vector form
% reshape them to 2D if needed

toc

%% save
if save_results == 1
    save('channel_pred.mat', 'channel_true','pre_mean','pre_var','region','g_x','g_y','idx_o');
end
