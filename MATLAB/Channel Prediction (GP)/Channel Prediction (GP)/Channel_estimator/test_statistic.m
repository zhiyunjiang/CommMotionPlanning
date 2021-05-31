% clc
% clear all
% close all
% routes_par=[];
% R2: CNL1 with the following routes [1 3 8 12 13 15 18 21]
% R4: CNL5 with the following routes [2]
% R1: CNL4 with the following routes [1 3 4 7 8]
% R3: CNL2 with the following routes [1 5 7 10 12]
% x_mat = [];
% y_mat = [];
% d_mat = [];
% sig_mat = [];
% for i = [5]
%     in = CNL_gen_2(i);
%     x_mat = [x_mat;in{1}(:)];
%     y_mat = [y_mat;in{2}(:)];
%     d_mat = [d_mat;in{3}(:)];
%     sig_mat = [sig_mat;in{4}(:)];
% end
x_b = 5;
y_b = 5;
x_mat = g_x;
y_mat = g_y;
d_mat = sqrt((g_x - x_b).^2 + (g_y - y_b).^2);
sig_mat = gamma_TOT_dB;

d_mat(find(d_mat)==0) = 0.01;
%-----------------channel estimation---------------------------------------
% commented Yuan Yan Jan. 29th 2013
% editted by Yuan Yan Jan. 20th 2013 based on Mehrzad's original
% test_statistics code
%--------------------------------------------------------------------------

% weight_filter is used to improve the performance of the fit of the
% correlation function of shadowing, because we have less measurements for
% larger distances
global weight_filter
weight_filter = 1;

% per is percentage of measurements
% res is the distance between two samples (meter)
% most of our routes (experiment data) have resolution of 1cm
per = 0.05;
res = .2;

alpha = 20;
beta = 10;
est_PL_par = [-12.89; 3];
[rho,~] = calc_power_mp_given_K_ric(10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%PL estimation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[M, N] = size(sig_mat);
S = round(M*N*per);
idx_sample = sort(randsample(M*N, S));
L = length(sig_mat(:));

% use polyfit2D to extract path loss parameters by using LS estimation
% use ML is not feasible due to the computational efficiency, so we use LS
% here.

PL_matrix = [ones(L,1) -10*log10(d_mat(:))];
est_PL_par = polyfit2D(d_mat(idx_sample), sig_mat(idx_sample));

% estimation of the path loss component
est_PL_comp = PL_matrix*est_PL_par;
% shadowing and MP components
sig_SHMP = sig_mat(:) - est_PL_comp(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%shadowing and MP estimation%%%%%%%%%%%%%%%%%%%%
% function [alpha, beta, rho] = SHMP_par_estimation(x_sample, y_sample, sig_SHMP, res, per, thr_sample_no, thr_SNR)
% 
% this function uses the channel samples to extract the channel parameters

[alpha, beta, rho] = SHMP_par_estimation(x_mat(idx_sample), y_mat(idx_sample), sig_SHMP(idx_sample), res, per, 40, 5);

% function out_MSE = GP2D_tot(ind_sample, sig, est_PL, x_axis, y_axis, alpha, beta, rho, M, N)
%
% this function is used to estimated the channel qualities at all the
% univisted locations in the workspace

[est_mean, est_var] = GP2D_tot_low_memory(idx_sample, sig_mat(idx_sample), est_PL_comp, x_mat, y_mat, alpha, beta, rho, M, N);

figure
surf(sig_mat);
figure
surf(est_mean);


