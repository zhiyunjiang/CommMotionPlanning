%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel Simulator
% usage
%-------------------------
% [gamma_TOT_dB, gamma_PL_SH_dB, gamma_PL_dB, gamma_SH_dB, gamma_MP_LIN, g_x, g_y] 
% = channel_simulator(region, q_b, K_PL, n_PL, alpha, beta, N_sin, PSD_at_f_c, lambda, K_ric, res, corr_mp)
% 
% inputs
%-------------------------
% region     : a vector containing the boundary of the target rectangular region 
% q_b        : position of the base station (remote station or transmitter)
% K_PL, n_PL : path-loss parameters (gamma_PL_dB = K_PL - 10*n_PL * log10(d)) 
% alpha      : power of the shadowing (in dB)
% beta       : decorrelation distance (in meter)
% N_sin      : # of sinusoids to use 
% PSD_at_f_c : the amplitude difference (in dB) between PSD of shadowing at
%              cutoff frequency and at frequency 0
%              A more detailed description of this variable can be found in
%              the Documentation.
% lambda     : the wavelength of transmission (in meter)
% K_ric      : Rician K factor (K_ric = 0 results in Rayleigh distribution.)
% res        : the resolution of the grid (in samples/m) = 1/(dis between 2 sample along x or y axis)
% corr_mp    : corr_mp = 1 -> correlated multipath and corr_mp = 0 -> uncorrelated multipath
% 
% outputs
%-------------------------
% gamma_TOT_dB   : total channel (path loss + shadowing + multipath) (dB)
% gamma_PL_SH_dB : path loss + shadowing (dB)
% gamma_PL_dB    : path-loss only (dB) 
% gamma_SH_dB    : shadowing only (dB)
% gamma_MP_LIN   : multipath fading component in linear domain 
% g_x, g_y       : the corresponding 2D grid 
%
% Original codes courtesy of Alireza Ghaffarkhah and Alejandro Gonzalez-Ruiz
% Updated by Herbert Cai (April 2016)
%
% If you have any questions or comments regarding the codes,
% please contact Herbert Cai at hcai@ece.ucsb.edu
%
% Dept. of Electrical and Computer Engineering
% University of California, Santa Barbara
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [gamma_TOT_dB, gamma_PL_SH_dB, gamma_PL_dB, gamma_SH_dB, gamma_MP_LIN, g_x, g_y] = channel_simulator(region, ...
    q_b, K_PL, n_PL, alpha, beta, N_sin, PSD_at_f_c, lambda, K_ric, res, corr_mp)


% the rectangular region specifying the environment
x_max = region(1);
x_min = region(2);
y_max = region(3);
y_min = region(4);

% the position of the base station
x_b = q_b(1);
y_b = q_b(2);
                         

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating the grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[g_x, g_y] = meshgrid(x_min:1/res:x_max,y_min:1/res:y_max);
[M, N] = size(g_x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path loss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\ngrid size = %d pixels\n', M*N*res)
fprintf('generating path loss...\n')

g_d = sqrt((g_x - x_b).^2 + (g_y - y_b).^2);

% prevent the path loss to become very large if the samples are very close
% to the base station
g_d(g_d < 1/res) = 1/res;

% generating the path loss
gamma_PL_dB = K_PL - 10*n_PL*log10(g_d);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% shadowing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('generating shadowing...\n')
gamma_SH_dB = generate_shadowing(alpha, beta, N_sin, PSD_at_f_c, g_x, g_y, res);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multipath fading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('generating multipath fading...\n')
gamma_MP_LIN = generate_multipath(lambda, g_x, g_y, res, K_ric, corr_mp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% overall channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gamma_PL_SH_dB = gamma_PL_dB + gamma_SH_dB;
gamma_PL_SH_LIN = 10.^(gamma_PL_SH_dB / 10);
gamma_TOT_LIN = gamma_PL_SH_LIN.*gamma_MP_LIN;
gamma_TOT_dB = 10*log10(gamma_TOT_LIN);





