%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel_simulator test file
% A sample channel is generated and tested to be a good match (first and second 
% moments as well as the PDF are tested). The correlation functions are tested
% in two separate scripts: test_corr_shadow.m and test_corr_multipath.m.
% (See the Documenation file for an overview of the test scripts.)
%
% References used in the comments =========================================
% 1. "JOR":(by Gonzales-Ruiz et al.)
% A Comprehensive Overview and Characterization of Wireless Channels for Networked Robotic and Control Systems 
% This is our paper published in the Journal of Robotics, which has the
% mathematical details of the simulation (Section 6).
%
% 2. "Goldsmith":
% Wireless Communications by A. Goldsmith
% This is the textbook that we will often refer to.
% =========================================================================
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

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The rectangular region specifying the environment, indicating the upper
% right and lower left corners of the workspace
x_max = 6;
x_min = -4;
y_max = 5.5;
y_min = -5;
region = [x_max x_min y_max y_min];

% Position of the base station (remote station or transmitter)
q_b = [0 0];

% Path loss parameters
% The resulting path loss component in dB is 
% gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
% K_PL is the path loss constant in dB and n_PL is the path loss exponent.
n_PL = 3;  
K_PL = -12.89;    


% Shadowing parameters
% alpha : power of the shadowing (in dB)
% beta : decorrelation distance (in meter). The shadowing correlation
%        model we use is: alpha*exp(-distance/beta).
% N_sin : # of sinusoids to use, this is arbitrary. The higher the number
%         the more accurate the generation of the correlated shadowing will
%         be. However, if it is too high, the computational cost is
%         significant. We have seen that 5000 provides a good balance
%         between accuracy and computational time.
% PSD_at_f_c : The amplitude difference (in dB) between PSD of shadowing at
%              cutoff frequency and at frequency 0, which is used to find
%              the cutoff frequency. The cutoff frequency is used to approximate
%              the desired shadowing PSD, which has a "low pass" shape
%              (see JOR Eq. 24 for the mathematical expression of the shadowing
%              PSD). The portion of the PSD beyond the cutoff frequency will 
%              be ignored in the approximation.
%              See the JOR paper (Section 6) for more details.

alpha = 10;            
beta = 2;            
N_sin = 5000;                    
PSD_at_f_c = 30;                  


% Multipath fading parameters
% lambda: the wavelength of transmission (in meter)
% res = the resolution of the grid (in samples/m) = 1/(dis between 2 samples along x or y axis)
% K_ric = parameter of the Rician distribution (see Goldsmith page 79)

% lambda = 0.125 corresponds to a carrier frequency of 2.4 GHz.
lambda = 0.125;

% ss_decorr is the decorrelation distance for multipath fading.
% This is the point where J_0(2*pi*ss_decorr/lambda) is equal to 0 (see Eq. 3.26 of Goldsmith).
% Note that ss_decorr = 0.4*lambda.
% This variable is not an input to the channel simulator, but is used to 
% to provide a guideline for choosing the simulation resoltuion.
ss_decorr = 0.05;         

% This makes 1/res = ss_decorr/2.
% Make sure 1/res < ss_decorr so that we can see some correlation.
res = 2/ss_decorr; 

K_ric = 10;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('generating the channel...\n')


% corr_mp = 1 -> correlated multipath 
% corr_mp = 0 -> uncorrelated multipath 
corr_mp = 0;
[gamma_TOT_dB, gamma_PL_SH_dB, gamma_PL_dB, gamma_SH_dB, gamma_MP_LIN, g_x, g_y] = channel_simulator(region, ...
    q_b, K_PL, n_PL, alpha, beta, N_sin, PSD_at_f_c, lambda, K_ric, res, corr_mp);


fnt_siz = 25;

f = figure;


surf(g_x, g_y, gamma_TOT_dB, 'EdgeColor','none');
% light
% shading interp
xlabel('x (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
ylabel('y (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
zlabel('Received power (PL + SH + MP) (dBm)','FontSize', fnt_siz ,  'FontWeight','bold');
axis tight
grid on
set(gca, 'FontSize', fnt_siz, 'FontWeight', 'bold');

maximize(f)


f = figure;
surf(g_x, g_y, gamma_PL_SH_dB, 'EdgeColor','none');
% light
% shading interp
xlabel('x (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
ylabel('y (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
zlabel('Received power (PL + SH) (dBm)','FontSize', fnt_siz,  'FontWeight','bold');
axis tight
grid on
set(gca, 'FontSize', fnt_siz, 'FontWeight', 'bold');

maximize(f)


f = figure;
surf(g_x, g_y, gamma_SH_dB, 'EdgeColor','none');
% light
% shading interp
xlabel('x (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
ylabel('y (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
zlabel('Received power (SH only) (dBm)','FontSize', fnt_siz,  'FontWeight','bold');
axis tight
grid on
set(gca, 'FontSize', fnt_siz, 'FontWeight', 'bold');

maximize(f)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test the Gaussian match for the generated 
% shadowing component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\ntesting the gaussian match for the generated shadowing...\n')


% Fit a normal dist. to the generated shadowing 
[mean_SH_num, sigma_SH_num] = normfit(gamma_SH_dB(:));

fprintf('\nThe mean of the generated shadowing component in dB : %5.2f\n', mean_SH_num)
fprintf('The var of the generated shadowing component and its desired value: %5.2f, %5.2f\n', sigma_SH_num^2, alpha)


f = figure;
histfit(gamma_SH_dB(:), 50)
title('pdf match for the generated shadowing','FontSize', fnt_siz,  'FontWeight','bold');
set(gca, 'FontSize', fnt_siz, 'FontWeight', 'bold');

maximize(f)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test the Rician match for the generated 
% multipath component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\ntesting the rician match for the generated multipath...\n')

gamma_MP_dB = 10*log10(gamma_MP_LIN);

% theoretical (desired) values of the mean and variance of multipath component in dB
[sigma_2_MP_theo, mean_MP_theo] = calc_power_mp_given_K_ric(K_ric);

% numerical values of the mean and variance of multipath component in dB
mean_MP_num = mean(gamma_MP_dB(:));
sigma_2_MP_num = var(gamma_MP_dB(:));

fprintf('\nThe mean of the generated multipath component in dB and its desired value : %5.2f, %5.2f\n', mean_MP_num, mean_MP_theo)
fprintf('The var of the generated multipath component and its desired value: %5.2f, %5.2f\n', sigma_2_MP_num, sigma_2_MP_theo)

f = figure;
histfit(gamma_MP_LIN(:), 50, 'rician')
title('pdf match for the generated multipath','FontSize', fnt_siz,  'FontWeight','bold');
set(gca, 'FontSize', fnt_siz, 'FontWeight', 'bold');

maximize(f)




                         