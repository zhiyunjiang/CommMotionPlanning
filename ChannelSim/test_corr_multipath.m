%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing the ACF of the multipath
% 
% In this file, the ACF of the inphase and quadrature components of the multipath
% fading are calculated and compared with the desired (theoretical) ones.
%
% References used in the comments =========================================
% 1. "Goldsmith":
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
% parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the rectangular region specifying the environment
x_max = 1;
x_min = -1;
y_max = 1;
y_min = -1;

% Multipath fading parameters
% lambda: the wavelength of transmission (in meter)
% res: the resolution of the grid (in samples/m) = 1/(dis between 2 samples along x or y axis)
% K_ric: parameter of the Rician distribution (see Goldsmith page 79)

% lambda = 0.125 corresponds to a carrier frequency of 2.4 GHz.
lambda = 0.125;

% ss_decorr is the decorrelation distance for multipath fading.
% This is the point where J_0(2*pi*ss_decorr/lambda) is equal to 0 (see Eq. 3.26 of Goldsmith).
% Note that ss_decorr = 0.4*lambda.
% This variable is not an input to the channel simulator, but is used to 
% to provide a guideline for choosing the simulation resoltuion.
ss_decorr = 0.05;         

% This makes 1/res = ss_decorr/10 to better visualize the correlation.
res = 10/ss_decorr; 

K_ric = 10;  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating the grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[g_x, g_y] = meshgrid(x_min:1/res:x_max,y_min:1/res:y_max);
[M, N] = size(g_x);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating correlated multipath (see function generate_multipath.m for
% the details and philosophy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the impulse response of the filter

[I, J] = meshgrid(1:N, 1:M);

I = I - round((N+1)/2);
J = J - round((M+1)/2);

d = sqrt(I.^2 + J.^2)/res;

h = besselj(0, 2*pi*d/lambda);  

% inphase and quad components
a = corr_gaussian_2D(0, 1, M, N, h);
b = corr_gaussian_2D(0, 1, M, N, h);


% 2D correlation of the inphase and quad components
C_inph = xcorr2(a)/(M*N);       %calculating the autocorrelation matrix of the inphase part
C_quad = xcorr2(b)/(M*N);       %calculating the autocorrelation matrix of the quadrature part


% to make the calculations quicker we just pick a slice of the field 
% specifically, we just consider the elements along x and y axis, 
% otherwise the calculation takes too long
[P, Q] = size(C_inph);

% ACF in x dir
ACF_inph_x = C_inph((P+1)/2,(Q+1)/2:end);       %pick the middle
ACF_quad_x = C_quad((P+1)/2,(Q+1)/2:end);   
d_x = (0:length(ACF_inph_x) - 1)/res;

% ACF in y dir
ACF_inph_y = C_inph((P+1)/2:end, (Q+1)/2);      %pick the middle
ACF_quad_y = C_quad((P+1)/2:end, (Q+1)/2);
d_y = (0:length(ACF_inph_y) - 1)/res;



% the theoretical ACF
ACF_T_x = besselj(0, 2*pi*d_x/lambda);  
ACF_T_y = besselj(0, 2*pi*d_y/lambda);  


%**************************************************************************
% show the plots
f = figure;
subplot(2,2,1)
plot(d_x, ACF_inph_x, '-b', 'LineWidth', 5);
hold on
plot(d_x, ACF_T_x, '--r', 'LineWidth', 5);
ylabel('ACF of inph (x dir)', 'FontSize', 30)
legend('Sim', 'Theory')
set(gca, 'FontSize', 30)
axis tight

subplot(2,2,2)
plot(d_y, ACF_inph_y, '-b', 'LineWidth', 5);
hold on
plot(d_y, ACF_T_y, '--r', 'LineWidth', 5);
ylabel('ACF of inph (y dir)', 'FontSize', 30)
legend('Sim', 'Theory')
set(gca, 'FontSize', 30)
axis tight

subplot(2,2,3)
plot(d_x, ACF_quad_x, '-b', 'LineWidth', 5);
hold on
plot(d_x, ACF_T_x, '--r', 'LineWidth', 5);
xlabel('d (m)', 'FontSize', 30)
ylabel('ACF of quad (x dir)', 'FontSize', 30)
legend('Sim', 'Theory')
set(gca, 'FontSize', 30)
axis tight

subplot(2,2,4)
plot(d_y, ACF_quad_y, '-b', 'LineWidth', 5);
hold on
plot(d_y, ACF_T_y, '--r', 'LineWidth', 5);
xlabel('d (m)', 'FontSize', 30)
ylabel('ACF of quad (y dir)', 'FontSize', 30)
legend('Sim', 'Theory')
set(gca, 'FontSize', 30)
axis tight

maximize(f)







