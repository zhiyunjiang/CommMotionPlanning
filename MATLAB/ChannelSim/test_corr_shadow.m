%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing the ACF of the shadowing
%
% Two methods are used for calculating the ACF of the resulting shadowing
% signal. One method (the faster one) uses the xcorr2 function of MATLAB to
% calculate the 2D ACF of the signal and considers the values along both x
% and y axis. The second method (the slower but more accurate one) uses a
% Monte Carlo method to calculate the ACF. In the Monte Carlo method, a set
% of random samples are chosen and the value of the ACF at distance d is 
% calculated by considering all the points that are on a circle (w/ radius 
% d) around each random sample. 
%
% References used in the comments =========================================
% 1. "JOR":(by Gonzales-Ruiz et al.)
% A Comprehensive Overview and Characterization of Wireless Channels for Networked Robotic and Control Systems 
% This is our paper published in the Journal of Robotics, which has the
% mathematical details of the simulation (Section 6).
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
x_max = 20;
x_min = -10;
y_max = 40;
y_min = 10;

% res = the resolution of the grid (in samples/m) = 1/(dis between 2 sample
% along x or y axis)
res = 5; 


% Shadowing parameters
% alpha : power of the shadowing (in dB)
% beta : decorrelation distance (in meter). The shadowing correlation
%        model we use is: alpha*exp(-distance/beta)
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
alpha = 16;            
beta = 5;            
N_sin = 5000;                    
PSD_at_f_c = 30;           



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating the grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[g_x, g_y] = meshgrid(x_min:1/res:x_max,y_min:1/res:y_max);
[M, N] = size(g_x);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test the correlation of shadowing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



disp('calculating the ACF using xcorr2 function of MATLAB ...')



gamma_sh_dB = generate_shadowing(alpha, beta, N_sin, PSD_at_f_c, g_x, g_y, res);
C = xcorr2(gamma_sh_dB)/(M*N);          %autocorrelation for the whole field  

    
% Since we only want to show a 1D plot, we just consider one row along 
% x-axis and one column along y-axis.

[P, Q] = size(C);

% ACF in x dir
ACF_x = C((P+1)/2,(Q+1)/2:end);     %pick the middle row
d_x = (0:length(ACF_x) - 1)/res;

% ACF in y dir
ACF_y = C((P+1)/2:end, (Q+1)/2);    %pick the middle column
d_y = (0:length(ACF_y) - 1)/res;

% the theoretical ACF
ACF_T_x = alpha*exp(-d_x./beta);
ACF_T_y = alpha*exp(-d_y./beta);


pause(0.2)

% in case one wants a more precise comparison of ACFs. 
% Here ACF_p(i) is the precise value of the ACF at distance d_p(i), which is 
% calculated by considering all the samples that are at distance d_p(i) 
% and averaging over them

if length(d_x) > length(d_y)        %we pick the longest axis
    d_p = d_x;
else
    d_p = d_y;
end;

% Since ACF_2D takes a while to return the ACF, it is desirable to 
% average only over a few samples in a Monte Carlo type approach.
% The parameter 'perc' controls the percentage of the samples used
% for averaging

disp('calculating the precise ACF using Monte Carlo method ...')
perc = 10;
ACF_p = ACF_2D(gamma_sh_dB, g_x, g_y, d_p, res, perc);
ACF_T_p = alpha*exp(-d_y./beta);


f = figure;
subplot(3,1,1)
plot(d_x, ACF_x, '-b', 'LineWidth', 5);
hold on
plot(d_x, ACF_T_x, '--r', 'LineWidth', 5);
ylabel('ACF (x dir)', 'FontSize', 30)
legend('Sim', sprintf('Theory, alpha = %5.2f, beta = %5.2f m', alpha, beta))
set(gca, 'FontSize', 30)

subplot(3,1,2)
plot(d_y, ACF_y, '-b', 'LineWidth', 5);
hold on
plot(d_y, ACF_T_y, '--r', 'LineWidth', 5);
ylabel('ACF (y dir)', 'FontSize', 30)
legend('Sim', sprintf('Theory, alpha = %5.2f, beta = %5.2f m', alpha, beta))
set(gca, 'FontSize', 30)

subplot(3,1,3)
plot(d_x, ACF_x, '-b', 'LineWidth', 5);
hold on
plot(d_x, ACF_T_x, '--r', 'LineWidth', 5);
xlabel('d (m)', 'FontSize', 30)
ylabel('ACF (precise)', 'FontSize', 30)
legend('Sim', sprintf('Theory, alpha = %5.2f, beta = %5.2f m', alpha, beta))
set(gca, 'FontSize', 30)









