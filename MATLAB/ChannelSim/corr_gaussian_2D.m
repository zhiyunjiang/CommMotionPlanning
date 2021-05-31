%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel Simulator
% This function generates correlated Gaussian samples
%
% usage 
%-----------------------
% g = corr_gaussian_2D(mu, sigma, M, N, h)
%
% inputs
%-------------------------
% mu, sigma : the mean and std (standard deviation) of the Gaussian 
% M, N      : the desired size of the output matrix
% h         : the impulse response of the filter (used to generate gaussian with arbitrary correlation)
%
% outputs
%-------------------------
% g : the resulting Gaussian matrix with desired correlation
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

function g = corr_gaussian_2D(mu, sigma, M, N, h)

a = normrnd(0, 1, M, N);            %generate an uncorrelated input

% make it correlated by filtering, we multiply the samples in the frequency domain
% with the appropriate filter and transform the result back to space domain
% see Eq. 27 of JOR for details

absH = sqrt(abs(fft2(h)));             
A_f = absH.*fft2(a);                             
a_f = ifft2(A_f);
g = a_f*sigma/std(a_f(:)) + mu;     %it gives the result the desired mean and std

