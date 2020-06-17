%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculating mean and variance of Rician in dB
%
% usage 
%-----------------------
% [sigma_2_mp, mean_mp] = calc_power_mp_given_K_ric(K_ric)
%
% The mean and variance are calculated by direct integration 
%
% inputs
%-------------------------
% K_ric      : rician K factor (k = 0 results in Rayleigh dist.)
% 
% outputs
%-------------------------
% sigma_2_mp : the variance of Rician fading in dB domain
% mean_mp    : the mean of Rician fading in dB domain
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

function [sigma_2_mp, mean_mp] = calc_power_mp_given_K_ric(K_ric)

x = [0.01:0.01:10];

% the PDF of Rician fading in linear domain
p = (1+K_ric)*exp(-K_ric - (K_ric + 1)*x).* besseli(0, 2*sqrt(x*K_ric*(K_ric + 1)));

% trapz integration 
y1 = 100*(log10(x)).^2.*p;
y2 = 10*log10(x).*p;

%trapz does a numerical integration via the trapezoidal method
T1 = trapz(x,y1);
T2 = trapz(x,y2);

mean_mp = T2;
sigma_2_mp = T1 - T2^2;