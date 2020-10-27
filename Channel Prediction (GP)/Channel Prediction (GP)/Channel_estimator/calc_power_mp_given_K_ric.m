%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculating mean and var of rician in dB v 05162010
%
% usage 
%-----------------------
% [sigma_2_mp, mean_mp] = calc_power_mp_given_K_ric(K_ric)
%
% The mean and var are calulated by direct integration 
%
% inputs
%-------------------------
% K_ric      : rician K factor (k = 0 results in Rayleigh dist)
% 
% outputs
%-------------------------
% sigma_2_mp : the var of rician fading in dB domain
% mean_mp    : the mean of rician fading in dB domain
%
% Copyright 2010 Alireza Ghaffarkhah (alinem@ece.unm.edu) and Alejandro
% Gonzalez-Ruiz (agon@unm.edu)
% Cooperative Network Lab, 
% Dept. of Electrical and Computer Engineering
% University of New Mexico, Albuqueurque NM, 87131
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [sigma_2_mp, mean_mp] = calc_power_mp_given_K_ric(K_ric)

x = [0.01:0.01:10];

% the pdf of rician fading in lin domain
p = (1+K_ric)*exp(-K_ric - (K_ric + 1)*x).* besseli(0, 2*sqrt(x*K_ric*(K_ric + 1)));

% trapz integration 
y1 = 100*(log10(x)).^2.*p;
y2 = 10*log10(x).*p;

T1 = trapz(x,y1);
T2 = trapz(x,y2);

mean_mp = T2;
sigma_2_mp = T1 - T2^2;