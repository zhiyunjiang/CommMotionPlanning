%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% multipath fading generator
% This is generated using the filtering approach to generate correlated
% inphase and quadrature components assuming uniform angle of arrival, 
% which results in Jake's Bessel function for correlation.
%
% usage 
%-----------------------
% gamma_mp_lin = generate_multipath(lambda, g_x, g_y, res, K_ric, corr_mp)
%
% inputs 
%-----------------------
% lambda    : the wavelength of transmission (in meter)
% g_x, g_y  : 2D grid 
% res       : resolution in (samples / m)
% K_ric     : Rician K factor (k = 0 results in Rayleigh distribution)
% corr_mp   : corr_mp = 1 -> correlated multipath and corr_mp = 0 -> uncorrelated multipath
%
% outputs
%-----------------------
% gamma_mp_lin : multipath fading component in linear domain
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

function  gamma_mp_lin = generate_multipath(lambda, g_x, g_y, res, K_ric, corr_mp)

[M, N] = size(g_x);

% We generate this because for multipath fading we need relative distance,
% not distance between Tx and Rx.
[I, J] = meshgrid(1:N, 1:M);        

I = I - round((N+1)/2);
J = J - round((M+1)/2);

d = sqrt(I.^2 + J.^2)/res;

% This is the impulse response of the filter that is applied to the uncorrelated 
% Gaussians (used in corr_gaussian_2D to generate correlated Gaussian).
% generate a Bessel function (as shown in Fig. 3.5 of Goldsmith) in 2D
h = besselj(0, 2*pi*d/lambda);          

% inphase and quadrature components

if corr_mp == 1
    % correlated 
    a = corr_gaussian_2D(0, 1, M, N, h);
    b = corr_gaussian_2D(0, 1, M, N, h);
else
    % uncorrelated 
    a = normrnd(0, 1, M, N);
    b = normrnd(0, 1, M, N);
end

% Generating the correlated Rician
% The standard way to generate a Rician is
% to have R=sqrt(X^2 + Y^2) where X~(v,s^2) and Y~(0,s^2), where v and s
% can be found as a function of K_ric in page 79 of Goldsmith.
% Note that the average of resulting multipath fading power is one, as it should be.
v = sqrt(K_ric/(K_ric + 1));
s = 1/sqrt(2*(1+K_ric));

a = a*s + v;
b = b*s;

gamma_mp_lin = a.^2 + b.^2;

%gamma_mp_dB = 10*log10(gamma_mp_lin);


