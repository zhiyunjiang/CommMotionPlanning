%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shadow fading generator
%
% usage 
%-----------------------
% gamma_sh_dB = generate_shadowing(alpha, beta, N, PSD_at_f_c, g_x, g_y, res)
% 
% inputs 
%-----------------------
% alpha      : power of the shadowing (in dB)
% beta       : decorrelation distance (in meter)
% N          : # of sinusoids to use 
% PSD_at_f_c : the amplitude difference (in dB) between PSD of shadowing at
%              cutoff frequency and at frequency 0
%              A more detailed description of this variable can be found in
%              the Documentation.
% g_x, g_y   : 2D grid 
% res        : resolution in (samples / m)
%
% outputs
%-----------------------
% gamma_sh_dB : shadowing component (dB)
%
% References used in the comments =========================================
% This code follows the method described by:
% [1] X. Cai and G. B. Giannakis, “A two-dimensional channel
% simulation model for shadowing processes,IEEE Transactions
% on Vehicular Technology, vol. 52, no. 6, pp. 1558-1567, 2003.
%
% "JOR":(by Gonzales-Ruiz et al.)
% A Comprehensive Overview and Characterization of Wireless Channels for 
% Networked Robotic and Control Systems 
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

function gamma_sh_dB = generate_shadowing(alpha, beta, N, PSD_at_f_c, g_x, g_y, res)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generating the shadowing component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% standard deviation of shadowing
sigma = sqrt(alpha);        

a = 1/beta;

% M is the number of radial frequencies, which is related to the number of sinusoids through: N = 2*M^2;
M = round(sqrt(N/2));       
                            
% For mathematical details see [1] or the JOR paper (paper info. at the beginning of this file)

PSD_at_f_c_lin = 10^(-PSD_at_f_c/10);                               %see definition of epsilon on page 16 of JOR
f_c = sqrt( 1/(4*pi^2)*((a^3/(PSD_at_f_c_lin))^(2/3)-a^2) );        %f_c is the cuttoff frequency, see definition of f_{r,c} on pager 16 of JOR
P = 1 - a/(sqrt(a^2+4*pi^2*f_c^2));                                 %P is the mean power within the cutoff frequency, we calculated this expression from Eq. 15 in [1]
cn = sqrt(2/N);                                                     %coefficient (c_j) in Eq. 28 of JOR, the value coms from page 4 of [1]


frm = zeros(1, M+1);  

% compute the frequency frm recursively (Eq. 17 of [1])
for m = 2:M+1
    frm(m) = 1/(2*pi)*sqrt((P/(M*a)-1/sqrt(a^2 + 4*pi^2*frm(m-1)^2))^(-2)-a^2);  
end


delta_phi = 2*pi/(4*M);         %see explanation below Eq. 29 of JOR

% 2M angles uniformly distributed in ( -pi/2 , pi/2 )
phi = linspace(-pi/2, pi/2-delta_phi/2 , 2*M);               

fxn = zeros(M, 2*M);        %sampling frequencies along x and y
fyn = zeros(M, 2*M);

for m = 1:M
    fxn(m,:) = (frm(m) + frm(m+1)) * cos(phi) /2;
    fyn(m,:) = (frm(m) + frm(m+1)) * sin(phi) /2;
end

theta = rand(size(fxn))*2*pi;           %the phase should be uniform between 0 and 2*pi

sn = zeros(size(g_x));

for i = 1:size(fxn,1)-1
    for j = 1:size(fxn,2)-1
        sn =  sn + cn* cos(2*pi*(fxn(i,j)*g_x + fyn(i,j) *g_y) + theta(i,j));      %this is generating the actual shadowing 
                                                                                   %by summing the sinusoids with appropriate frequencies and phase (see Eq. 28 of JOR)
    end
end

gamma_sh_dB = (sn-mean(sn(:)))/std(sn(:))*sigma;            %giving the shadowing the appropriate power
