%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getConstParams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets parameters used in simulations for our Communication Aware RRT*
% Paper. See inline comments for description of outputs

function [n_PL, K_PL, sh_decorr, sigma_SH, PSD_at_f_c, mp_decorr, lambda, K_ric,...
    corr_mp, sigma_mp, region, N_sin, res, n_samples, n_sims, gamma_TH] = getConstParams()

    %Path loss parameters
    % The resulting path loss component in dB is 
    % gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
    % K_PL is the path loss constant in dB and n_PL is the path loss exponent.
    % Based on real workld measurements. See:
    % M. Malmirchegini and Y. Mostofi. On the spatial predictability  of 
    % communication channels. IEEE Transactions on Wireless Comms., 2012
    n_PL = 3.86;
    K_PL = -41.34;
    
    
    %Shadowing Parameters
    % Shadowing decorrelation distance
    sh_decorr = 3.09;   
    % shadowing power
    alpha = 10.24;%+2;
    % standard deviation of shadowing and power conn
    sigma_SH = sqrt(alpha);
    PSD_at_f_c = 30;

    %Multipath Parameters
    % Multipath decorrelation, used to set grid resolution
    % not relevant here since we're modeling as uncorrelated and log-normal
    mp_decorr = 0.05;         
    lambda = 0.125;
    K_ric = 14.56;     
    % corr_mp = 1 -> correlated multipath 
    % corr_mp = 0 -> uncorrelated multipath 
    corr_mp = 0;
    % Note - from measurements, 3.02 is a good value
    % sigma_mp = 1;%=0 -> no mp when modeling MP fading as zero-mean, lognormal
    sigma_mp = 3.02;
    % channel simulation parameter
    N_sin = 5000;
    % minimum required channel gain in dB
    gamma_TH = -101;
    
    %Comm Region
    % 50m x 50m region with 0.2 m resolution gives us a 250x250 grid, i.e. 
    % 62500 elements
    x_max = 50;
    x_min = 0;
    y_max = 50;
    y_min = 0;
    region = [x_max x_min y_max y_min];
    % res = 2/mp_decorr; % good value if working with correlated MP
    res = 5;%0.2 m
    
    % sampling percentage: 500/62500 = 0.8% 
    n_samples = 500;
    
    % How many times to run the experiment for each setting
    n_sims = 20;
    
end

