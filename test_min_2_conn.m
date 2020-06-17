%Find Min Dist to Connectivity
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_max = 6;
x_min = -4;
y_max = 5.5;
y_min = -5;
region = [x_max x_min y_max y_min];

source = [-4, 5];
dest = [6, -5];

% maze = [RectObs([-1, -1.5, 3, -5.5]), RectObs([0, -0.5, 6, -2]), ...
%             RectObs([1, 0.5, 3, -5.5]), RectObs([2, 1.5, 6, -2])];
% simple = [CircObs(1.5, [0,4]), RectObs([4, 2, 0, -2 ]), CircObs(1, [2,2])];
% 
% dense = RectObs.empty;
% dx = 0.3;
% dy = 0.3;
% for i=1:floor(x_max - x_min)
%     
%     for j = 1:floor(y_max - y_min)
%         y_jitter = 2*dy*rand(1); 
%         x_jitter = 2*dx*rand(1); 
%         obs = RectObs([i+x_min+dx + x_jitter, i+x_min + x_jitter, j+y_min+dy + y_jitter, j+y_min+y_jitter]);
%         dense((i-1)*floor(y_max - y_min) + j) = obs;
%     end 
% end

obs_mod = ObstacleMod([]);
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Connectivity Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position of the base station (remote station or transmitter)
q_b = [3 5.5];

% Path loss parameters
% The resulting path loss component in dB is 
% gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
% K_PL is the path loss constant in dB and n_PL is the path loss exponent.
n_PL = 3;  
K_PL = -12.89;
%Shadowing decorrelation distance
beta = 2;    

%Multipath decorrelation, used to set grid resolution
ss_decorr = 0.05;         
res = 2/ss_decorr;

%shadowing power
alpha = 10;
%standard deviation of shadowing and power conn
sigma_SH = sqrt(alpha);

N_sin = 5000;                    
PSD_at_f_c = 30;
lambda = 0.125;
K_ric = 10;     

% corr_mp = 1 -> correlated multipath 
% corr_mp = 0 -> uncorrelated multipath 
corr_mp = 0;

[gamma_TOT_dB, ~, gamma_PL_dB, ~, ~, g_x, g_y] = channel_simulator(region, ...
    q_b, K_PL, n_PL, alpha, beta, N_sin, PSD_at_f_c, lambda, K_ric, res, corr_mp);
%%
%connectivity priors
%find probability of being connected - dobule integral in the form of a
%convolution