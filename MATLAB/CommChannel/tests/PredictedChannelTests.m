%PredictedChannelTests
%%
%Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Connectivity Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position of the base station (remote station or transmitter)
%q_b = [3 5.5];
q_b = [-0.2, 1.345];
% Path loss parameters
% The resulting path loss component in dB is 
% gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
% K_PL is the path loss constagrd_sznt in dB and n_PL is the path loss exponent.
n_PL = 4.2;  
K_PL = -45;

%Shadowing Parameters
%Shadowing decorrelation distance
sh_decorr = 6;    
%shadowing power
alpha = 8.41;
%standard deviation of shadowing and power conn
sigma_SH = sqrt(alpha);
PSD_at_f_c = 30;

%Multipath Parameters
%Multipath decorrelation, used to set grid resolution
%not relevant if modeling mp as uncorrelated
mp_decorr = 0.05;         

lambda = 0.125;
K_ric = 1.59;     

% corr_mp = 1 -> correlated multipath 
% corr_mp = 0 -> uncorrelated multipath 
corr_mp = 0;

%
sigma_mp = 0;%=0 -> no mp when modeling MP fading as zero-mean, lognormal

cp = ChannelParams(q_b, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Comm Channel Simulation & Generation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%keep it small for testing
x_max = 30;
x_min = -10;
y_max = 10;
y_min = -30;
region = [x_max x_min y_max y_min];
%simulation parameter
N_sin = 5000;

res = 5;
cc = CommChannel(cp, N_sin, region, res);
cc.simulateSH(); cc.simulateMP();
%plot the true channel for comparison
cc.plot(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the Channel Analyzer With Observations, View Posterior Probability of Connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sample the "true" channel
n_samples = 100;
[sample_pos, sample_vals] = cc.randSample(n_samples);

gamma_TH = -90;
pc = PredictedChannel(cp, cc, sample_pos, sample_vals);
%Test plots
figure()
receiver_noise = 1e-10; R = 2; BER = 1e-6;
qos = QoSParams(BER, R, receiver_noise);
req_power = pc.getEReqTXPowerW( qos);
pc.plotRequiredTXPower2D(req_power);
figure()
pc.plotMeans2D();
figure()
pc.plotConnected2D(0.8, gamma_TH);
figure()
pc.plotPosteriors2D(gamma_TH);

%%
% Runt tests for this predicted channel
test_count = 0;
fail_count = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% scaledBeta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_count  = test_count + 1;
if res*sh_decorr ~= pc.scaledBeta()
    fail_count = fail_count + 1;
    fprintf('Failed scaledBeta test');
end





