%Performance and Unit Test for ChannelAnalyzer

%%
%Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Connectivity Characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position of the base station (remote station or transmitter)
%q_b = [3 5.5];
q_b = [-450 0];
% Path loss parameters
% The resulting path loss component in dB is 
% gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
% K_PL is the path loss constant in dB and n_PL is the path loss exponent.
n_PL = 4.2;  
K_PL = 27;

%Shadowing Parameters
%Shadowing decorrelation distance
sh_decorr = 12.92;    
%shadowing power
alpha = 8.41;
%standard deviation of shadowing and power conn
sigma_SH = sqrt(alpha);
PSD_at_f_c = 30;

%Multipath Parameters
%Multipath decorrelation, used to set grid resolution
mp_decorr = 0.05;         

lambda = 0.125;
K_ric = 1.59;     

% corr_mp = 1 -> correlated multipath 
% corr_mp = 0 -> uncorrelated multipath 
corr_mp = 0;

cp = ChannelParams(q_b, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Comm Channel Simulation & Generation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_max = 100;
x_min = 0;
y_max = 100;
y_min = 0;
region = [x_max x_min y_max y_min];
%simulation parameter
N_sin = 5000;

%we're not going to be looking at MP component for a quick minute, so
%setting the resolution to something more meaningful
%res = 2/mp_decorr;
res = 10;
cc = CommChannel(cp, N_sin, region, res);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the Channel Analyzer, Preview Prior Probability of Connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%minimum power required for connectivity
gamma_TH = -80;
no_mp = 0;
ca = ChannelAnalyzer(cc, gamma_TH, no_mp);

%Performance test computing the expected distance
pathy = zeros([100*res, 1]);
pathx = cumsum(ones([100*res, 1]));

path = [flip(pathx), pathy];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Arjun's recursive approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Test Performance of Recursive Markovian Approximation of Expected FPD\n');
exp4 = @() ca.ExpectedFPD(path,4); 
exp4_runtime = timeit(exp4);
fprintf('Expected Runtime for a 100-point path is %5f\n', exp4_runtime);
%%
%verify we accurately calculate straighline-path expected FPD
[approx_PMF_4, distances] = ca.ApproxFPDPMF(path,4);
fprintf('Found PMF 4\n');
approx_CDF_4 = cumsum(approx_PMF_4);

step_size = 1/res;
figure(5)
clf;
hold on
plot(cumsum(distances)*step_size, approx_CDF_4);

%verify by simulating several channels and calculating expected distance
%to connectivity

all_dists = ca.simulateFPD(path, 1000, 1);
histogram(all_dists*step_size, 'Normalization', 'cdf', 'BinEdges',(cumsum(distances))*step_size);
legend('Markovian Approximation', 'sims');

ylabel('Probability')
xlabel('Distance (m)')
title('FPD CDF for Straight Line Path')

%%
%Testing with non-straighline paths

root = [100, 1];
deltax = ones([99, 1]);%deltax(1:2:end) = 0;
deltay = zeros([99, 1]);deltay(2:3:end) = 1;
pathx = 150 - cumsum(deltax);
pathy = 1 + cumsum(deltay);
path = [root;pathx, pathy];

%verify we accurately calculate straighline-path expected FPD
[approx_PMF_4, distances] = ca.ApproxFPDPMF(path,4);
fprintf('Found PMF 4\n');
approx_CDF_4 = cumsum(approx_PMF_4);

step_size = 1/res;
figure(6)
clf;
hold on
plot(cumsum(distances)*step_size, approx_CDF_4);

%verify by simulating several channels and calculating expected distance
%to connectivity

all_dists = ca.simulateFPD(path, 1000, 1);
histogram(all_dists*step_size, 'Normalization', 'cdf', 'BinEdges',(cumsum(distances))*step_size);
legend('Markovian Approximation', 'sims');

ylabel('Probability')
xlabel('Distance (m)')
title('FPD CDF for Stair-Step Path')


