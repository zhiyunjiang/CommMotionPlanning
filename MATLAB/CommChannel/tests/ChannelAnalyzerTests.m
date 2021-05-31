%%
%Setup
q_b = [-450 0];
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
pathy = zeros([200, 1]);
pathx = cumsum(ones([200, 1]))/2;

path = [flip(pathx), pathy];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test Arjun's recursive approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Test Performance of Recursive Markovian Approximation of Expected FPD\n');
exp4 = @() ca.ExpectedFPD(path,4); 
exp4_runtime = timeit(exp4);
fprintf('Expected Runtime for a 200-point path is %5f\n', exp4_runtime);
%%
%verify we accurately calculate straighline-path expected FPD
[approx_PMF_1, distances_1] = ca.ApproxFPDPMF(path,1);
[approx_PMF_2, distances_2] = ca.ApproxFPDPMF(path,2);
fprintf('Found PMFs\n');
approx_CDF_1 = cumsum(approx_PMF_1);
approx_CDF_2 = cumsum(approx_PMF_2);
%%
figure()
clf;
hold on
plot(cumsum(distances_1), approx_CDF_1);
plot(cumsum(distances_2), approx_CDF_2);
%verify by simulating several channels and calculating expected distance
%to connectivity

all_dists = ca.simulateFPD(path, 100, 0);
histogram(all_dists, 'Normalization', 'cdf', 'BinEdges',(cumsum(distances_1)));
legend( 'Arjun full', 'Arjun no mp' , 'sims');

ylabel('Probability')
xlabel('Distance (m)')
title('FPD CDF for Straight Line Path')

%%
%Testing with non-straighline paths

root = [100, 1];
deltax = ones([99, 1]);deltax(2:2:end) = 0;
deltay = zeros([99, 1]);deltay(2:2:end) = 1;
pathx = 100 - cumsum(deltax);
pathy = 1 + cumsum(deltay);
path = [root;pathx, pathy];

[approx_PMF_1, distances_1] = ca.ApproxFPDPMF(path,1);
[approx_PMF_2, distances_2] = ca.ApproxFPDPMF(path,2);
fprintf('Found PMFs\n');
approx_CDF_1 = cumsum(approx_PMF_1);
approx_CDF_2 = cumsum(approx_PMF_2);

figure()
clf;
hold on
plot(cumsum(distances_1), approx_CDF_1);
plot(cumsum(distances_2), approx_CDF_2);

%verify by simulating several channels and calculating expected distance
%to connectivity

all_dists = ca.simulateFPD(path, 20, 0);
histogram(all_dists, 'Normalization', 'cdf', 'BinEdges',(cumsum(distances_1)));
legend('Markovian Approximation', 'sims');

ylabel('Probability')
xlabel('Distance (m)')
title('FPD CDF for Stair-Step Path')


