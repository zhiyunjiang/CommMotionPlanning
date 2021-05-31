%demonstrate how to use solver
%%
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


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Comm Channel Simulation & Generation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_max = 25;
x_min = 0;
y_max = 25;
y_min = 0;
region = [x_max x_min y_max y_min];
%simulation parameter
N_sin = 5000;

%res = 2/mp_decorr;
res = 1;
cc = CommChannel(cp, N_sin, region, res);
%cc.simulateShadowing();
cc.plot()
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the Channel Analyzer, Preview Prior Probability of Connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%minimum power required for connectivity
gamma_TH = -84.5;
no_mp = 1;
ca = ChannelAnalyzer(cc, gamma_TH, no_mp);
% figure(2)
% clf
% ca.plot_priors()

%%
%mininum probability of channel power being greater than or equal to
%gamma_TH require to count spot as a connected spot
p_TH = 0.95;
% sz = size(cc.getGammaPLdB());
% connected = [];
% for i = 1:sz(1)
%     y_loc = (i-1)/res + y_min;
%     for j = 1:sz(2)
%         if p_TH <= 1-ca.NoConnectionPrior([i,j])
%             x_loc = (j-1)/res + x_min;
%             connected = [connected; [x_loc, y_loc]];
%         end
%     end
% end
connected = [1,1];
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path Plan to Connected Spots Using Expected distance cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%setup obstacles
source = [24, 24];

maze = [RectObs([-1, -1.5, 3, -5.5]), RectObs([0, -0.5, 6, -2]), ...
            RectObs([1, 0.5, 3, -5.5]), RectObs([2, 1.5, 6, -2])];
simple = [CircObs(2, [5,5]), CircObs(2, [5, 10]),...
    CircObs(2, [5,15]), CircObs(2, [5, 20]),...
    CircObs(2, [15,7.5]), CircObs(2, [15, 12.5]),...
    CircObs(2, [15,17.5]), CircObs(2, [15, 22.5]),...
    CircObs(2, [10, 7.5]), CircObs(2, [10, 12.5]),...
    CircObs(2, [10,17.5]), CircObs(2, [10, 2.5])];
obs_mod = ObstacleMod(simple);
%cost function
cost_fxn = @(n1, n2, path) ExpectedDistance(n1, n2, ca, 2);
path_res = res;
problem_instance = PathPlanningProblem(region, path_res, source, connected, obs_mod, cost_fxn);

figure (3)
clf
hold on
problem_instance.plotProb(1)
%%
pgs = PotentialGameSolver(problem_instance, ca);

path = pgs.solve(100*(25^2));
hold on
plot(path(:,1), path(:,2));