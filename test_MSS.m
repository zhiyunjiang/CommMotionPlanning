%%
%Test MidSetSolver
%Deterministic Comm Path Planning
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

maze = [RectObs([-1, -1.5, 3, -5.5]), RectObs([0, -0.5, 6, -2]), ...
            RectObs([1, 0.5, 3, -5.5]), RectObs([2, 1.5, 6, -2])];
simple = [CircObs(1.5, [0,4]), RectObs([4, 2, 0, -2 ]), CircObs(1, [2,2])];

dense = RectObs.empty;
dx = 0.3;
dy = 0.3;
for i=1:floor(x_max - x_min)
    
    for j = 1:floor(y_max - y_min)
        y_jitter = 2*dy*rand(1); 
        x_jitter = 2*dx*rand(1); 
        obs = RectObs([i+x_min+dx + x_jitter, i+x_min + x_jitter, j+y_min+dy + y_jitter, j+y_min+y_jitter]);
        dense((i-1)*floor(y_max - y_min) + j) = obs;
    end 
end

obs_mod = ObstacleMod(dense);
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

[gamma_TOT_dB, gamma_PL_SH_dB, gamma_PL_dB, gamma_SH_dB, gamma_MP_LIN, g_x, g_y] = channel_simulator(region, ...
    q_b, K_PL, n_PL, alpha, beta, N_sin, PSD_at_f_c, lambda, K_ric, res, corr_mp);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find Set of Connected Points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Channel power connectiivty threshold (dB)
gamma_TH = -5;
sz = size(gamma_TOT_dB);
connected = [];
for i = 1:sz(1)
    y_loc = (i-1)/res + y_min;
    for j = 1:sz(2)
        if gamma_TOT_dB(i, j) >= gamma_TH
            x_loc = (j-1)/res + x_min;
            connected = [connected; [x_loc, y_loc]];
        end
    end
end

% scatter(connected(:,1), connected(:,2)); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Path Plan to Connected Spots Using Manhattan distance cost function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cost function
cost_fxn = @(n1, n2, path) ManhattanDistance(n1, n2, path);
path_res = 10;
problem_instance = PathPlanningProblem(region, path_res, source, connected, obs_mod, cost_fxn);
problem_instance.setSource2(dest);
%%
figure(1)
clf
hold on
problem_instance.plotProb(1)
%%
%stop on first solution, stop after 2 seconds, unless a solution hasn't
%been found
solution_not_required = 0;
max_iterations = Inf;
%max run time in seconds
max_run_time = 60;
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);

%1 - uniform random (RRT), 2 - deterministic sequence (RDT)
type = 1;
dest_freq = 10;
sequence = [];%sampler.sequence;%use previous run's sequence
sampler = Sampler(type, dest_freq, sequence);

do_rewire = 1;
steer_rad = 5;

%theta in [0, 1] - how much weight to put on total distance (as opposed to
%distance to connectivity)
%theta = 0 - minimize distance to connectivity ignore distance from
%connectivity to goal
%theta = 1 - minimize total distance traveled

theta = 1;

mss = MidSetSolver(sampler, stop_criteria, do_rewire, steer_rad);
tic
mss.solve(problem_instance, theta);
toc
%%
best_bridge = mss.getBBSF();
figure(1)
hold on
best_bridge.pathToRoot(1);
legend('Path from Source', 'Path to Destination');
title('Feasible Path Found via Two Rapidly Exploring Random Trees');

figure(2)
series = mss.getHistBBSFCost();
plot(series(:,1), series(:,2))
xlabel('Elapsed Time');
ylabel('Minimum Cost so Far');
title('Minimum Cost over Time');
%%
thetas = 0:0.1:1;
n = 100;
sums = zeros(size(thetas));

%stop on first solution, stop after 2 seconds, unless a solution hasn't
%been found
solution_not_required = 0;
max_iterations = Inf;
%max run time in seconds
max_run_time = 20;
stop_criteria = StopCriteria(solution_not_required, max_iterations, max_run_time);

for i = 1:11
    theta = thetas(i);
    for j = 1:n
        mss = MidSetSolver(sampler, stop_criteria, do_rewire, steer_rad);
        mss.solve(problem_instance, theta);
        best_bridge = mss.getBBSF();
        sums(i) = sums(i) + best_bridge.costFromSource();
    end
    fprintf('Completed theta set\n');
end

averages = sums/n;

plot(thetas, averages);

