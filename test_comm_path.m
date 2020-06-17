%Will combine elements of test_channel and RDT_test to generate a path in
%through the area modeled by the channel simulator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_max = 6;
x_min = -4;
y_max = 5.5;
y_min = -5;
region = [x_max x_min y_max y_min];

source = [0, 5];
dest = [0, -5];

% Position of the base station (remote station or transmitter)
q_b = [0 0];

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
star = 1;
sequence = [];
%use a lower resolution for path planning, else problem becomes unwieldy
%plues, don't need same resolution
rdt_res = 10;%tick marks/meter
%obstacles = [RectObs([-1, -1.5, 3, -5]), RectObs([0, -0.5, 5.5, -2]), ...
            % RectObs([1, 0.5, 3, -5]), RectObs([2, 1.5, 5.5, -2])];
obstacles = []; 
obstacle_mod = ObstacleMod(obstacles);
dest_freq = 10;

fprintf('Finding Path with RDT ...\n')
RDT = PathPlanningProblem(region, rdt_res, obstacles);
[V, goal_vertex, ~] = RDT.FindPath(source, dest, star, sequence, dest_freq);
fprintf('RDT path found!\n')

path = goal_vertex.pathToRoot(0);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate the path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Path distance can be pulled directly from the leaf node
C_dist = goal_vertex.distToHere;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Represent Connectivity at each node as a random 
% variable ~N(path loss power dB, shadowing power)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%mean of the RV
gamma_PL_dB = generate_pathloss(region, q_b, res, K_PL, n_PL);

%shadowing power
alpha = 10;
%standard deviation of shadowing and power conn
sigma_SH = sqrt(alpha);

%Channel power connectiivty threshold (dB)
gamma_TH = -20;

N_sin = 5000;                    
PSD_at_f_c = 30;
lambda = 0.125;
K_ric = 10;     

% corr_mp = 1 -> correlated multipath 
% corr_mp = 0 -> uncorrelated multipath 
corr_mp = 0;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predict FPD based on gamma_PL and sigma_SH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Communication base cost will be distance traveled until the threshold
%power (dB) is met, i.e. distance traveled until connectivity. 

%However, there may not be a point along our path where we're gaurantted to
%be connected, so we calculate the distance at which the probability of
%connectivity is above a threshold value.

fprintf('Calculating Prior FPD CDF\n');
%[PMF,CDF] = pathPMF(path, 1/rdt_res, gamma_TH, x_min, y_min, res, gamma_PL_dB, sigma_SH, beta);

[PMF1, CDF1] = pathPMFComplete(path, 1/rdt_res, beta, sigma_SH, K_ric, gamma_TH, gamma_PL_dB, x_min, y_min, res );
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate several the channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now generate the rest of the channel. For detailed explanations of parameters and
%reference, see ChannelSim/test_channel.m


%Generate several channels and see how the predicted distribution compares
%to the simulated
n_sim = 30;

[g_x, g_y] = meshgrid(x_min:1/res:x_max,y_min:1/res:y_max);
[M, N] = size(g_x);

gamma_SH_dBs = [];
gamma_MP_LINs = [];

for i=1:n_sim
    gamma_SH_dBs(:,:,i) = generate_shadowing(alpha, beta, N_sin, PSD_at_f_c, g_x, g_y, res);
    gamma_MP_LINs(:,:,i) = generate_multipath(lambda, g_x, g_y, res, K_ric, corr_mp);
    
    if mod(i,10)==0
       fprintf('%i iterations complete...\n',i); 
    end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute Empirical FPD CPF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Calculating Simulated FPD CDF using %d simulations\n', n_sim);
sim_locations = zeros(n_sim, 2);
sim_costs = zeros(n_sim, 1);
for i = 1:n_sim
    gamma_SH_dB = gamma_SH_dBs(:,:,i);
    gamma_MP_LIN = gamma_MP_LINs(:,:,i);
    
    gamma_PL_SH_dB = gamma_PL_dB + gamma_SH_dB;
    gamma_PL_SH_LIN = 10.^(gamma_PL_SH_dB / 10);
    gamma_TOT_LIN = gamma_PL_SH_LIN.*gamma_MP_LIN;
    gamma_TOT_dB = 10*log10(gamma_TOT_LIN);
  
    %compare to the channel with shadowing and multipath effects
    sim_costs(i) = distToConnectivity(path, 1/rdt_res, gamma_TH, x_min, y_min, res, gamma_TOT_dB);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the Channel and path accross the channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pretty plot the last channel generated. Code taken from ChannelSim/test_channel.m
figure(1);
clf

fnt_siz = 25;

surf(g_x, g_y, gamma_TOT_dB, 'EdgeColor','none');
% light
% shading interp
xlabel('x (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
ylabel('y (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
zlabel('Received power (PL + SH + MP) (dBm)','FontSize', fnt_siz ,  'FontWeight','bold');
axis tight
grid on

set(gca, 'FontSize', fnt_siz, 'FontWeight', 'bold');

%maximize(f)
hold on

path = goal_vertex.pathToRoot(1);
hold on
for i=1:length(obstacles)
    obstacles(i).plotObstacle();
    hold on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optionally plot the entire tree generated by RRT*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%f = figure();
%V(1).plotTree('r');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now compare the predicted FPD to that found via Monte Carlo methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot histogram of distances
figure(2);
clf

histogram(sim_costs,'Normalization', 'cdf', 'BinWidth', 1/rdt_res);
hold on
plot(CDF(:,1), CDF(:,2));
plot(CDF1(:,1), CDF1(:,2));
xlabel('Distance Along Path (m)');
ylabel('Cumulative Probability of Connectivity');
%xlim([4,10])
legend('Simulated FPD','Predicted FPD');

%%
function [PMF,CDF] = pathPMF(path, rdt_step, gamma_th, x_min, y_min, res, gamma_PL_dB, sigma, beta)
    PMF = zeros(size(path));
    CDF = zeros(size(path));
    PL = 0;
    p_no_conn_to_here = 1;
    CP = 0;
    for i = 1:length(path)
       %find probability that FPD is here
       vertex = path(i,:);
        if i == 1
            PL = getPL(vertex, x_min, y_min, res, gamma_PL_dB);
            p_no_conn_here = normcdf(gamma_th, PL, sigma);
        else
            PL_prev = PL;

            PL = getPL(vertex, x_min, y_min, res, gamma_PL_dB);
            integrand = @(g) normcdf(gamma_th, PL+exp(-rdt_step/beta)*(g - PL_prev),...
                sigma*sqrt((1 - exp(-2*rdt_step/beta)))) .* normpdf(g, PL_prev, sigma);
            %probability of no connection here given no connection
            %previously
            p_no_conn_here = integral(integrand, -Inf, gamma_th)*normcdf(gamma_th, PL_prev, sigma);
        end
        p_fpd_here = (1-p_no_conn_here)*p_no_conn_to_here;
        CP = CP + p_fpd_here;
        PMF(i,:) = [(i-1)*rdt_step, p_fpd_here];
        CDF(i,:) = [(i-1)*rdt_step, CP];
        p_no_conn_to_here = p_no_conn_to_here*p_no_conn_here;
   end
end


function [PMF,CDF] = pathPMFComplete(path, rdt_step, beta, sigma_SH, K_ric, gamma_TH, gamma_PL_dB, x_min, y_min, res)
    PMF = zeros(size(path));
    CDF = zeros(size(path));
    
    
    ro = exp(-rdt_step/beta);
    sigma_SH_c = sigma_SH*sqrt(1-ro);
    
    nu_MP = sqrt(K_ric/(K_ric + 1));%non-centrality parameter
    sigma_MP = 1/sqrt(2*(1+K_ric));%spread parameter
    ric_MP = makedist('Rician','s',nu_MP,'sigma',sigma_MP);
    
    gamma_PL_0 = getPL(path(1,:), x_min, y_min, res, gamma_PL_dB);
    
    %use linear (not dB) power when working with MP
    J = @(w) normpdf(w, 0, sigma_SH).*cdf(ric_MP, toLin(gamma_TH - gamma_PL_0 - w));
    
    %denominator to be used in a number of calculations
    p_Y0_lessthan_TH = integral(J, -Inf, Inf);
    % p_no_prior_connection given not connection at start position
    p_no_prior_connection = 1;
    
    for i = 2:length(path)
        vertex = path(i,:);
        %Calculate the probability of first passage distance at this vertex
        %These calculations assume straighline paths
        %Using the fomrulation given in "Statistics of the Distance 
        %Traveled until Connectviity for Unmanned Vehicles",
        %Muralidharan and Mostofi, 2018
        
        J_prev = J;
        
        gamma_PL = getPL(vertex, x_min, y_min, res, gamma_PL_dB); 
        
        J = @(w) (cdf(ric_MP, toLin(gamma_TH - gamma_PL - w)) / ro) .* ...
            arrayfun( @(wk) integral(@(s) normpdf(s, wk,sigma_SH_c) .* J_prev(s/ro), -Inf , Inf), w );
        
        p_no_connection_to_here = integral(J, -Inf, Inf)/p_Y0_lessthan_TH;
        i
        p_FPD_here = p_no_prior_connection - p_no_connection_to_here
        
        PMF(i,:) = [(i-1)*rdt_step, p_FPD_here];
        CDF(i,:) = [(i-1)*rdt_step, CDF(i-1,2) + p_FPD_here ];
        
        p_no_prior_connection = p_no_connection_to_here;
    end
    
end

function lin = toLin(val_db)
    lin = 10.^(val_db/10);
end

function pl = getPL(vertex, x_min, y_min, res, gamma_PL_dB)
   x_index = round((vertex(1) - x_min)*res) + 1;
   y_index = round((vertex(2) - y_min)*res) + 1;
   pl = gamma_PL_dB(y_index, x_index); 
end

function cost = distToConnectivity(path, rdt_step, gamma_th, x_min, y_min, res, gamma_TOT_dB)
    cost = Inf;
    for i = 1:length(path)
        node = path(i,:);
        x_index = round((node(1) - x_min)*res) + 1;
        y_index = round((node(2) - y_min)*res) + 1;
        if gamma_TOT_dB(y_index, x_index)>=gamma_th
            cost = (i-1)*rdt_step;
            break 
        end
    end
end

