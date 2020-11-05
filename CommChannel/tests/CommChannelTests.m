%CommChannel Tests and Examples
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel 1 - MP Correlated, Rician
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q_b = [0 0];
n_PL = 4.2;  
K_PL = -40;

sh_decorr = 6;    
alpha = 8.41;
sigma_SH = sqrt(alpha);
PSD_at_f_c = 30;

lambda = 0.125;
K_ric = 1.59;
mp_decorr = 0.4*lambda;         
corr_mp = 1;
sigma_mp = 1;

cp = ChannelParams(q_b, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, ...
                    corr_mp, PSD_at_f_c, sigma_mp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Comm Channel Simulation & Generation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%keep it small for testing
x_max = 13.56;
x_min = -2.5;
y_max = 1.5;
y_min = -15.36;
region = [x_max x_min y_max y_min];
%simulation parameter
N_sin = 5000;

res = 2/mp_decorr;
cc = CommChannel(cp, N_sin, region, res);
cc.simulateSH();cc.simulateMP();
[pos, vals] = cc.randSample(5);

gammas = cc.getGammaTOTdBAtPt(pos);

error = vals - gammas;
if sum(error) ~= 0
   warning('error should be zero');
end

%test plotting functions
cc.plot(0);
figure();
cc.plotConnected(-65);
figure()
receiver_noise = 1e-10; R = 2; BER = 1e-6;
qos = QoSParams(BER, R, receiver_noise);
cc.plotRequiredTXPower(qos)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel 2 - MP Uncorrelated, Rician
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Position of the base station (remote station or transmitter)
q_b = [2, 3.4];         

% corr_mp = 0 -> uncorrelated multipath; corr_mp = 1 -> correlated multipath 
corr_mp = 0;

cp = ChannelParams(q_b, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, ...
                    corr_mp, PSD_at_f_c, sigma_mp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Comm Channel Simulation & Generation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%keep it small for testing
x_max = 13.56;
x_min = -2.5;
y_max = 1.5;
y_min = -15.36;
region = [x_max x_min y_max y_min];
%simulation parameter
N_sin = 5000;

res = 2/mp_decorr;
cc = CommChannel(cp, N_sin, region, res);
cc.simulateSH();cc.simulateMP();
[pos, vals] = cc.randSample(5);

gammas = cc.getGammaTOTdBAtPt(pos);

error = vals - gammas;
if sum(error) ~= 0
   warning('error should be zero');
end

%test plotting functions
cc.plot(0);
figure();
cc.plotConnected(-75);
figure()
receiver_noise = 1e-10; R = 2; BER = 1e-6;
qos = QoSParams(BER, R, receiver_noise);
cc.plotRequiredTXPower(qos)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Channel 3 - MP Uncorrelated, lognormal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%keep it small for testing
x_max = 13.56;
x_min = -2.5;
y_max = 1.5;
y_min = -15.36;
region = [x_max x_min y_max y_min];
%simulation parameter
N_sin = 5000;

res = 2/mp_decorr;
cc = CommChannel(cp, N_sin, region, res);
cc.simulateSH();cc.simulateMP(2);
[pos, vals] = cc.randSample(5);

gammas = cc.getGammaTOTdBAtPt(pos);

error = vals - gammas;
if sum(error) ~= 0
   warning('error should be zero');
end

%test plotting functions
cc.plot(0);
figure();
cc.plotConnected(-65);
figure()
receiver_noise = 1e-10; R = 2; BER = 1e-6;
qos = QoSParams(BER, R, receiver_noise);
cc.plotRequiredTXPower(qos)