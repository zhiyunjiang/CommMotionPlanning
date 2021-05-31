%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup our robot with reasonable stats. 
% We assume wifi-like ranges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency - typical wifi frequencies are 2.4, 5, & 6 GHz
F5GHz = 5;
F2_4GHz = 2.4;
% 1 Watt seems to be close to an absolute maximum, 100-200 mW (20-23 dBm)
% more reasonable
MOBILE_TX_POWER = TXPwr(23);
MOBILE_TX_POWER.view();

%QoS related params
receiver_noise = 1e-10;  BER = 1e-6;
R = 6;% MQAM, M = 64;
qos = QoSParams(BER, R, receiver_noise);
% Calculate Channel Gain Threshold
gamma_TH = qos.thresholdChannelGain(MOBILE_TX_POWER.W);
% Set chance constraint
%p_th = 0.94;
p_th = 0.9;


[n_PL, K_PL, sh_decorr, sigma_SH, PSD_at_f_c, mp_decorr, lambda, K_ric,...
    corr_mp, sigma_mp, region, N_sin, res, n_samples, n_sims, ~] = getConstParams(); 

%Setup first pair of points
%qb1 = [4,3];
qb1 = [25,26];
cp1 = ChannelParams(qb1, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);
 
%qb2 = [23, 28];
qb2 = [2, 2];
cp2 = ChannelParams(qb2, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);

%Multipath model: 1 - rician, 2 - log normal
MP_Mdl = 2;%use log normal model for simulating multipath

cc1 = CommChannel(cp1, N_sin, region, res);
cc1.simulateSH(); cc1.simulateMP(MP_Mdl);

cc2 = CommChannel(cp2, N_sin, region, res);
cc2.simulateSH(); cc2.simulateMP(MP_Mdl);

%Sample
[obs_pos, obs_vals] = cc1.randSample(n_samples);
[obs_pos_poi, obs_vals_poi] = cc2.randSample(n_samples);


%Create Predicted Channel
pc1 = PredictedChannel(cp1, cc1, obs_pos, obs_vals);
pc2 = PredictedChannel(cp2, cc2, obs_pos_poi, obs_vals_poi);


%Setup second pair of points
%qb3 = [45,37];
qb3 = [26,24];
cp3 = ChannelParams(qb3, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);
 
%qb4 = [32, 17];
qb4 = [48, 2];
cp4 = ChannelParams(qb4, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);

%Multipath model: 1 - rician, 2 - log normal
MP_Mdl = 2;%use log normal model for simulating multipath

cc3 = CommChannel(cp3, N_sin, region, res);
cc3.simulateSH(); cc3.simulateMP(MP_Mdl);

cc4 = CommChannel(cp4, N_sin, region, res);
cc4.simulateSH(); cc4.simulateMP(MP_Mdl);

%Sample
[obs_pos, obs_vals] = cc3.randSample(n_samples);
[obs_pos_poi, obs_vals_poi] = cc4.randSample(n_samples);


%Create Predicted Channel
pc3 = PredictedChannel(cp3, cc3, obs_pos, obs_vals);
pc4 = PredictedChannel(cp4, cc4, obs_pos_poi, obs_vals_poi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup Third pair of points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%qb5 = [10,45];
qb5 = [24,24];
cp5 = ChannelParams(qb5, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);
 
%qb6 = [28, 28];
qb6 = [24, 49];
cp6 = ChannelParams(qb6, n_PL, K_PL, sigma_SH, sh_decorr, lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp);

%Multipath model: 1 - rician, 2 - log normal
MP_Mdl = 2;%use log normal model for simulating multipath

cc5 = CommChannel(cp5, N_sin, region, res);
cc5.simulateSH(); cc5.simulateMP(MP_Mdl);

cc6 = CommChannel(cp6, N_sin, region, res);
cc6.simulateSH(); cc6.simulateMP(MP_Mdl);

%Sample
[obs_pos, obs_vals] = cc5.randSample(n_samples);
[obs_pos_poi, obs_vals_poi] = cc6.randSample(n_samples);


%Create Predicted Channel
pc5 = PredictedChannel(cp5, cc5, obs_pos, obs_vals);
pc6 = PredictedChannel(cp6, cc6, obs_pos_poi, obs_vals_poi);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot region where both have channel gain < threshold

conn_field_1 = pc1.getConnectionField(p_th, gamma_TH).*pc2.getConnectionField(p_th, gamma_TH);

conn_field_2 = 2*pc3.getConnectionField(p_th, gamma_TH).*pc4.getConnectionField(p_th, gamma_TH);
conn_field_3 = 4*pc5.getConnectionField(p_th, gamma_TH).*pc6.getConnectionField(p_th, gamma_TH);


conn_field_total = conn_field_1 + conn_field_2 + conn_field_3;
cc1.plotField(conn_field_total, 'Region of Joint Connectivity', 'Connectivity Indicator')
hold on
%might as well plot the base stations, too
plotRS(qb1,[0.35, 0.55, 0.91]); plotRS(qb2,[0.35, 0.55, 0.91]);
plotRS(qb3,'c'); plotRS(qb4,'c');
plotRS(qb5,'y'); plotRS(qb6,'y');

%%
%Vornoi Diagrams
% figure
% [x1,y1] = find(conn_field_1);
% voronoi(x1,y1)
% hold on
% [x2,y2] = find(conn_field_2);
% voronoi(x2,y2)
% hold on
% [x3,y3] = find(conn_field_3);
% voronoi(x3,y3)
% hold off
% 
% %set colors
% myfig = gcf;
% axesObj = myfig.Children;
% lineObjs = axesObj.Children;
% line3 = lineObjs(3);
% line3.Color='r';
% line4 = lineObjs(4);
% line4.Color='r';
% line5 = lineObjs(5);
% line5.Color='g';
% line6 = lineObjs(6);
% line6.Color='g';