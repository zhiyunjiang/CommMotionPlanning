function [alpha, beta, rho] = SHMP_par_estimation(x_sample, y_sample, sig_SHMP, res, per, thr_sample_no, thr_SNR)
% function [alpha, beta, rho] = SHMP_par_estimation(x_sample, y_sample, sig_SHMP, res, per, thr_sample_no, thr_SNR)
% 
% this function uses the channel samples to extract the channel parameters
%
% input:
% - x_sample and y_sample are vectors of coordinations of samples
% - sig_SHMP is the vector of shadowing and MP components (by substracting 
% PL component from the original signal)
% - res is between two samples (meter)
% - per is percentage of measurements
% - thr_sample_no is the minimum number of pairs of samples for specific distances. We do
% not consider the statistics is accurate if the pairs of samples are
% smaller than thr_sample_no.
% - thr_SNR is the threshold to accept/reject a fit. If the difference between
% fit curve and numerical ACF is large than thr_SNR, then we
% reject this fit, and consider the channel as uncorrelated.
% 
% output:
% - alpha: shadowing power (with a square, i.e. alpha in Mehrzad's thesis)
% - beta: correlation distance of shadowing
% - rho: power of MP (with a square, i.e. sigma^2 in Mehrzad's thesis)

% This code is editted by Yuan Yan based on Mehrzad's original
% WLS_decorr_2D_SHSS code

global weight_filter

% ave_power_SHMP is the average power of the shadowing and MP components
ave_power_SHMP = mean(sig_SHMP.^2);

M = length(x_sample);
x_min = min(x_sample);
x_max = max(x_sample);
y_min = min(y_sample);
y_max = max(y_sample);

N = ceil(norm([x_max-x_min y_max-y_min])/res);

count_ACF = zeros(N + 1,1);
numerical_ACF_raw = zeros(N + 1,1);

for i = 1:M
    % distance between [x(i),y(i)] to all the points in x and y
    d = sqrt((x_sample - x_sample(i)).^2 + (y_sample - y_sample(i)).^2);
    % change distance to index
    d = round(d/res);
    % used to calculate autocorrelation numerically
    numerical_ACF_raw(d+1) = numerical_ACF_raw(d+1) + sig_SHMP(i)*sig_SHMP(:);
    count_ACF(d+1) = count_ACF(d+1) + 1;
end

% remove all the entries in count_ACF that are smaller than thr_sample_no.
% This removes all the points that the number of pairs of samples for
% specific distances are smaller than thr_sample_no for accuracy. 
ind_temp1 = find(count_ACF >= thr_sample_no);
% remove all the points that the autocorrelation is larger than the
% autocorrelation at 0 point, because autocorrelation  suppose to be
% maximized at 0 point.
ind_temp2 = find(numerical_ACF_raw(1)./numerical_ACF_raw >= 1);

% calculate autocorrelation numerically
ind_valid = intersect(ind_temp1, ind_temp2);
numerical_ACF = numerical_ACF_raw(ind_valid)./count_ACF(ind_valid);

if isempty(numerical_ACF)
    alpha = ave_power_SHMP;
    beta = eps;
    rho = 0;
    return;
end

d_ACF = (ind_valid(:) - 1)*res;
% weight is weighing vector, per is the percentage of samples. Note that the
% idea here is that, if per is small, there are not that many
% samples, so the numerical estimation of ACF for large distance
% pairs is not that accurate. As a result, we put less weight on
% those points. The design of this function is just based on this
% simple intuition, there is no mathematical reason behind it.
weight = exp(-weight_filter*(d_ACF - d_ACF(1))/per);

% myfunc defines objective function to find the best fit for the power and 
% decorrelation distance of shadowing. Since at point numerical_ACF(1), there
% suppose to contain the power of MP, so we fit from numerical_ACF(2) to the
% rest
myfunc=@(x) sum(weight(2:end).*(numerical_ACF(2:end) - x(1)*exp(-d_ACF(2:end)/x(2))).^2);
% out = patternsearch(myfunc,[ave_power_SHMP 10]);
% out = fmincon(myfunc,[ave_power_SHMP-1 20],[-1 0;1 0;0 -1;0 1],[0;50;0;10]);
out = fmincon(myfunc,[ave_power_SHMP-1 10],[-1 0;1 0;0 -1;0 1],[0;50;0;5]);

%        out = search_ACF(myfunc);
% find the corresponding parameters
alpha = out(1);
beta = out(2);
% rho = numerical_ACF(1) - alpha
rho = ave_power_SHMP - alpha;
if rho < 0
    rho = 0.01;
end

% the estimated ACF
est_ACF = alpha*exp(-d_ACF/beta) + rho*zeros(length(d_ACF),1);
% calculate how far is the estimated fit from the numerical ACF, if
% they are too far (>thr_SNR), the reject this fit and return the
% channel as uncorrelated
est_SNR = S2N(weight.*numerical_ACF, weight.*est_ACF);
% if est_SNR < thr_SNR || alpha < 0 || beta < 0 || rho < 0
%     alpha = ave_power_SHMP;
%     beta = eps;
%     rho = 0;
% end

return;