%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PredictedChannel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a predicted model of the channel given observations
% Can either estimate channel parameters or use the true channel parameters
% given.

classdef PredictedChannel < handle
    
    properties (Access = public)
        cp;
        cc;
        
        obsCount;
        obsVals;
        obsPos;% in gridCoordinates
        obsPL;
        obsDiff;
        obsCov;
        UInv;
        
         means;
         vars;
         
         kPl_est;
         nPl_est;
         alpha_est;%estimated shadowing power (variance)
         beta_est;%shadowing decorrelation distance
         scaledBeta_est
         rho_est;%estimated MP power (variance). Assumes log normal distirbution of multipath
    end
    
    properties (Access = private)
       useEstimatedParams = 0; 
    end
    
    methods (Access = public)
        function this = PredictedChannel(channel_params, comm_channel, observation_pos,...
                                        channel_samples, use_estimates)
                                    
            if nargin < 5%by default, use the true parameters
                use_estimates = 0;
            end
            
            this.useEstimatedParams = use_estimates;
            this.cp = channel_params;
            this.cc = comm_channel;
            
            this.obsCount = length(channel_samples);
            this.obsPos = observation_pos;
            this.obsVals = channel_samples;
            if this.useEstimatedParams
                this.estimatePLParams();
            end
            this.obsPL = this.estimatePL(this.obsPos);
            this.obsDiff = this.obsVals - this.obsPL;

            this.obsCov = this.calcObsCov();
            this.UInv = (this.obsCov + this.rho()*eye(this.obsCount))^-1;
            
            
            this.setStats();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % estimatePL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Estimate the path loss gain at a point in grid coordinates
        % Input:
        % this - reference to the PredictedChannel object
        % pt - either a single point [x,y] or matrix of points 
        %       [x1,y1;...xn,yn]. Assumes points are in raw coordinates
        % Output:
        % gamma_PL_dB_at_point - returns the path loss gain in dB
        function gamma_PL_dB_at_point = estimatePL(this, pt)
            if this.useEstimatedParams
                gamma_PL_dB_at_point = this.kPl_est - 10*this.nPl_est*...
                            log10(sqrt(sum((this.cp.qBase - pt).^2, 2)));
            else
                gamma_PL_dB_at_point = this.cc.getGammaPLdBAtPt(pt);
            end
           
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % getMeanAtGridPoint
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get the mean (expected value) of the normally distributed channel
        % gain dB random variable.
        % Input:
        % this - reference to the PredictedChannel object
        % pt - a single point [x,y]
        %
        % Output:
        % mean_dB - returns expected channel gain dB
        function mean_dB = getMeanAtGridPoint(this, pt)
           grid_pt = this.cc.toGridCoordinate(pt);
           mean_dB = this.means(grid_pt(1), grid_pt(2)); 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % posteriorPConn
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculates the probability that the channel gain is above a giveb
        % channel gain threshold in dB.
        % Input:
        % this - reference to the PredictedChannel object
        % pt - a single point [x,y]
        % gamma_th - the minimum required channel gain for communication
        %
        % Output:
        % p_conn - the probability that the channel gain at pt is above
        %           gamma_th
        function p_conn = posteriorPConn(this, pt, gamma_th)
            grid_pt = this.cc.toGridCoordinate(pt);
            mean = this.means(grid_pt(1), grid_pt(2));
            variance = this.vars(grid_pt(1), grid_pt(2));
            
            p_conn = normcdf(gamma_th, mean, sqrt(variance), 'upper');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotPosteriors2D
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots the posterior probability of connectivity over the
        % workspace given a minimum required channel gain.
        % Input:
        % this - reference to the PredictedChannel object
        % gamma_th - the minimum required channel gain for communication
        function plotPosteriors2D(this, gamma_th)
           [gx, gy] = this.cc.getMeshGrid();
           grd_sz = size(gx);
           x_vals = gx(1,:);
           y_vals = gy(:,1);
           post_map = zeros(grd_sz);
           %imagesc plots the first index as the column number, second as
           %row number, so we save off data as a transposed array
           for i = 1:grd_sz(2)
              for j = 1:grd_sz(1)
                    p_conn = this.posteriorPConn([x_vals(i), y_vals(j)], gamma_th);
                    post_map(j,i) = p_conn;
              end
           end
           imagesc(x_vals, y_vals, post_map);
           set(gca, 'YDir', 'normal');
           c = colorbar;
           c.Label.String = 'Probability of Connectivity';
           xlabel('x (m)');
           ylabel('y (m)');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % getConnectionField
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fetches a matrix representing a binary field with 1 indicating
        % communication is possible and 0 indicating connection is not
        % possible. This is based on a minimum required channel gain, which
        % in turn can be calculated from other requirements (e.g. SNR, BER,
        % etc...)
        % Input:
        % this - reference to the PredictedChannel object
        % p_th - probability threshold. Since we're working with a
        %           predicted channel, this is the probability that the
        %           channel is good enough for communication
        % gamma_th - the minimum required channel gain for communication
        %
        % Output:
        % conn_field - matrix of the same dimension as the grid-scaled
        %               workspace. Represents the binary connection field.
        function conn_field = getConnectionField(this, p_th, gamma_th)
           [gx, gy] = this.cc.getMeshGrid();
           x_vals = gx(1,:);
           y_vals = gy(:,1);
           conn_field = zeros(flip(size(gx)));
           for i = 1:length(x_vals)
              for j = 1:length(y_vals)
                    conn = this.posteriorPConn([x_vals(i), y_vals(j)], gamma_th) >= p_th;
                    conn_field(i,j) = conn;
              end
           end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotConnected2D
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plots the binary connectivity field over the workspace given a
        % minimum required channel gain and threshold probability.
        % Input:
        % this - reference to the PredictedChannel object
        % p_th - probability threshold. Since we're working with a
        %           predicted channel, this is the probability that the
        %           channel is good enough for communication
        % gamma_th - the minimum required channel gain for communication
        function plotConnected2D(this, p_th, gamma_th)
           conn_field = this.getConnectionField(p_th, gamma_th);
           title = strcat('Connected Regions for \Gamma_{th} = ', sprintf('%d dB', gamma_th),...
                           ', p_{conn} >= ', sprintf('%g',p_th));
           this.cc.plotField(conn_field,title);
        end
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % setStats
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % For each point on the grid, calculates the mean and variance of
        % the normal RV's that model channel gain dB.
        % Input:
        % this - reference to the PredictedChannel object
        function setStats(this)
           [gx, gy] = this.cc.getMeshGrid();
           grd_sz = size(gx);
           x_vals = gx(1,:);
           y_vals = gy(:,1);
           this.means = zeros(flip(grd_sz));
           this.vars = zeros(flip(grd_sz));
           for i = 1:grd_sz(2)
              for j = 1:grd_sz(1)
                    pt = [x_vals(i), y_vals(j)];
                    this.means(i,j) = this.posteriorExpecteddB(pt);
                    this.vars(i,j) = this.calcVariance(pt);
              end
           end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotMeans2D
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot the expected value of the random variables represntning
        % channel gain dB at each grid location in the workspace.
        % Input:
        % this - reference to the PredictedChannel object
        function plotMeans2D(this)
           imagesc(this.cc.gx(1,:), this.cc.gy(:,1), this.means');
           set(gca, 'YDir', 'normal');
           c = colorbar;
           c.Label.String = 'Expected Channel Power (dBm)';
           xlabel('x (m)');
           ylabel('y (m)');
           title('Expected Channel Power');
        end
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % getEReqTXPowerW
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculates the expected required transmit power in Watts based on
        % certain quality of service and constellation parameters.
        % Input:
        % this - reference to the PredictedChannel object
        % qos - quality of service information, including BER, noise value,
        %       and spectral efficiency
        %
        % Output:
        % exp_req_tx_pwr - matrix of expected required transmit powers over
        %                   the workspace
        function exp_req_tx_pwr = getEReqTXPowerW(this, qos)
           [gx, ~] = this.cc.getMeshGrid();
           grid_sz = size(gx);
           exp_req_tx_pwr = zeros(flip(grid_sz));
           for i = 1:grid_sz(2)
              for j = 1:grid_sz(1)
                    expected_gamma = this.means(i,j);
                    variance = this.vars(i,j);
                    exp_req_tx_pwr(i,j) = qos.expReqTXPowerW(expected_gamma, variance);
              end
           end
 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotRequiredTXPower2D
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot the expected required transmit power, as calculated in
        % getEReqTXPowerW. Separating into another function to allow for
        % combining fields (e.g. at each location, take the min, max, or 
        % sum over multiple fields).
        %
        % Input:
        % this - reference to the PredictedChannel object
        % req_power - field of expected required power, as calculated in 
        %               getEReqTXPowerW 
        function plotRequiredTXPower2D(this, req_power)       
           imagesc(this.cc.gx(1,:), this.cc.gy(:,1), req_power');
           set(gca, 'YDir', 'normal');
           c = colorbar;
           c.Label.String = 'Expected Required TX Power, W';
           xlabel('x (m)');
           ylabel('y (m)');
           title('Expected Required TX Power');
        end
         
         
    end
    
    methods (Access = private)
        
        function mean = calcMean(this, K, pos)
            mean = this.estimatePL(pos) + K*this.obsDiff;
        end
        
        function [k, phi] = K(this, pos)
            phi = this.varPosWithObs(pos);
            k = phi'*this.UInv;
        end
        
        function var = calcVariance(this, pos)
            [k, phi] = this.K(pos);
            var = this.alpha() + this.rho()  - k*phi;
        end
        
        function cov = calcObsCov(this)
            obs_count = length(this.obsVals);
            
            diffs = zeros(obs_count);
            for i = 1:obs_count
               row = sqrt(sum((this.obsPos - this.obsPos(i,:)).^2,2));
               diffs(i,:) = row;
            end
            cov = this.cp.sigmaSH^2*exp(-diffs/this.beta());
        end
        
        function cov = varPosWithObs(this, pos)
            diffs = sqrt(sum((this.obsPos - pos).^2,2));
            cov = this.cp.sigmaSH^2*exp(-diffs/this.beta());
        end
        
        function mean = posteriorExpecteddB(this, pos)
            phi = this.varPosWithObs(pos);
            K = phi'*this.UInv;
            mean = this.calcMean(K, pos);
        end
        
        function estimatePLParams(this)
            obs_distances = sqrt(sum((this.obsPos - (this.cp.qBase)).^2, 2));
            est_PL_par = polyfit2D(obs_distances, this.obsVals);
            this.kPl_est = est_PL_par(1);
            this.nPl_est = est_PL_par(2);
            fprintf('Estimated PL params\n kPl: %.2f    nPl: %.2f', this.kPl_est, this.nPl_est);
        end
        
        function estimateSHMPParams(this)
            % weight_filter is used to improve the performance of the fit of the
            % correlation function of shadowing, because we have less measurements for
            % larger distances
            global weight_filter
            weight_filter = 1;
            [M,N] = size(this.cc.gx);
            per = this.obsCount/(M*N);
            [this.alpha_est, this.beta_est, this.rho_est] = SHMP_par_estimation(...
                this.obsPos(:,1), this.obsPos(:,2), this.obsDiff, 1, per, 40, 5);
            fprintf('Estimated SH and MP params with %.2f %% samples\n', per*100);
            fprintf('alpha: %.2f\n', this.alpha_est);
            fprintf('beta: %.2f\n', this.beta_est);
            fprintf('rho: %.2f\n', this.rho_est);
        end
        
        function b = beta(this)
           if this.useEstimatedParams
               b = this.beta_est;
           else
               b = this.cp.decorrSH;
           end
        end
        
        function a = alpha(this)
            if this.useEstimatedParams
               a = this.alpha_est;
           else
               a = this.cp.sigmaSH^2;
           end
        end
        
        function r = rho(this)
            if this.useEstimatedParams
               r = this.rho_est;
           else
               r = this.cp.sigmaMP^2;
           end
        end
    end
end

