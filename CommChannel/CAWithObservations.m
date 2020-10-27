%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predicted Channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a predicted model of the channel given observations

classdef CAWithObservations < handle
    
    properties (Access = public)
        cp;
        cc;
        
        obsCount;
        obsVals;
        obsGridPoss;% in gridCoordinates
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
        function this = CAWithObservations(channel_params, comm_channel, observation_pos,...
                                        channel_samples)
            this.cp = channel_params;
            this.cc = comm_channel;
            
            this.obsCount = length(channel_samples);
            this.obsGridPoss = observation_pos;
            this.obsVals = channel_samples;
            if this.useEstimatedParams
                this.estimatePLParams();
            end
            this.obsPL = this.estimatePL(this.obsGridPoss);
            this.obsDiff = this.obsVals - this.obsPL;
            if this.useEstimatedParams
                this.estimateSHMPParams();
            end
            this.obsCov = this.calcObsCov();
            this.UInv = (this.obsCov + this.rho()*eye(this.obsCount))^-1;
            
            this.setStats();
        end
        
        function b_res = scaledBeta(this)
           b_res = this.beta()*this.cc.res; 
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
        
        function gamma_PL_dBm_at_point = estimatePL(this, pt)
            if this.useEstimatedParams
                raw_point = this.cc.resPoint2Raw(pt);
                gamma_PL_dBm_at_point = this.kPl_est - 10*this.nPl_est*...
                            log10(sqrt(sum((this.cp.qBase - raw_point).^2, 2)));
            else
                gamma_PL_dBm_at_point = this.cc.getGammaPLdBAtPoint(pt);
            end
           
        end
        
        function mean = getMeanAtGridPoint(this, pt)
           mean = this.means(pt(1)+1, pt(2) + 1); 
        end
        
        function p_conn = posteriorPConn(this, pos, gamma_th)
            mean = this.means(pos(1) + 1, pos(2) + 1);
            variance = this.vars(pos(1) + 1, pos(2) + 1);
            
            p_conn = normcdf(gamma_th, mean, sqrt(variance), 'upper');
        end
        
        function plotPosteriors2D(this, gamma_th)
           res = this.cc.res;
           region = this.cc.region()*res;
           x_counts = region(1) - region(2) + 1;
           y_counts = region(3) - region(4) + 1;
           post_map = zeros([x_counts, y_counts]);
           for i = 1:x_counts
              for j = 1:y_counts
                    x_grid = i-1;
                    y_grid = j-1;
                    p_conn = this.posteriorPConn([x_grid, y_grid], gamma_th);
                    post_map(j,i) = p_conn;
              end
           end
         
           imagesc(this.cc.gx(1,:), this.cc.gy(:,1), post_map(:,:,3));
           set(gca, 'YDir', 'normal');
           c = colorbar;
           c.Label.String = 'Probability of Connectivity';
           xlabel('x (m)');
           ylabel('y (m)');
        end
        
        function conn_field = getConnectionField(this, p_th, gamma_th)
           res = this.cc.res;
           region = this.cc.region()*res;
           x_counts = region(1) - region(2) + 1;
           y_counts = region(3) - region(4) + 1;
           conn_field = zeros([x_counts, y_counts]);
           for i = 1:x_counts
              for j = 1:y_counts
                    x_grid = i-1;
                    y_grid = j-1;
                    conn = this.posteriorPConn([x_grid, y_grid], gamma_th) >= p_th;
                    conn_field(j,i) = conn;
              end
           end
        end
        
        function plotConnected2D(this, p_th, gamma_th)
           conn_field = this.getConnectionField(p_th);
           title = strcat('Connected Regions for \Gamma_{th} = ', sprintf('%d dBm', gamma_th),...
                           ', p_{conn} >= ', sprintf('%g',p_th));
           this.cc.plotField(conn_field,title);
        end
         
        function setStats(this)
           res = this.cc.res;
           region = this.cc.region()*res;
           x_counts = region(1) - region(2) + 1;
           y_counts = region(3) - region(4) + 1;
           this.means = zeros([x_counts, y_counts]);
           this.vars = zeros([x_counts, y_counts]);
           for i = 1:x_counts
              for j = 1:y_counts
                    x_grid = i-1;
                    y_grid = j-1;
                    this.means(i,j) = this.posteriorExpecteddB([x_grid, y_grid]);
                    this.vars(i,j) = this.calcVariance([x_grid, y_grid]);
              end
           end
        end
        
        function plotMeans2D(this)
           imagesc(this.cc.gx(1,:), this.cc.gy(:,1), this.means');
           set(gca, 'YDir', 'normal');
           c = colorbar;
           c.Label.String = 'Expected Channel Power (dBm)';
           xlabel('x (m)');
           ylabel('y (m)');
           title('Expected Channel Power');
        end
         
        function exp_req_tx_pwr = getEReqTXPowerW(this, qos)
           res = this.cc.res;
           region = this.cc.region()*res;
           x_counts = region(1) - region(2) + 1;
           y_counts = region(3) - region(4) + 1;
           exp_req_tx_pwr = zeros([x_counts, y_counts]);
           for i = 1:x_counts
              for j = 1:y_counts
                    expected_gamma = this.means(i,j);
                    variance = this.vars(i,j);
                    exp_req_tx_pwr(i,j) = qos.expReqTXPowerW(expected_gamma, variance);
              end
           end
 
        end
        
        function plotRequiredTXPower2D(this, req_power)
                    
           imagesc(this.cc.gx(1,:), this.cc.gy(:,1), req_power);
           set(gca, 'YDir', 'normal');
           c = colorbar;
           c.Label.String = 'Expected Required TX Power, mW';
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
               for j = i+1:obs_count
                   diff = norm(this.obsGridPoss(i,:) - this.obsGridPoss(j,:)); 
                   diffs(i,j) = diff;
                   diffs(j,i) = diff;
               end
            end
            cov = this.cp.sigmaSH^2*exp(-diffs/this.scaledBeta());
        end
        
        function cov = varPosWithObs(this, pos)
            obs_count = length(this.obsVals);
            diffs = zeros([obs_count,1]);
            for i = 1:obs_count
               diffs(i) = norm(this.obsGridPoss(i,:) - pos); 
            end
            
            cov = this.cp.sigmaSH^2*exp(-diffs/this.scaledBeta());
        end
        
        function mean = posteriorExpecteddB(this, pos)
            phi = this.varPosWithObs(pos);
            K = phi'*this.UInv;
            mean = this.calcMean(K, pos);
        end
        
        function estimatePLParams(this)
            res = this.cc.res;
            obs_distances = sqrt(sum((this.obsGridPoss - (this.cp.qBase*res)).^2, 2))/res;
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
            res = this.cc.res;
            [M,N] = size(this.cc.gx);
            per = this.obsCount/(M*N);
            [this.alpha_est, this.beta_est, this.rho_est] = SHMP_par_estimation(...
                this.obsGridPoss(:,1), this.obsGridPoss(:,2), this.obsDiff, 1/res, per, 40, 5);
            fprintf('Estimated SH and MP params with %.2f %% samples\n', per*100);
            fprintf('alpha: %.2f\n', this.alpha_est);
            fprintf('beta: %.2f\n', this.beta_est);
            fprintf('rho: %.2f\n', this.rho_est);
        end
        
            
    end
end

