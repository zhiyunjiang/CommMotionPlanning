classdef CAWithObservations< handle
    %CAWithOBservations  - All pos corrdinates in grid, not raw coordinates
    
    properties (Access = public)
        cp;
        cc;
        scaledDecorrSH;
        
        obsCount;
        obsVals;
        obsPoss;
        obsPL;
        obsDiff;
        obsCov;
        UInv;
        
        gammaTH;
        
         means;
    end
    
    methods (Access = public)
        function this = CAWithObservations(channel_params, comm_channel, observation_pos,...
                                        observation_vals, gamma_th)
            this.cp = channel_params;
            this.cc = comm_channel;
            this.scaledDecorrSH = this.cp.decorrSH*this.cc.res;
            
            this.obsCount = length(observation_vals);
            this.obsPoss = observation_pos;
            this.obsVals = observation_vals;
            this.obsPL = this.cc.getGammaPLdBAtPoint(this.obsPoss);
            this.obsDiff = this.obsVals - this.obsPL;
            this.obsCov = this.calcObsCov();
            this.UInv = (this.obsCov + this.cp.sigmaMP*eye(this.obsCount))^-1;
            
            this.gammaTH = gamma_th;
            
            this.setMeans();
        end
        
        function mean = getMeanAtGridPoint(this, pt)
           mean = this.means(pt(1)+1, pt(2) + 1); 
        end
        
        function p_conn = posteriorPConn(this, pos)
            phi = this.varPosWithObs(pos);
            K = phi'*this.UInv;
            mean = this.means(pos(1) + 1, pos(2) + 1);
            variance = this.cp.sigmaSH^2 + this.cp.sigmaMP^2  - K*phi;
            
            p_conn = normcdf(this.gammaTH, mean, sqrt(variance), 'upper');
        end
        
        function plotPosteriors2D(this)
           res = this.cc.res;
           region = this.cc.region()*res;
           x_counts = region(1) - region(2) + 1;
           y_counts = region(3) - region(4) + 1;
           post_map = zeros([x_counts, y_counts,3]);
           for i = 1:x_counts
              for j = 1:y_counts
                    x_grid = i-1;
                    y_grid = j-1;
                    p_conn = this.posteriorPConn([x_grid, y_grid]);
                    post_map(j,i,:) = [x_grid, y_grid, p_conn];
              end
           end
         
           imagesc(this.cc.gx(1,:), this.cc.gy(:,1), post_map(:,:,3));
           set(gca, 'YDir', 'normal');
           c = colorbar;
           c.Label.String = 'Probability of Connectivity';
           xlabel('x (m)');
           ylabel('y (m)');
        end
        
        function conn_field = getConnectionField(this, p_th)
           res = this.cc.res;
           region = this.cc.region()*res;
           x_counts = region(1) - region(2) + 1;
           y_counts = region(3) - region(4) + 1;
           conn_field = zeros([x_counts, y_counts]);
           for i = 1:x_counts
              for j = 1:y_counts
                    x_grid = i-1;
                    y_grid = j-1;
                    conn = this.posteriorPConn([x_grid, y_grid]) >= p_th;
                    conn_field(j,i) = conn;
              end
           end
        end
        
        function plotConnected2D(this, p_th)
           conn_field = this.getConnectionField(p_th);
           title = strcat('Connected Regions for \Gamma_{th} = ', sprintf('%d dBm', this.gammaTH),...
                           ', p_{conn} >= ', sprintf('%g',p_th));
           this.cc.plotField(conn_field,title);
        end
         
        function setMeans(this)
           res = this.cc.res;
           region = this.cc.region()*res;
           x_counts = region(1) - region(2) + 1;
           y_counts = region(3) - region(4) + 1;
           this.means = zeros([x_counts, y_counts]);
           for i = 1:x_counts
              for j = 1:y_counts
                    x_grid = i-1;
                    y_grid = j-1;
                    p_conn = this.posteriorExpecteddB([x_grid, y_grid]);
                    this.means(j,i) = p_conn;
              end
           end
        end
        
        function plotMeans2D(this)
           imagesc(this.cc.gx(1,:), this.cc.gy(:,1), this.means);
           set(gca, 'YDir', 'normal');
           c = colorbar;
           c.Label.String = 'Expected Channel Power (dBm)';
           xlabel('x (m)');
           ylabel('y (m)');
           title('Expected Channel Power');
        end
         
        function exp_req_tx_pwr = calcReqTXPwr(this, receiver_noise, BER, R)
           K =  -1.5/log(5*BER);
           res = this.cc.res;
           region = this.cc.region()*res;
           x_counts = region(1) - region(2) + 1;
           y_counts = region(3) - region(4) + 1;
           exp_req_tx_pwr = zeros([x_counts, y_counts]);
           for i = 1:x_counts
              for j = 1:y_counts
                    expected_gamma = this.means(i,j);
                    CNR_lin = 10.^(expected_gamma/10) / receiver_noise;
                    req_power = ((2^R - 1)/K)*(1./CNR_lin);
                    exp_req_tx_pwr(j,i) = req_power;
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
            mean = this.cc.getGammaPLdBAtPoint(pos) + K*this.obsDiff;
        end
        
        function cov = calcObsCov(this)
            obs_count = length(this.obsVals);
            diffs = zeros(obs_count);
            for i = 1:obs_count
               for j = i+1:obs_count
                   diff = norm(this.obsPoss(i,:) - this.obsPoss(j,:)); 
                   diffs(i,j) = diff;
                   diffs(j,i) = diff;
               end
            end
            cov = this.cp.sigmaSH^2*exp(-diffs/this.scaledDecorrSH);
        end
        
        function cov = varPosWithObs(this, pos)
            obs_count = length(this.obsVals);
            diffs = zeros([obs_count,1]);
            for i = 1:obs_count
               diffs(i) = norm(this.obsPoss(i,:) - pos); 
            end
            
            cov = this.cp.sigmaSH^2*exp(-diffs/this.scaledDecorrSH);
        end
        
        function mean = posteriorExpecteddB(this, pos)
            phi = this.varPosWithObs(pos);
            K = phi'*this.UInv;
            mean = this.calcMean(K, pos);
        end
    end
end

