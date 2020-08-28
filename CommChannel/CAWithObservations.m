classdef CAWithObservations< handle
    %CAWithOBservations  - All pos corrdinates in grid, not raw coordinates
    
    properties
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
        
        gammaTH
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
        end
        
        function p_conn = posteriorPConn(this, pos)
            phi = this.varPosWithObs(pos);
            K = phi'*this.UInv;
            mean = this.calcMean(K, pos);
            variance = this.cp.sigmaSH^2 + this.cp.sigmaMP^2  - K*phi;
            
            p_conn = normcdf(this.gammaTH, mean, sqrt(variance), 'upper');
        end
        
        function mean = posteriorExpecteddB(this, pos)
            phi = this.varPosWithObs(pos);
            K = phi'*this.UInv;
            mean = this.calcMean(K, pos);
        end
        
        function plotPosteriors(this)
           
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
                    post_map(i,j,:) = [x_grid, y_grid, p_conn];
              end
           end
         
           h = imagesc(post_map(:,:,3));
           set(h, 'XData', [0, x_counts-1]);
           set(h, 'YData', [0, y_counts-1]);
           colorbar;
           zlabel('Posterior Probability of Connectivity')
           
        end
        
         function plotMeans(this)
           
           res = this.cc.res;
           region = this.cc.region()*res;
           x_counts = region(1) - region(2) + 1;
           y_counts = region(3) - region(4) + 1;
           post_map = zeros([x_counts, y_counts,3]);
           for i = 1:x_counts
              for j = 1:y_counts
                    x_grid = i-1;
                    y_grid = j-1;
                    p_conn = this.posteriorExpecteddB([x_grid, y_grid]);
                    post_map(i,j,:) = [x_grid, y_grid, p_conn];
              end
           end
         
           h = imagesc(post_map(:,:,3)');
           set(h, 'XData', [0, x_counts-1]);
           set(h, 'YData', [0, y_counts-1]);
           colorbar;
           zlabel('Posterior Expected PL dB')
           
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
    end
end

