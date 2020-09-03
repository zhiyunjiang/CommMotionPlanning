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
         
           imagesc(0,0, post_map(:,:,3));
           c = colorbar;
           c.Label.String = 'Probability of Connectivity';
           xlabel(sprintf('x (%g m)', 1/this.cc.res));
           ylabel(sprintf('y (%g m)', 1/this.cc.res));
        end
        
         function plotConnected2D(this, p_th)
           res = this.cc.res;
           region = this.cc.region()*res;
           x_counts = region(1) - region(2) + 1;
           y_counts = region(3) - region(4) + 1;
           post_map = zeros([x_counts, y_counts,3]);
           for i = 1:x_counts
              for j = 1:y_counts
                    x_grid = i-1;
                    y_grid = j-1;
                    conn = this.posteriorPConn([x_grid, y_grid]) >= p_th;
                    post_map(j,i,:) = [x_grid, y_grid, conn];
              end
           end
         
           imagesc(0,0, post_map(:,:,3));
           colormap([0.5 0.5 0.5; 0 0 0.5]);
           colorbar;
           colorbar('YTick',[0.25 0.75],'YTicklabel',{'Disconnected', 'Connected'},...
                    'FontSize', 7, 'FontName', 'Calibri')
           title(strcat('Connected Regions for \Gamma_{th} = ', sprintf('%d dBm', this.gammaTH),...
                           ', p_{conn} >= ', sprintf('%g',p_th))) ;
           xlabel(sprintf('x (%g m)', 1/this.cc.res));
           ylabel(sprintf('y (%g m)', 1/this.cc.res));
         end
        
         function plotMeans2D(this)
           
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
                    post_map(j,i,:) = [x_grid, y_grid, p_conn];
              end
           end
         
           imagesc(0, 0, post_map(:,:,3));
           c = colorbar;
           c.Label.String = 'Expected Channel Power, dBm';
           xlabel(sprintf('x (%g m)', 1/this.cc.res));
           ylabel(sprintf('y (%g m)', 1/this.cc.res));
           title('Expected Channel Power');
         end
        
         function plotRequiredTXPower2D(this, receiver_noise, BER, R)
           K =  -1.5/log(5*BER);
           res = this.cc.res;
           region = this.cc.region()*res;
           x_counts = region(1) - region(2) + 1;
           y_counts = region(3) - region(4) + 1;
           post_map = zeros([x_counts, y_counts,3]);
           for i = 1:x_counts
              for j = 1:y_counts
                    x_grid = i-1;
                    y_grid = j-1;
                    expected_gamma = this.posteriorExpecteddB([x_grid, y_grid]);
                    CNR_lin = 10.^(expected_gamma/10) / receiver_noise;
                    req_power = ((2^R - 1)/K)*(1./CNR_lin);
                    post_map(j,i,:) = [x_grid, y_grid, req_power];
              end
           end
         
           imagesc(0, 0, post_map(:,:,3));
           c = colorbar;
           c.Label.String = 'Expected Required TX Power, mW';
           xlabel(sprintf('x (%g m)', 1/this.cc.res));
           ylabel(sprintf('y (%g m)', 1/this.cc.res));
           title(sprintf('Expected Required TX Power, BER = %d', BER))
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

