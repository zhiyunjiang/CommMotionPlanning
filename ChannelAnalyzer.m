classdef ChannelAnalyzer < handle
    
    properties
        commChannel;
        channelParams;
        gammaTH;
        noMP;
        ricDist;
        scaledDecorrSH;
        
        %from arjun_python/fpd_with_mp.py
        g;
        u;
        delU;
        gammaSHVec;
    end
    
    methods
        function this = ChannelAnalyzer(comm_channel, gamma_TH, no_mp)
           this.commChannel = comm_channel;
           this.channelParams = comm_channel.channelParams;
           this.gammaTH = gamma_TH;
           this.scaledDecorrSH = this.channelParams.decorrSH*this.commChannel.res;
           
           if nargin == 2
               no_mp = 0;
           end
           this.noMP = no_mp;
           K_ric = this.channelParams.kRic;
           this.ricDist = makedist('Rician', 's', sqrt(K_ric/(K_ric + 1)), 'sigma',1/sqrt(2*(1+K_ric)));
           
           %Values that will needed to calculate probabilities later
           this.g = 2^16;
           this.u = 8*this.channelParams.sigmaSH;
           this.delU = 2*this.u/this.g;
           vec = 0:this.u-1;
           this.gammaSHVec = -this.u + this.delU*vec;
        end
        
        function prior = NoConnectionPrior(this, point, epsilon)
            
            if nargin == 2
               epsilon = 0; 
            end
            
           gamma_gap = this.gammaTH - this.commChannel.getGammaPLdBAtPoint(point) - epsilon;
           if this.noMP
               %then we're just dealing with path loss and shadowing. Check
               %probability that pathloss + shadowing < threshold ->
               %gamma_SH ~ N(0, sigma_SH)
               prior = normcdf(gamma_gap,0, this.channelParams.sigmaSH );
           else
               %now looking at the probabiltiy that gamma_SH + gamma_MP <
               %gamma_TH - gamma_PL. Will compute with fast fourier
               %transform, similar to Arjun's paper
               
           end
        end
        
        function no_conn_to_here = NoConnectionOnPath(this, path, method)
            if nargin == 2
                method = 1;
            end
            
            if this.noMP
                no_conn_to_here = NoConnectionOnPathNoMP(this, path, method);
            else
               %TODO - what should this look like for paths of arbitrary
               %shape?
               no_conn_to_here = 1;
            end
        end
        
        function exp_dist = ExpectedFPD(this, path, method)
           %can be implemented more efficiently, but will leave as this for now
           %expectation is conditional expectation, conditioned on no
           %connection at initial location (up-crossing first passage distance)
           if method <= 2
               %Multivariate-Gaussian base approach
               exp_dist = 0;

               path_dim = length(path);

               %handle the case of a path consisting of single point
               if min(size(path)) == 1
                    path_dim = 1;
               end

               no_conn_at_root = this.NoConnectionPrior(path(1,:));

               for i = 1:path_dim-1
                   current_dist = norm(path(i,:) - path(i+1,:));

                   p_no_conn_cond = this.NoConnectionOnPath(path(1:i,:), method)/no_conn_at_root;

                   exp_dist = p_no_conn_cond*current_dist + exp_dist;
               end
           elseif method == 3
               [g_pdf, d] = this.FPDPDFStraightlineNoMP(path);
               exp_dist = (g_pdf.*d)'*cumsum(d);
           end
        end
        
        function [approx_PMF, distances] = ApproxFPDPMF(this, path, method)
            path_dim = length(path);
            
           %handle the case of a path consisting of single point
           if min(size(path)) == 1
                path_dim = 1;
           end
           if method <= 2
           %methods based on joint pdf, use mvncdf of various
           %implementations
               no_conn_at_root = this.NoConnectionPrior(path(1,:));
               approx_PMF = zeros([path_dim, 1]);
               distances = zeros([path_dim, 1]);
               p_no_conn_prev = 1;

               for i = 2:path_dim
                   distances(i) = norm(path(i-1,:) - path(i,:));
                    p_no_conn_cond = this.NoConnectionOnPath(path(1:i,:), method)/no_conn_at_root;
                   approx_PMF(i) = p_no_conn_prev - p_no_conn_cond;
                   p_no_conn_prev = p_no_conn_cond;
               end
           elseif method == 3
               [g_pdf, distances] = this.FPDPDFStraightlineNoMP(path);
               approx_PMF = g_pdf.*distances;
           end
        end
        
        function [approx_PMF, distances] = ApproxFPDPMF2(this, step, num_steps, dsrc, thetasrc, thetaos, root)
            [g_pdf, distances] = this.FPDPDFStraightlineNoMP2(step, num_steps, dsrc, thetasrc, thetaos, root);
             approx_PMF = g_pdf.*distances;
        end
        
        function [g_aug, dists_aug] =  IterativeFPDPDFSLNoMP2(this, d, step, dsrc, thetasrc, thetaos, root, g, dists, epsilon)
            n = length(dists);
            
            %the root of the path
            if n == 0
                dists_aug = 0;
                g_new = 0;
            else
                d_new = step;
                dists_aug = [dists; d_new];
                
                d_l = d - 0.001*(step);
                w_simpson = SimpsonWeights(n, dists_aug(2:end));
                w_gnew = w_simpson(end);
                
                denom = (1 - 2*w_gnew*this.psiVolterra2(d, d_l, this.gammaTH, thetasrc, thetaos, dsrc, root));

                if n == 1 %first step along the path  no enough points to do numerical approximation of integral
                    g_new = (-2*this.psiVolterraUp2(d, thetasrc, thetaos, dsrc, root, epsilon))/denom;
                else
                    psi_current = zeros([n, 1]);
                    for i = 1:n
                        psi_current(i) = this.psiVolterra2(d, step*(i-1), this.gammaTH, thetasrc, thetaos, dsrc, root);
                    end
                    w_integral = w_simpson(1:end-1);
                    integrand_vals = psi_current.*g;
                    
                    g_new = (-2*this.psiVolterraUp2(d, thetasrc, thetaos, dsrc, root, epsilon) + 2*w_integral*integrand_vals)...
                        /denom;
                end
            end
            
            if g_new < 0
                g_new = 0;
            end
            g_aug = [g; g_new];
        end
        
        function plot_priors(this)
           cc = this.commChannel;
           gamma_gap = this.gammaTH - cc.getGammaPLdB(); 
           p_conn = 1 - normcdf(gamma_gap,0, this.channelParams.sigmaSH);
           [gx, gy] = cc.getMeshGrid();
           surf(gx, gy, p_conn, 'EdgeColor','none');
           
            % light
            % shading interp
            fnt_siz = 16;
            xlabel('x (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
            ylabel('y (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
            zlabel('P(connection)','FontSize', fnt_siz ,  'FontWeight','bold');
            axis tight
            grid on
            set(gca, 'FontSize', fnt_siz, 'FontWeight', 'bold');
        end
        
        
         function [g_aug, d_aug] = IterativeFPDPDFSLNoMP(this, g, d, root, prevs, current, epsilon)
            n = length(d);
            
            if n == 0
                %All probabilites conditioned on no connection at first
                %point
                d_aug = 0;
                g_new = 0;
            else
            
                d_new = norm(prevs(end,:) - current);
                d_aug = [d; d_new];
                
                d_prev = current - 0.01*(prevs(end,:) - current);
                w_simpson = SimpsonWeights(n, d_aug(2:end));
                w_gnew = w_simpson(end);
                
                denom = (1 - 2*w_gnew*this.psiVolterra(current, d_prev, this.gammaTH) );

                if n == 1 
                    g_new = (-2*this.psiVolterraUp(root, current, prevs(end,:), epsilon))/denom;
                else
                    psi_current = zeros([length(prevs), 1]);
                    for i = 1:length(prevs)
                        psi_current(i) = this.psiVolterra(current, prevs(i,:), this.gammaTH);
                    end
                    w_integral = w_simpson(1:end-1);
                    integrand_vals = flip(psi_current).*g;
                    
                    g_new = (-2*this.psiVolterraUp(root, current, prevs(end,:), epsilon) + 2*w_integral*integrand_vals)...
                        /denom;
                end
            end
            
            g_aug = [g; g_new];
            
        end
    end
    
    methods (Access = private)
        function no_conn_to_here = NoConnectionOnPathNoMP(this, path, method)
            if nargin == 2
                method = 1;
            end
            
            %handle the case of a path consisting of single point
            if min(size(path)) == 1
                gamma_gap = this.gammaTH - this.commChannel.getGammaPLdBAtPoint(path(1,:)); 
                no_conn_to_here = normcdf(gamma_gap, 0, this.channelParams.sigmaSH );
            elseif method <= 2
                %methods that try to calculate using the joint pdf of the
                %jointly gaussian RV's
                
                %matlab's multivariate normal CDF (mvncdf) can only handle 25
                %dimensions to only take the 25 most recent path spots. Eventually,
                %I'll want a better way to handle this
                %Don't need to do this with botev!
                
                path_dim = length(path);
                if path_dim > 25 && method == 2
                   path = path(path_dim - 24:path_dim,:);
                   path_dim = 25; 
                end

                gamma_gap = zeros([ path_dim, 1]);

                %build the covariance matrix
                pwr_sh = this.channelParams.sigmaSH^2;
                path_cov = diag(pwr_sh*ones([path_dim, 1]));

                for i = 1:path_dim
                   gamma_gap(i) = this.gammaTH - this.commChannel.getGammaPLdBAtPoint(path(i,:)); 
                   for j = i+1:path_dim
                      cov_entry =  pwr_sh * exp(-norm(path(i,:) - path(j,:)) /...
                                            (this.scaledDecorrSH));
                      path_cov(i,j) = cov_entry;
                      path_cov(j,i) = cov_entry;
                   end
                end
                if method == 1
                    %using Botev's improved mvncf
                    lower = -Inf*ones([path_dim,1]);
                    n_samples = path_dim*(5*10^4);
                    est = mvncdfBotev(lower, gamma_gap, path_cov, n_samples);
                    no_conn_to_here = est.prob;
                elseif method == 2
                    %use matlab's mvncdf, truncate path to 25 most recent
                    no_conn_to_here = mvncdf(gamma_gap, [], path_cov);
                end
            elseif method ==3
                [g_pdf, d] = this.FPDPDFStraightlineNoMP(path);
                g_cdf = cumsum(g_pdf.*d);
                no_conn_to_here = 1 - g_cdf(end);
            end
        end
        
        %FPDPDFStraightlineNoMP - calculate the PMF for epsilon-upcrossing
        %FPD. Based on Di Nardo.
        %%%%%%%%%%%%%%%%%%%%%%%%
        %Inputs
        %%%%%%%%%%%%%%%%%%%%%%%%
        % this - current ChannelAnalyzer object
        % path - Path object
        %%%%%%%%%%%%%%%%%%%%%%%%
        %Outputs
        %%%%%%%%%%%%%%%%%%%%%%%%
        % fpd_pmf - n-1 dimensional vector, with fpd_pmf(i) being the
        % probability of fpd ocurring at the path's i+1th point
        % distances - n-1 dimensional vector, with distances(i) giving the
        % distance traveled upon reaching the path's i+1th point
        function [g_pdf, d] = FPDPDFStraightlineNoMP(this, path, epsilon)
            if nargin == 2
               epsilon = 0.0001; 
            end
            
            path_dim = length(path);
            if min(size(path)) == 1
                path_dim = 1;
            end
            
            root = path(1,:);
            [g_pdf,d] = this.IterativeFPDPDFSLNoMP([], [], root,...
                    [], [], epsilon);
            for i = 2:path_dim
                [g_pdf, d] = this.IterativeFPDPDFSLNoMP(g_pdf, d, root,...
                    path(1:i-1,:), path(i,:), epsilon);
            end
        end
        
        function [gs, dists]  = FPDPDFStraightlineNoMP2(this, step, num_steps, dsrc, thetasrc, thetaos, root, epsilon)
            if nargin == 7
               epsilon = 0.0001; 
            end
            
            
            gs = []; dists = [];
            for i = 1:num_steps + 1
                d = step*(i-1);
                [gs, dists] =  this.IterativeFPDPDFSLNoMP2(d, step, dsrc, thetasrc, thetaos, root, gs, dists, epsilon);
            end
        end
        
        function p_v = psiVolterraUp(this, root, current_point, prev_point, epsilon)
            cc = this.commChannel; cp = this.channelParams;
            
            m_r = cc.getGammaPLdBAtPoint(root);
            m_d = cc.getGammaPLdBAtPoint(current_point);
            d2r = norm(root - current_point);
            
            g_th = this.gammaTH;
            m_prime = this.plPrime( current_point, prev_point);
            p_no_conn_at_root = this.NoConnectionPrior(root, epsilon);
            sigma_sh = cp.sigmaSH;
            pwr_sh = sigma_sh^2;
            
            
            f1 = normpdf(g_th - m_r - epsilon, 0, sigma_sh);
            m_cond = m_d + exp(-d2r/this.scaledDecorrSH)*(g_th - epsilon - m_r);
            var_cond = pwr_sh*(1 - exp(-2*d2r/this.scaledDecorrSH));
            if var_cond == 0
                f2 = 0;
            else
                f2 = normpdf(g_th, m_cond, sqrt(var_cond));
            end
            f3 = normpdf(g_th - m_d, 0, cp.sigmaSH);
            
            ups = (g_th - epsilon - m_r - exp(-d2r/this.scaledDecorrSH)*(g_th - m_d))/...
                sqrt(2*pwr_sh*(1 - exp(-2*d2r/this.scaledDecorrSH)));
            
            p_v = ( 1/(2*p_no_conn_at_root) )*( (-2*pwr_sh/this.scaledDecorrSH)*exp(-d2r/this.scaledDecorrSH)*f1*f2...
            + 0.5*f3*(1+erf(ups))* (-m_prime - (g_th - m_d)/this.scaledDecorrSH));
            
            %dinardo = this.psiVolterraUpdiNardo(root, current_point, prev_point, epsilon);
            %p_v = dinardo;
        end
        
        function pvu =  psiVolterraUpdiNardo(this, root, current_point, prev_point, epsilon)
           cc = this.commChannel;
           d = norm(root - current_point);
           h1_r = this.h1Cov(0);
           h2_r = this.h2Cov(0);
           h1_d = this.h1Cov(d);
           h2_d = this.h2Cov(d);
           h1p_d = this.h1CovPrime(d);
           h2p_d = this.h2CovPrime(d);
           S = this.gammaTH;
           
           m_r = cc.getGammaPLdBAtPoint(root);
           m_d = cc.getGammaPLdBAtPoint(current_point);
           
           sigma_sh = this.channelParams.sigmaSH;
           pwr_sh = sigma_sh^2;
           m_prime = this.plPrime( current_point, prev_point);
           m_cond = m_d + exp(-d/this.scaledDecorrSH)*(S - epsilon - m_r);
           var_cond = pwr_sh*(1 - exp(-2*d/this.scaledDecorrSH));
           
           
           erf_arg = (S-epsilon - m_r)/sqrt(2*h1_r*h2_r);
           erf_arg2 = sqrt(h1_d/(2*h1_r*(h1_d*h2_r - h1_r*h2_d)))...
               *(S - epsilon - m_r - (S-m_d)*(h1_r/h1_d));
           
           pvu = (1/(1+erf(erf_arg))) * ( (h1_r/h1_d)*(h1_d*h2p_d - h1p_d*h2_d)...
               * normpdf(S-epsilon, m_r, sigma_sh) * normpdf(S, m_cond, sqrt(var_cond)) ...
               + 0.5*normpdf(S, m_d, sigma_sh)*(1+erf(erf_arg2))*(-m_prime - (h1p_d/h1_d)*(S - m_d)));
        end
        
        function p_v = psiVolterra(this, current_point, prev_point, y)
            cc = this.commChannel;
            l = 0;%preivous points "distance"
            %current "distance" is just distance from previous to current
            d = norm(prev_point - current_point);

            m_d = cc.getGammaPLdBAtPoint(current_point);
            m_l = cc.getGammaPLdBAtPoint(prev_point);
            
            m_prime = this.plPrime( current_point, prev_point);
            S = this.gammaTH; S_prime = 0;
            
            h1_d = this.h1Cov(d);
            h1p_d = this.h1CovPrime(d);
            h2_d = this.h2Cov(d);
            h2p_d = this.h2CovPrime(d);
            h1_l = this.h1Cov(l);
            h2_l = this.h2Cov(l);
            m_cond = m_d + (h2_d/h2_l)*(y - m_l);
            var_cond = h2_d*(h1_d - (h2_d/h2_l))*h1_l;
            
            f = @(s) normpdf(s, m_cond, sqrt(var_cond));
            
            p_v = ((S_prime - m_prime)/2 - ((S-m_d)/2)*((h1p_d*h2_l-h2p_d*h1_l)/...
                (h1_d*h2_l-h2_d*h1_l)) - ((y - m_l)/2) * ((h2p_d*h1_d-h2_d*h1p_d) /...
                (h1_d*h2_l-h2_d*h1_l))) * f(S);
        end
        
        function p_v = psiVolterraUp2(this, d, thetasrc, thetaos, dsrc, root, epsilon)
           cc = this.commChannel;
          
           h1_r = this.h1Cov(0);
           h2_r = this.h2Cov(0);
           h1_d = this.h1Cov(d);
           h2_d = this.h2Cov(d);
           h1p_d = this.h1CovPrime(d);
           h2p_d = this.h2CovPrime(d);
           S = this.gammaTH;
           
           m_r = cc.getGammaPLdBAtPoint(root);
           current_point = root - d*[cos(thetasrc + thetaos), sin(thetasrc + thetaos)];
           m_d = cc.getGammaPLdBAtPoint(current_point);
            
           sigma_sh = this.channelParams.sigmaSH;
           pwr_sh = sigma_sh^2;
           m_prime = this.plPrime2( d, thetasrc, dsrc);
           m_cond = m_d + exp(-d/this.scaledDecorrSH)*(S - epsilon - m_r);
           var_cond = pwr_sh*(1 - exp(-2*d/this.scaledDecorrSH));
           
           
           erf_arg = (S-epsilon - m_r)/sqrt(2*h1_r*h2_r);
           erf_arg2 = sqrt(h1_d/(2*h1_r*(h1_d*h2_r - h1_r*h2_d)))...
               *(S - epsilon - m_r - (S-m_d)*(h1_r/h1_d));
           
           p_v = (1/(1+erf(erf_arg))) * ( (h1_r/h1_d)*(h1_d*h2p_d - h1p_d*h2_d)...
               * normpdf(S-epsilon, m_r, sigma_sh) * normpdf(S, m_cond, sqrt(var_cond)) ...
               + 0.5*normpdf(S, m_d, sigma_sh)*(1+erf(erf_arg2))*(-m_prime - (h1p_d/h1_d)*(S - m_d)));
        end
        
        function p_v =  psiVolterra2(this, d, l, y, thetasrc, thetaos, dsrc, root)
            cc = this.commChannel;
            
            current_point = root + d*[-cos(thetasrc - thetaos), sin(thetasrc - thetaos)];
            m_d = cc.getGammaPLdBAtPoint(current_point);
            prev_point = root + l*[-cos(thetasrc - thetaos), sin(thetasrc - thetaos)];
            m_l = cc.getGammaPLdBAtPoint(prev_point);
            
            m_prime = this.plPrime2( d, thetasrc, dsrc);
            S = this.gammaTH; S_prime = 0;
            
            h1_d = this.h1Cov(d);
            h1p_d = this.h1CovPrime(d);
            h2_d = this.h2Cov(d);
            h2p_d = this.h2CovPrime(d);
            h1_l = this.h1Cov(l);
            h2_l = this.h2Cov(l);
            m_cond = m_d + (h2_d/h2_l)*(y - m_l);
            var_cond = h2_d*(h1_d - (h2_d/h2_l))*h1_l;
            
            f = @(s) normpdf(s, m_cond, sqrt(var_cond));
            
            p_v = ((S_prime - m_prime)/2 - ((S-m_d)/2)*((h1p_d*h2_l-h2p_d*h1_l)/...
                (h1_d*h2_l-h2_d*h1_l)) - ((y - m_l)/2) * ((h2p_d*h1_d-h2_d*h1p_d) /...
                (h1_d*h2_l-h2_d*h1_l))) * f(S);
        end
        
        function m_prime = plPrime(this, current_point, prev_point)
            cp = this.channelParams;
            qb = (cp.qBase - ...
                [this.commChannel.region(2), this.commChannel.region(4)])*this.commChannel.res;
            d2b = norm(current_point - qb); 
            
            v1 = qb - prev_point; v2 = qb - current_point - prev_point;
            cos_theta = -max(min(dot(v1,v2)/(norm(v1)*norm(v2)),1),-1);
            m_prime = -10*cp.nPL/(d2b*log(10))*cos_theta;
        end
        
        function m_prime = plPrime2(this, d, thetasrc, dsrc)
            cp = this.channelParams;
           m_prime = -10*cp.nPL*log10(exp(1))*(d - dsrc*cos(thetasrc))...
               /(dsrc^2 + d^2 - 2*dsrc*d*cos(thetasrc));
        end
        
        function h1_cov = h1Cov(this, d)
            cp = this.channelParams;
            h1_cov = cp.sigmaSH*exp(d/this.scaledDecorrSH);
        end
        
        function h1_cov_prime = h1CovPrime(this, d)
            h1_cov_prime = this.h1Cov(d)/this.scaledDecorrSH;
        end
        
        function h2_cov = h2Cov(this, d)
            cp = this.channelParams;
            h2_cov = cp.sigmaSH*exp(-d/this.scaledDecorrSH);
        end
        
        function h2_cov_prime = h2CovPrime(this, d)
            h2_cov_prime = -this.h2Cov(d)/this.scaledDecorrSH;
        end
        
        %%%%%%%
        %To Troubleshoot, implementFPDPDF for a path parameterized as in
        %oroginal paper - starting point x0 and angle theta w.r.t to source
        
        
        
    end
    
end

