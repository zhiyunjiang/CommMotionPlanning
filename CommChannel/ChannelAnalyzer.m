%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ChannelAnalyzer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates firts passage distance related values for path through a given
% channel
%
% Methods
% 1 - full implementation of Arjun's proposed iterative algorithm for
%       calculating FPD - HANDLES MULTIPATH
% 2 - Iterative calculation akin to Arjun's implementation, based on Di
%       Nardo's paper - does not account for multipath
classdef ChannelAnalyzer < handle
    
    properties
        commChannel;
        channelParams;
        gammaTH;
        noMP;
        ricDist;
        riceCDFLUT;%Rician dist look up table
        LUTRes = 1000; % the number of entries to build out in the LUT
        LUTMax;
        
        %from arjun_python/fpd_with_mp.py
        g;
        u;
        delU;
        gammaSHVec;
    end
    
    methods
        function this = ChannelAnalyzer(comm_channel, gamma_TH, no_mp)
           this.commChannel = comm_channel;
           this.channelParams = comm_channel.cp;
           this.gammaTH = gamma_TH;
           
           if nargin == 2
               no_mp = 0;
           end
           this.noMP = no_mp;
           K_ric = this.channelParams.kRic;
           this.ricDist = makedist('Rician', 's', sqrt(K_ric/(K_ric + 1)), 'sigma',1/sqrt(2*(1+K_ric)));
           this.computeRiceCDFLUT();
           %Values that will needed to perform numerical integration later
           %number of samples to use
           %this.g = 2^16;
           this.g = 2^14;
           %samples will run from 8 std dev below to 8 above mean (0)
           this.u = 8*this.channelParams.sigmaSH;
           %distance between each sample
           this.delU = 2*this.u/this.g;
           vec = 0:this.g-1;
           %sampling values of gammaSH
           this.gammaSHVec = -this.u + this.delU*vec;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % simulateFPD
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculates the first passage distance along a path for several
        % realizations of a comm channel
        % Input:
        % this - reference to the ChannelAnalyzer object
        % path - nX2 matrix with each row representing a waypoint along the
        %           path
        % n_sims - the number of times to simulate the channel
        % is_markov - if true, treat the channel along the path as
        %               Markovian. This allows for faster simulation, but 
        %               is a poorer model of reality
        %
        % Output:
        % fpd_dists - vector of n_sims distances, the first passage distance
        %               for each simulation
        function fpd_dists = simulateFPD(this, path, n_sims, is_markov)
            cc = this.commChannel;
            sims_run = 0;
            fpd_dists = zeros([n_sims, 1]);

            while sims_run < n_sims
               %generate a new realization of the path
               gamma_sim = cc.simulatePath(path, this.noMP, is_markov); 
               %check if root is not connected (event we're conditioned on)
               if gamma_sim(1) > this.gammaTH
                   %if not, we're not interested, simulate again.
                   continue;
               end

               current_dist = 0;
               for j=1:length(path)-1
                  if gamma_sim(j) < this.gammaTH
                    current_dist = current_dist + norm(path(j,:) - path(j+1,:));
                  else
                      break;
                  end
               end
               
               sims_run = sims_run + 1;
               if mod(sims_run, 20) == 0
                  fprintf('%d simulations complete\n', sims_run); 
               end
               fpd_dists(sims_run) = current_dist;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NoConnectionPrior
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find the marginal pronability of no connection at a point
        % Input:
        % this - reference to the ChannelAnalyzer object
        % pt - x,y point of interest
        % epsilon - arbitrarily small value used in theoretical derivation
        %
        % Output:
        % fpd_dists - vector of n_sims distances, the first passage distance
        %               for each simulation
        function marginal = NoConnectionMarginal(this, pt, epsilon)
            
            if nargin == 2
               epsilon = 0; 
            end
            
           gamma_gap = this.gammaTH - this.commChannel.getGammaPLdBAtPt(pt) - epsilon;
           if this.noMP
               %then we're just dealing with path loss and shadowing. Check
               %probability that pathloss + shadowing < threshold ->
               %gamma_SH ~ N(0, sigma_SH)
               marginal = normcdf(gamma_gap,0, this.channelParams.sigmaSH );
           else
               %now looking at the probabiltiy that gamma_SH + gamma_MP <
               %gamma_TH - gamma_PL. Will compute with fast fourier
               %transform, similar to Arjun's paper
               J0 = this.J0(pt);
               marginal = this.IntegrateJ(J0);
           end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % NoConnectionOnPath
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find the probability that connection never established while
        % traversing the path.
        % Input:
        % this - reference to the ChannelAnalyzer object
        % path - the path being analyzed
        % method - 1: Arjun's iterative approximation using FFT,
        %          2: Iterative Arjun's w/o multipath
        %
        % Output:
        % no_conn - probability of no connection to along the path
        function no_conn = NoConnectionOnPath(this, path, method)
            if nargin == 2
                method = 1;
            end

            [approx_PMF, d] = ApproxFPDPMF(this, path, method);
            no_conn = approx_PMF'*d;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ExpectedFPD
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculates the first passage distance
        % Input:
        % this - reference to the ChannelAnalyzer object
        % path - the path being analyzed
        % method - 1: Arjun's iterative approximation using FFT,
        %          2: Iterative Arjun's w/o multipath
        %
        % Output:
        % exp_dist - Expected first passage distance
        function exp_dist = ExpectedFPD(this, path, method)
            if nargin == 2
                method = 1;
            end
           [approx_PMF, distances] = ApproxFPDPMF(this, path, method);
           exp_dist = approx_PMF'*distances;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ApproxFPDPMF
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate the discretized PDF (i.e. PMF) of first pasage distance
        % Input:
        % this - reference to the ChannelAnalyzer object
        % path - the path being analyzed
        % method - 1: Arjun's iterative approximation using FFT,
        %          2: Iterative Arjun's w/o multipath
        %
        % Output:
        % approx_PMF - Probability mass at each distance
        % distances - distance associated with each probability mass in
        %               approx_PMF. First value is always 0.
        function [approx_PMF, distances] = ApproxFPDPMF(this, path, method)
            if nargin == 2
                method = 1;
            end
            
           if method == 1
               [approx_PMF, distances] = this.FPDPMFStraightline(path);
           else
               [approx_PMF, distances] = this.FPDPDFStraightlineNoMP(path);
           end
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
        
    end
    
    methods (Access = private)
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
                    integrand_vals = psi_current.*g;
                    
                    g_new = (-2*this.psiVolterraUp(root, current, prevs(end,:), epsilon) + 2*w_integral*integrand_vals)...
                        /denom;
                end
            end
            
            if g_new < 0
                g_new =0;
            end
            
            g_aug = [g; g_new];
            
        end
        
        function J_0 = J0(this, root)
            %Approximate integral from -Inf to Inf  of P(gamma_PL + gamma_SH + gamma_MP < gamma_TH)
            %using recursive J formulation as in Arjun's paper
            gamma_pl_root = this.commChannel.getGammaPLdBAtPt(root);
            J_0 = normpdf(this.gammaSHVec, 0, this.channelParams.sigmaSH)...
                .*this.mpCDF(this.gammaTH - gamma_pl_root - this.gammaSHVec);
        end
        
        function J_next = ItterativeJNextFromPoint(this, J_prev, point, step_size)
            gamma_pl = this.commChannel.getGammaPLdBAtPt(point);
           J_next = this.ItterativeJNext(J_prev, gamma_pl, step_size); 
        end
        
        function J_next = ItterativeJNext(this, J_prev, gamma_pl, step_size)
            cp = this.channelParams;
            rho = exp(-step_size/this.channelParams.decorrSH);
            sig = sqrt(cp.sigmaSH^2 *(1 - rho^2));
            phi = normpdf(this.gammaSHVec, 0, sig);
            G = this.g;
            
            J_prev_contracted = zeros([this.g, 1]);
            gfloor = floor(G/2);
            for i = 1:G
               J_prev_index = round((1 - rho^(-1) )*gfloor + (i/rho) ) + 1;
               if J_prev_index >= 1 && J_prev_index <= G 
                  J_prev_contracted(i) = J_prev(J_prev_index);
               end
            end
            
            J_next = (1/rho)*this.delU*conv(phi, J_prev_contracted, 'same')...
                .*this.mpCDF(this.gammaTH - gamma_pl - this.gammaSHVec);
        end
        
        function result = IntegrateJ(this, J)
            %estimate integral using Riemann sum
            result = sum(J)*this.delU;
        end
        
        function [p_no_conn_to_here, d] = NoConnAlongSL(this, path)
            path_dim = length(path);
            if min(size(path)) == 1
                path_dim = 1;
            end
            p_no_conn_to_here = zeros([path_dim, 1]);
            d = zeros([path_dim, 1]);
            J_0 = this.J0(path(1,:));

            p_no_conn_to_here(1) = this.IntegrateJ(J_0);
            d(1) = 0;
            
            J_prev = J_0;
            for i = 2:path_dim
                gamma_pl = this.commChannel.getGammaPLdBAtPt(path(i,:));
                d(i) = norm(path(i-1,:) - path(i,:));
                J_prev = this.ItterativeJNext(J_prev, gamma_pl, d(i));
                p_no_conn_to_here(i) = this.IntegrateJ(J_prev);
            end
            p_no_conn_to_here = p_no_conn_to_here/p_no_conn_to_here(1);
        end
        
        function [g, d] = FPDPMFStraightline(this, path)
            
           [p_no_conn_to_here, d] = this.NoConnAlongSL(path);
            
            g = [0; p_no_conn_to_here(1:end-1) - p_no_conn_to_here(2:end)];
        
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
        
        function p_v = psiVolterraUp(this, root, current_point, prev_point, epsilon)
            cc = this.commChannel; cp = this.channelParams;
            beta = this.channelParams.decorrSH;
            
            m_r = cc.getGammaPLdBAtPt(root);
            m_d = cc.getGammaPLdBAtPt(current_point);
            d2r = norm(root - current_point);
            
            g_th = this.gammaTH;
            m_prime = this.plPrime( current_point, prev_point);
            p_no_conn_at_root = this.NoConnectionMarginal(root, epsilon);
            sigma_sh = cp.sigmaSH;
            pwr_sh = sigma_sh^2;
            
            
            f1 = normpdf(g_th - m_r - epsilon, 0, sigma_sh);
            m_cond = m_d + exp(-d2r/beta)*(g_th - epsilon - m_r);
            var_cond = pwr_sh*(1 - exp(-2*d2r/beta));
            if var_cond == 0
                f2 = 0;
            else
                f2 = normpdf(g_th, m_cond, sqrt(var_cond));
            end
            f3 = normpdf(g_th - m_d, 0, cp.sigmaSH);
            
            ups = (g_th - epsilon - m_r - exp(-d2r/beta)*(g_th - m_d))/...
                sqrt(2*pwr_sh*(1 - exp(-2*d2r/beta)));
            
            p_v = ( 1/(2*p_no_conn_at_root) )*( (-2*pwr_sh/beta)*exp(-d2r/beta)*f1*f2...
            + 0.5*f3*(1+erf(ups))* (-m_prime - (g_th - m_d)/beta));
            
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
           beta = this.channelParams.decorrSH;
           
           m_r = cc.getGammaPLdBAtPt(root);
           m_d = cc.getGammaPLdBAtPt(current_point);
           
           sigma_sh = this.channelParams.sigmaSH;
           pwr_sh = sigma_sh^2;
           m_prime = this.plPrime( current_point, prev_point);
           m_cond = m_d + exp(-d/beta)*(S - epsilon - m_r);
           var_cond = pwr_sh*(1 - exp(-2*d/beta));
           
           
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

            m_d = cc.getGammaPLdBAtPt(current_point);
            m_l = cc.getGammaPLdBAtPt(prev_point);
            
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

        function m_prime = plPrime(this, current_point, prev_point)
            cp = this.channelParams;
            qb = (cp.qBase - ...
                [this.commChannel.region(2), this.commChannel.region(4)]);
            d2b = norm(current_point - qb); 
            
            v1 = qb - prev_point; v2 = qb - current_point - prev_point;
            cos_theta = -max(min(dot(v1,v2)/(norm(v1)*norm(v2)),1),-1);
            m_prime = -10*cp.nPL/(d2b*log(10))*cos_theta;
        end
        
        function h1_cov = h1Cov(this, d)
            cp = this.channelParams;
            h1_cov = cp.sigmaSH*exp(d/cp.decorrSH);
        end
        
        function h1_cov_prime = h1CovPrime(this, d)
            cp = this.channelParams;
            h1_cov_prime = this.h1Cov(d)/cp.decorrSH;
        end
        
        function h2_cov = h2Cov(this, d)
            cp = this.channelParams;
            h2_cov = cp.sigmaSH*exp(-d/cp.decorrSH);
        end
        
        function h2_cov_prime = h2CovPrime(this, d)
            cp = this.channelParams;
            h2_cov_prime = -this.h2Cov(d)/cp.decorrSH;
        end
        
        function cp = mpCDF(this, x)
           if this.noMP
               cp = (x>=0);
           else
               x_lin = 10.^(x/20);
               
               %find closest entry in LUT
               x_lut_index = max(min(round(x_lin/(this.LUTMax/this.LUTRes)), this.LUTRes), 0) + 1;
               cp = this.riceCDFLUT(x_lut_index);
           end
        end
        
        function computeRiceCDFLUT(this)
            %hard code for now
            this.LUTMax = 7*this.ricDist.std();
            
            %now build out the table
            x = 0:this.LUTMax/this.LUTRes:this.LUTMax;
            
            this.riceCDFLUT = this.ricDist.cdf(x);
        end
        
        function path_dim = getPathDim(this, path)
            path_dim = length(path);
            if min(size(path)) == 1
                path_dim = 1;
            end 
        end
        
    end
    
end

