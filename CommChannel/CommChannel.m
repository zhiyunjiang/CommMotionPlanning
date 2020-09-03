classdef CommChannel < handle
    %ChannelGenerator - Handles simulating/generating a channel
    
    properties (Access = public)
        channelParams;
        nSin;
        region;
        res;
    end
    
    properties(Access = private)
        gx; gy;
        gammaPLdB;
        gammaSHdB;
        gammaMPdB;
        ricDist;
    end
    
    methods
        function this = CommChannel(cp, N_sin, region, res, gamma_PL)
            this.channelParams = cp;
            this.nSin = N_sin;
            this.region = region;
            this.res = res;
            
            % the rectangular region specifying the environment
            x_max = this.region(1);
            x_min = this.region(2);
            y_max = this.region(3);
            y_min = this.region(4);
            
            [this.gx, this.gy] = meshgrid(x_min:1/res:x_max,y_min:1/res:y_max);
            
            if nargin == 4
              %If PL not given, we need to generate
              this.generatePL();
            else
                this.gammaPLdB = gamma_PL;
            end
            
            %initialize shadowing and multipath to 0
            this.gammaSHdB = 0;
            this.gammaMPdB = 0;
            
           K_ric = this.channelParams.kRic;
           this.ricDist = makedist('Rician', 's', sqrt(K_ric/(K_ric + 1)), 'sigma',1/sqrt(2*(1+K_ric)));
        end
        
        function [g_x, g_y] = getMeshGrid(this)
            g_x = this.gx;
            g_y = this.gy;
        end
        
        function gamma_PL_dB = generatePL(this)
            % generating the path loss
            gamma_PL_dB = generate_pathloss(this.region, this.channelParams.qBase,...
                this.res, this.channelParams.kPL, this.channelParams.nPL);
            
            this.gammaPLdB = gamma_PL_dB;
        end
        
        function gamma_SH_dB = simulateSH(this)
            gamma_SH_dB = generate_shadowing(this.channelParams.sigmaSH^2, this.channelParams.decorrSH,...
                this.nSin, this.channelParams.psdAtFC, this.gx, this.gy, this.res);
            
            this.gammaSHdB = gamma_SH_dB;
        end
        
        function gamma_MP_dB = simulateMP(this, mp_model)
            if nargin == 1
                mp_model = 1;%default to Rician
            end
            
            if mp_model == 1%Rician MP
                gamma_MP_lin = generate_multipath(this.channelParams.lambda, this.gx, this.gy, ...
                    this.res, this.channelParams.kRic, this.channelParams.corrMP);
                gamma_MP_dB = 10*log10(gamma_MP_lin);
                this.gammaMPdB = gamma_MP_dB;
            elseif mp_model == 2%log linear model
                %get size
                sz = size(this.getGammaPLdB());
                this.gammaMPdB = this.channelParams.sigmaMP*randn(sz);
            end
            
        end
        
        function gamma_TOT_dB = getGammaTOTdB(this)
            gamma_TOT_dB = this.gammaPLdB + this.gammaSHdB + this.gammaMPdB;
        end
        
        function gamma_PL_dB = getGammaPLdB(this)
            gamma_PL_dB = this.gammaPLdB;
        end
        
        function gamma_SH_dB = getGammaSHdB(this)
            gamma_SH_dB = this.gammaSHdB;
        end
        
        function gamma_MP_dB = getGammaMPdB(this)
            gamma_MP_dB = this.gammaMPdB;
        end
        
        function [sample_pos, sample_vals] = randSample(this, n)
           xs = randi(round(this.region(1)*this.res) ,[n,1]);
           ys = randi(round(this.region(3)*this.res) , [n,1]);
           sample_pos = [xs,ys];
           sample_vals = this.getGammaTOTdBAtPoint(sample_pos);
        end
        
        function gamma_PL_dB_at_point = getGammaPLdBAtPoint(this, point)
           cp = this.channelParams;
           raw_point = this.resPoint2Raw(point);
           gamma_PL_dB_at_point = cp.kPL - 10*cp.nPL*log10(sqrt(sum((cp.qBase - raw_point).^2, 2)));
        end
        
        function q_BRes = getGridBase(this)
            cp = this.channelParams;
            q_BRes = this.rawPoint2Res(cp.qBase);
        end
        
        function gamma_TOT_dB_at_point = getGammaTOTdBAtPoint(this, point)
            
            gamma_TOT_dB = this.getGammaTOTdB();
            %take into account that grids are 0 indexed, but gamma arrays
            %are 1 indexed
            point = point + 1;%shift over by one to match matrix indices
            lin_idx = sub2ind(size(gamma_TOT_dB), point(:,2), point(:,1));
            gamma_TOT_dB_at_point = gamma_TOT_dB(lin_idx);
        end
        
        %INCOMPLETE
        function gamma = simulatePath(this, path, no_mp, is_markov)
           path_dim = length(path);
           
           if min(size(path)) == 1
               path_dim = 1;
           end
           
           cp = this.channelParams;
           sigmaSH = cp.sigmaSH;
           decorr_sh = cp.decorrSH*this.res;
           
           gamma = zeros([path_dim, 1]);
           if ~is_markov
              this.simulateMP();
              this.simulateSH();
              gamma_TOT = this.getGammaTOTdB();
              gamma(1) = gamma_TOT(path(1,2) + 1, path(1,1) + 1);
           else
               gamma_mps = this.simMP(no_mp, path_dim);
               z_shs = sigmaSH*randn([path_dim, 1]);
               gamma_shs = zeros([path_dim,1]);
               gamma_pl_prev = this.getGammaPLdBAtPoint(path(1,:));
               gamma_shs(1) = z_shs(1);
               gamma(1) = gamma_pl_prev + gamma_shs(1)...
               + gamma_mps(1);
           end
           
           for i = 2:path_dim
               if is_markov
                   dist = norm(path(i-1,:) - path(i,:));
                   rho = exp(-dist/decorr_sh);
                   gamma_pl_cur = this.getGammaPLdBAtPoint(path(i,:));
                   gamma_shs(i) = rho*gamma_shs(i-1) + sqrt(1-rho^2)*z_shs(i);
                   gamma(i) = gamma_pl_cur + gamma_shs(i) + gamma_mps(i);
                   
               else
                   %TODO implement for non-markovian paths (will have to
                   %build out covariance matrices, etc...)
                   %for now, brute force
                   gamma(i) = gamma_TOT(path(i,2) + 1, path(i,1) + 1);
               end
           end
        end  
        
        function plot(this, components, plot2d)
            if nargin < 3
                plot2d = 0;
            end
            
            if nargin==1 || components == 0
                gamma = this.getGammaTOTdB();
                component_labels = 'PL + SH + MP';
            elseif components == 1
                %pathloss only
                gamma = this.gammaPLdB;
                component_labels = 'PL';
            elseif components == 2
                %shadowing only
                gamma = this.gammaSHdB;
                component_labels = 'SH';
            elseif components == 3
                %multipath only
                gamma = this.gammaMPdB;
                component_labels = 'MP';
            elseif components == 4
                %pathloss and shadowing
                gamma = this.gammaPLdB + this.gammaSHdB;
                component_labels = 'PL + SH';
            end
            
            f = figure;
            fnt_siz = 16;
            if ~plot2d
                
                surf(this.gx, this.gy, gamma, 'EdgeColor','none');
                % light
                % shading interp
                xlabel('x (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
                ylabel('y (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
                z_label = strcat('Received Power (', component_labels, ') (dBm)');
                zlabel(z_label,'FontSize', fnt_siz ,  'FontWeight','bold');
                axis tight
                grid on
                set(gca, 'FontSize', fnt_siz, 'FontWeight', 'bold');
            else
                imagesc(0,0, gamma);
                
                c = colorbar;
                c.Label.String = 'Received Power (dBm)';
                xlabel('x (m)');
                ylabel('y (m)');
                title(strcat('Received Power (', component_labels, ')'));
            end
            

        end
        
        function plotConnected(this, gamma_th)
            gamma = this.getGammaTOTdB();
            conn = double((gamma >= gamma_th));
            
            imagesc(0,0, conn);
            colormap([0.5 0.5 0.5; 0 0 0.75]);
            colorbar('YTick',[0.25 0.75],'YTicklabel',{'Disconnected', 'Connected'},'FontSize', 7, 'FontName', 'Calibri')
            title(strcat('Connected Regions for \Gamma_{th} = ', sprintf('%d dBm', gamma_th))) ;
            xlabel('x (m)');
            ylabel('y (m)');
        end
        
        function plotRequiredTXPower(this, receiver_noise, BER, R)
            K =  -1.5/log(5*BER);
            gamma = this.getGammaTOTdB();
            CNR_lin = 10.^(gamma/10) / receiver_noise;
            req_power = ((2^R - 1)/K)*(1./CNR_lin);
            imagesc(0, 0, req_power);
            c = colorbar;
            c.Label.String = 'Required TX Power (mW)';
            title(sprintf('Required TX Power, BER = %d', BER));
        end
        
    end
    
    methods (Access = private)
        function gamma_mp = simMP(this, no_mp, n)
            if nargin == 2
                n = 1;
            end
            
            %for now assume uncorrelated MP
            if no_mp
                gamma_mp = zeros([n, 1]);
            else
                gamma_mp_lin = this.ricDist.random([n,1]);
                gamma_mp = 20*log10(gamma_mp_lin);
            end
        end
        
        function [x_index, y_index] = point2Index(this, point)
            x_index = round(this.res*(point(1) - this.region(2)));
            y_index = round(this.res*(point(2) - this.region(4)));
        end
        
        function res_pt = rawPoint2Res(this, point)
            res_pt = this.res*(point - [this.region(2), this.region(4)]);
        end
        
        function raw_pt = resPoint2Raw(this, point)
            raw_pt = (point/this.res) + [this.region(2), this.region(4)];
        end
    end
end

