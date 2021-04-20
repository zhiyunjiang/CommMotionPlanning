%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CommChannel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Object representing a simulated communication channel.

classdef CommChannel < handle
    %ChannelGenerator - Handles simulating/generating a channel
    
    properties (Access = public)
        cp;%ChannelParams objects
        nSin;%Simulation channel. See ChannelSims folder for details
        region;% [x max, x min, y max, y min] in meters
        res;% resolution (i.e. 1/step)
        gx; gy;
    end
    
    properties(Access = private)
        gammaPLdB;% as generated by ChannelSim package, first index corresponds to y value, second to x
        gammaSHdB;% as generated by ChannelSim package, first index corresponds to y value, second to x
        gammaMPdB;% as generated by ChannelSim package, first index corresponds to y value, second to x
        ricDist;
        xImax;
        yImax;
    end
    
    methods
        function this = CommChannel(cp, N_sin, region, res, gamma_PL)
            this.cp = cp;
            this.nSin = N_sin;
            this.region = region;
            this.res = res;
            
            % the rectangular region specifying the environment
            x_max = this.region(1);
            x_min = this.region(2);
            y_max = this.region(3);
            y_min = this.region(4);
            if x_max <= x_min || y_max <= y_min
               error('Region''s max values must be greater than its min values'); 
            end
            
            [this.gx, this.gy] = meshgrid(x_min:1/res:x_max,y_min:1/res:y_max);
            [this.yImax, this.xImax] = size(this.gy); 
            
            if nargin == 4
              %If PL not given, we need to generate
              this.generatePL();
            else
                this.gammaPLdB = gamma_PL;
            end
            
            %initialize shadowing and multipath to 0
            this.gammaSHdB = 0;
            this.gammaMPdB = 0;
            
            K_ric = this.cp.kRic;
            this.ricDist = makedist('Rician', 's', sqrt(K_ric/(K_ric + 1)), 'sigma',1/sqrt(2*(1+K_ric)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % generatePL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate the PL component of the channel based on the channel
        % params
        % Input:
        % this - reference to the CommChannel object
        %
        % Output:
        % gamma_PL_dB - matrix of PL powers (in dB) over the region
        function gamma_PL_dB = generatePL(this)
            % generating the path loss
            gamma_PL_dB = generate_pathloss(this.region, this.cp.qBase,...
                this.res, this.cp.kPL, this.cp.nPL);
            
            this.gammaPLdB = gamma_PL_dB;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % simulateSH
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate a realization of the channel's shadowing component over
        % the entire region
        % Input:
        % this - reference to the CommChannel object
        %
        % Output:
        % gamma_SH_dB - matrix of SH powers (in dB) over the region
        function gamma_SH_dB = simulateSH(this)
            gamma_SH_dB = generate_shadowing(this.cp.sigmaSH^2, this.cp.decorrSH,...
                this.nSin, this.cp.psdAtFC, this.gx, this.gy, this.res);
            
            this.gammaSHdB = gamma_SH_dB;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % simulateMP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generate a realization of the channel's multipath fading 
        % component over the entire region
        % Input:
        % this - reference to the CommChannel object
        % mp_model - flag indicating which multipath model to use. 1
        %               indicates Rician (default). 2 indicates log normal.
        %
        % Output:
        % gamma_MP_dB - matrix of MP power (in dB) over the region
        function gamma_MP_dB = simulateMP(this, mp_model)
            if nargin == 1
                mp_model = 1;%default to Rician
            end
            
            if mp_model == 1%Rician MP
                gamma_MP_lin = generate_multipath(this.cp.lambda, this.gx, this.gy, ...
                    this.res, this.cp.kRic, this.cp.corrMP);
                gamma_MP_dB = 10*log10(gamma_MP_lin);
                this.gammaMPdB = gamma_MP_dB;
            elseif mp_model == 2%log normal model
                %get size
                sz = size(this.getGammaPLdB());
                this.gammaMPdB = this.cp.sigmaMP*randn(sz);
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % getGammaTOTdB
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Retreives the total power across the region. Always includes PL.
        % Includes most recent SH, MP components generated. If none have
        % been generated, assumed to be 0.
        % Input:
        % this - reference to the CommChannel object
        %
        % Output:
        % gamma_TOT_dB - matrix of channel power (in dB) over the region
        function gamma_TOT_dB = getGammaTOTdB(this)
            gamma_TOT_dB = this.gammaPLdB + this.gammaSHdB + this.gammaMPdB;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % setKPl
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sets a new value for the path loss constant and updates the
        % path loss power across the region.
        % Input:
        % this - reference to the CommChannel object
        % kpl - the new path loss exponenet value
        function setKPl(this,kpl)
           diff = kpl - this.cp.kPL; 
           this.gammaPLdB = this.gammaPLdB + diff;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % randSample
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sample the channel
        % Input:
        % this - reference to the CommChannel object
        % n - number of samples
        %
        % Output:
        % sample_pos = n X 2 matrix with the x,y location of the samples in
        %               continuous coordinates
        function [sample_pos, sample_vals] = randSample(this, n)
           % sample from the discretized space
           xs = randi( this.xImax, [n,1] );
           ys = randi( this.yImax, [n,1] );
           sample_pos = toRawFromGrid(this.region, this.res, [xs,ys]);
           sample_vals = this.getGammaTOTdBAtPt(sample_pos);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % getGammaPLdBAtPt
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fetches the path loss component power in dB at a point. 
        % Input:
        % this - reference to the CommChannel object
        % pt - point in coordinates in continuous space
        %
        % Output:
        % gamma_PL_dB_at_point = path loss power at the given point
        function gamma_PL_dB_at_point = getGammaPLdBAtPt(this, pt)
           c_p = this.cp;
           gamma_PL_dB_at_point = c_p.kPL - 10*c_p.nPL*log10(sqrt(sum((c_p.qBase - pt).^2, 2)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % getGammaTOTdBAtPt
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Fetches the path loss component power in dB at a point. Snaps
        % point to nearest location on the discretized, noramlize
        % workspace.
        % Input:
        % this - reference to the CommChannel object
        % pt - point in coordinates in continuous space
        %
        % Output:
        % gamma_PL_dB_at_point = path loss power at the given point
        function gamma_TOT_dB_at_point = getGammaTOTdBAtPt(this, pt)
            grid_pt = toGridFromRaw(this.region, this.res, pt);
            gamma_TOT_dB_at_point = this.getGammaTOTdBAtGridPoint( grid_pt);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % getReqTXPowerW
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculates the required TX power at each point in the workspace
        % for a given BER, spectral efficiency, and noise power.
        %
        % Input:
        % this - reference to the CommChannel object
        % qp - QoSParams object containiong BER, spectral efficiency, and
        %       noise power
        % get_linear - if 0, returns power in dBW, else returns in Watts
        %
        % Output:
        % req_power - matrix with entires equal to requried power
        function req_power = getReqTXPowerW(this, qp, get_linear)
            if nargin == 2
                get_linear = 0;
            end
            
            gamma = this.getGammaTOTdB();
            req_power = qp.reqTXPower(gamma);
            
            if ~get_linear
                req_power = 10*log10(req_power);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % getConnectionField
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Computes binary field of connectivity for a given channel power
        % threshold.
        %
        % Input:
        % this - reference to the CommChannel object
        % gamma_th - minimum channel power required for communication (dB)
        %
        % Output:
        % connection_field - matrix with entires equal to 1 if reliable
        %                       communication possible at that point, 0
        %                       otherwise.
        function connection_field = getConnectionField(this, gamma_th)
            gamma = this.getGammaTOTdB();
            connection_field = double((gamma >= gamma_th));
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % simulatePath
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculate the channel power along a path with a new realization
        % of the channel
        %
        % Input:
        % this - reference to the CommChannel object
        % path - nX2 matrix with each row representing a waypoint along the
        %           path
        % no_mp - if set to 
        % is_markov - if true, treat the channel along the path as
        %               Markovian. This allows for faster simulation, but 
        %               is a poorer model of reality
        %
        % Output:
        % gamma - n-dimensional vector with the channel power at each point
        %           along the path.
        function gamma = simulatePath(this, path, no_mp, is_markov)
           path_dim = length(path);
           
           if min(size(path)) == 1
               path_dim = 1;
           end
           
           sigmaSH = this.cp.sigmaSH;
           decorr_sh = this.cp.decorrSH;
           
           gamma = zeros([path_dim, 1]);
           if ~is_markov
               if ~no_mp
                   this.simulateMP();
               else
                    this.gammaMPdB = 0;
               end
              
              this.simulateSH();
              gamma(1) = this.getGammaTOTdBAtPt(path(1,:));
           else%markovian setup
               gamma_mps = this.simMP(no_mp, path_dim);
               z_shs = sigmaSH*randn([path_dim, 1]);
               gamma_shs = zeros([path_dim,1]);
               gamma_pl_prev = this.getGammaPLdBAtPt(path(1,:));
               gamma_shs(1) = z_shs(1);
               gamma(1) = gamma_pl_prev + gamma_shs(1)...
               + gamma_mps(1);
           end
           
           for i = 2:path_dim
               if is_markov
                   dist = norm(path(i-1,:) - path(i,:));
                   rho = exp(-dist/decorr_sh);
                   gamma_pl_cur = this.getGammaPLdBAtPt(path(i,:));
                   gamma_shs(i) = rho*gamma_shs(i-1) + sqrt(1-rho^2)*z_shs(i);
                   gamma(i) = gamma_pl_cur + gamma_shs(i) + gamma_mps(i);
                   
               else
                   gamma(i) = this.getGammaTOTdBAtPt(path(i,:));
               end
           end
        end  
        
        function grid_pt = toGridCoordinate(this, pt)
            grid_pt = toGridFromRaw(this.region, this.res, pt);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % plotting
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
                
                this.plotField(gamma, strcat('Received Power (', component_labels, ')'), ...
                    'Received Power (dBm)');
            end
            

        end
        
        function plotConnected(this, gamma_th)
            connection_field = this.getConnectionField(gamma_th);
            this.plotField(connection_field, strcat('Connected Regions for \Gamma_{th} = ', sprintf('%d dBm', gamma_th)));
        end
        
        function plotRequiredTXPower(this,  qp, use_linear)
            if nargin == 2
                use_linear = 0;
            end
            
            if use_linear
                c_bar_label = 'Required TX Power (W)';
            else
                c_bar_label = 'Required TX Power (dB)';
            end
            req_power = this.getReqTXPowerW(qp, use_linear);
            this.plotField(req_power,'Required TX Power',c_bar_label);  
        end
        
        function plotDoubleConnectField(this, field, field_title)
            %imagesc is a weird one. Tranpose field before plotting
            imagesc(this.gx(1,:), this.gy(:,1), field');
            set(gca, 'YDir', 'normal');
            title(field_title);
            fnt_size = 12;
            xlabel('x (m)', 'FontSize', fnt_size);
            
            ylabel('y (m)', 'FontSize', fnt_size);
            
            no_conn = [0.5 0.5 0.5]; a_conn =[0 0.5 0]; b_conn = [0 0 0.5]; both_conn = [1 1 1];
            colormap([no_conn; a_conn; b_conn; both_conn]);
            colorbar_sec_size = 3/4;
            first = colorbar_sec_size/2; second = first + colorbar_sec_size;
            third = second + colorbar_sec_size; fourth = third + colorbar_sec_size;
            
            colorbar('YTick',[first, second, third, fourth],'YTicklabel',...
                {'Disconnected', 'A Only', 'B Only' 'Both'},...
                'FontSize', 9, 'FontName', 'Calibri');
        end
        
        function plotField(this, field, field_title, c_bar_label)
            imagesc(this.gx(1,:), this.gy(:,1), field');
            set(gca, 'YDir', 'normal');
            title(field_title);
            xlabel('x (m)');
            ylabel('y (m)');
            
            c = colorbar;
            if nargin==4
                c.Label.String = c_bar_label;
            else
                colormap([0.5 0.5 0.5; 0 0 0.75]);
                colorbar('YTick',[0.25 0.75],'YTicklabel',{'Disconnected', 'Connected'},'FontSize', 7, 'FontName', 'Calibri')
            end
        end
        
        function plotField3D(this, field, field_title, c_bar_label)
            
            surf(this.gx, this.gy, field, 'EdgeColor','none');
            title(field_title)
            % light
            % shading interp
            xlabel('x (m)');
            ylabel('y (m)');
            zlabel(c_bar_label);
            axis tight
            grid on
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % getters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [g_x, g_y] = getMeshGrid(this)
            g_x = this.gx;
            g_y = this.gy;
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
        
        function q_BRes = getGridBase(this)
            q_BRes = toGridFromRaw(rhis.region, this.res, this.cp.qBase);
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
        
        function gamma_TOT_dB_at_point = getGammaTOTdBAtGridPoint(this, grid_pt)
            gamma_TOT_dB = this.getGammaTOTdB();
            lin_idx = sub2ind(size(gamma_TOT_dB), grid_pt(:,2), grid_pt(:,1));
            gamma_TOT_dB_at_point = gamma_TOT_dB(lin_idx);
        end
        
    end
end

