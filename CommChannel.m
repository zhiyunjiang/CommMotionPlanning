classdef CommChannel < handle
    %ChannelGenerator - Handles simulating/generating a channel
    
    properties (Access = public)
        channelParams;
        nSin;
        region;
        res
    end
    
    properties(Access = private)
        gx; gy;
        gammaPLdB;
        gammaSHdB;
        gammaMPdB;
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
        
        function gamma_SH_dB = simulateShadowing(this)
            gamma_SH_dB = generate_shadowing(this.channelParams.sigmaSH^2, this.channelParams.decorrSH,...
                this.nSin, this.channelParams.psdAtFC, this.gx, this.gy, this.res);
            
            this.gammaSHdB = gamma_SH_dB;
        end
        
        function gamma_MP_dB = simulateMP(this)
            gamma_MP_lin = generate_multipath(this.channelParams.lambda, this.gx, this.gy, ...
                this.res, this.channelParams.kRic, this.channelParams.corrMP);
            gamma_MP_dB = 10*log10(gamma_MP_lin);
            this.gammaMPdB = gamma_MP_dB;
        end
        
        function gamma_TOT_dB = getGammaTOTdB(this)
            gamma_TOT_dB = this.gammaPLdB + this.gammaSHdB + this.gammaMPdB;
        end
        
        function plot(this, components)
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
            y_label = strcat('Received Power (', component_labels, ') (dBm)');
            f = figure;
            fnt_siz = 16;
            
            surf(this.gx, this.gy, gamma, 'EdgeColor','none');
            % light
            % shading interp
            xlabel('x (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
            ylabel('y (m)', 'FontSize', fnt_siz,  'FontWeight', 'bold');
            zlabel(y_label,'FontSize', fnt_siz ,  'FontWeight','bold');
            axis tight
            grid on
            set(gca, 'FontSize', fnt_siz, 'FontWeight', 'bold');

            maximize(f)

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
        
        function gamma_PL_dB_at_point = getGammaPLdBAtPoint(this, point)
            %not sure why y and x indices end up flipped in gammaPLdB, but
            %they do, and since it's the inherited code that works that
            %way, will leave untouhced for now
            
            %take into account that grids are 0 indexed, but gamma arrays
            %are 1 indexed
            if sum(mod(point,1)) == 0 %round, whole numbers
                gamma_PL_dB_at_point = this.gammaPLdB(point(2)+1, point(1)+1);
            else
               %compute
               cp = this.channelParams;
               raw_point = this.resPoint2Raw(point);
               gamma_PL_dB_at_point = cp.kPL - 10*cp.nPL*log10(norm(cp.qBase - raw_point));
            end
        end
        
        function q_BRes = getGridBase(this)
            cp = this.channelParams;
            q_BRes = this.rawPoint2Res(cp.qBase);
        end
        
        function gamma_TOT_dB_at_point = getGammaTOTdBAtPoint(this, point)
            
            gamma_TOT_dB = this.getGammaTOTdB();
            %take into account that grids are 0 indexed, but gamma arrays
            %are 1 indexed
            gamma_TOT_dB_at_point = gamma_TOT_dB(point(2) + 1, point(1) + 1);
        end
        
        
    end
    
    methods (Access = private)
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

