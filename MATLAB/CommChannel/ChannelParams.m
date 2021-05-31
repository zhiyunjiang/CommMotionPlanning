%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predicted Channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for the probabilistic modeling of channel gain. Characterizes
% the path loss, shadowing, and multipath componenets 

classdef ChannelParams
    %ChannelParams - Contains parameters defining the communication channel
    properties
        qBase;%base station location in [x,y], generally assumed to be meters
        nPL;% path loss exponent
        kPL;% path loss constant
        sigmaSH;% shadowing standard deviation
        decorrSH;% shadowing decorrelation distance. Referred to beta in some papers
        decorrMP;% multipath decorrelation distance.
        lambda;% wavelength
        kRic;%  if modeling multipath as rician, kRic parameterizes the distribution
        corrMP% if corrMP = 1, treating multipath affects as correlated. Else treat as uncorrelated
        psdAtFC%
        sigmaMP;%if modeling shadowing as a zero-mean, log normal distribution (uncorrelated)
    end
    
    methods
        function this = ChannelParams(q_b, n_PL, K_PL, sigma_SH, sh_decorr,...
                                        lambda, K_ric, mp_decorr, corr_mp, PSD_at_f_c, sigma_mp)
           % Position of the base station (remote station or transmitter)
            this.qBase = q_b;

            % Path loss parameters
            % The resulting path loss component in dB is 
            % gamma_PL_dB = K_PL - 10 * n_PL * log10(d), where d is the distance to the base station.
            % K_PL is the path loss constant in dB and n_PL is the path loss exponent.
            this.nPL = n_PL;
            if this.nPL<0
                warning('Path loss exponent less than 0, which is unusual');
            end
            this.kPL = K_PL;


            % Shadowing parameters
            % alpha : power of the shadowing (in dB)
            % beta : decorrelation distance (in meter). The shadowing correlation
            %        model we use is: alpha*exp(-distance/beta).

            this.sigmaSH = sigma_SH;            
            this.decorrSH = sh_decorr;

            % Multipath fading parameters
            % lambda: the wavelength of transmission (in meter)
            % res = the resolution of the grid (in samples/m) = 1/(dis between 2 samples along x or y axis)
            % K_ric = parameter of the Rician distribution (see Goldsmith page 79)
            this.kRic = K_ric;
            
            % lambda = 0.125 corresponds to a carrier frequency of 2.4 GHz.
            this.lambda = lambda;

            % ss_decorr is the decorrelation distance for multipath fading.
            % This is the point where J_0(2*pi*ss_decorr/lambda) is equal to 0 (see Eq. 3.26 of Goldsmith).
            % Note that ss_decorr = 0.4*lambda.
            % This variable is not an input to the channel simulator, but is used to 
            % to provide a guideline for choosing the simulation resoltuion.
            this.decorrMP = mp_decorr;         
            
            this.corrMP = corr_mp;
            
            if nargin == 10
                sigma_mp = 1;
            end
            this.sigmaMP = sigma_mp;
            
            this.psdAtFC = PSD_at_f_c;
        end 
    end
end

