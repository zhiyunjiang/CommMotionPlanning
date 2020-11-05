%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QoSParams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contains parameters relevant to fining required transmit power for a
% given BER and receiver noise

classdef QoSParams < handle
    
    properties (SetAccess = private)
        BER;% target bir error rate
        R;% spectral efficiency
        RXNoise;%in mW;
        K;
    end
    
    methods (Access = public)
        function this = QoSParams(ber, r, rx_noise)
            if ber < 0 || ber >1
               error('QosParams: BER must be in range [0,1]'); 
            end
            this.BER = ber;
            if r<0
                error('QosParams: spectral efficiency must be non-negative');
            end
            this.R = r;
            if rx_noise < 0
                error('QoSParams: receiver noise power must be non-negative');
            end
            this.RXNoise = rx_noise;
            this.K = -1.5/log(5*this.BER);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reqTXPower
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculates the TX power required to meet the QoS requirement
        % given channel gain in db
        %
        % Input:
        % this - reference to the AStarSolver object
        % channel_gain_dB - channel gain in dB
        %
        % Output:
        % req_power_W - the required transmit power in Watts
        function req_power_W = reqTXPower(this, channel_gain_dB)
            CNR_lin = this.toLinCNR(channel_gain_dB);
            req_power_W = this.calcPower(1./CNR_lin);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % expReqTXPowerW
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculates the expected TX power required to meet the QoS 
        % requirement given channel gain in db
        %
        % Input:
        % this - reference to the AStarSolver object
        % exp_ch_gain_dB - expected channel gain in dB
        % exp_ch_gain_var - channel gain variance
        %
        % Output:
        % req_power_W - the required transmit power in Watts
        function exp_req_power_W = expReqTXPowerW(this, exp_ch_gain_dB, exp_ch_gain_var)
            exp_CNR_lin = this.toLinCNR(exp_ch_gain_dB);
            
            exp_CNR_lin_inverse = exp( (log(10)/10)^2 * exp_ch_gain_var/2)./exp_CNR_lin;
            
            exp_req_power_W = this.calcPower(exp_CNR_lin_inverse);
        end
            
    end
    
    methods (Access = private)
        function CNR_lin = toLinCNR(this, channel_power_dBm)
            CNR_lin = 10.^(channel_power_dBm/10) / this.RXNoise;
        end
        
        function power_W = calcPower(this, CNR_lin_inverse)
            %divide by 1000 to get Watts, not milliwatts
            power_W = (((2^this.R - 1)/this.K)*CNR_lin_inverse)/1000;
        end
    end
end

