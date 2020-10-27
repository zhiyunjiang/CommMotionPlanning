classdef QoSParams < handle
    
    properties (SetAccess = private)
        BER;
        R;
        RXNoise;%in mW;
        K;
    end
    
    methods (Access = public)
        function this = QoSParams(ber, r, rx_noise)
            this.BER = ber;
            this.R = r;
            this.RXNoise = rx_noise;
            this.K = -1.5/log(5*this.BER);
        end
        
        function req_power_W = reqTXPower(this, channel_power_dBm)
            CNR_lin = this.toLinCNR(channel_power_dBm);
            req_power_W = this.calcPower(1./CNR_lin);
        end
        
        function exp_req_power_W = expReqTXPowerW(this, exp_ch_pwr_dBm, exp_ch_pwr_var)
            exp_CNR_lin = this.toLinCNR(exp_ch_pwr_dBm);
            
            exp_CNR_lin_inverse = exp( (log(10)/10)^2 * exp_ch_pwr_var/2)./exp_CNR_lin;
            
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

