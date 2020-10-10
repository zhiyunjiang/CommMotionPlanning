classdef QoSParams < handle
    
    properties (SetAccess = private)
        BER;
        R;
        RXNoise;%in mW;
        K;
    end
    
    methods
        function this = QoSParams(ber, r, rx_noise)
            this.BER = ber;
            this.R = r;
            this.RXNoise = rx_noise;
            this.K = -1.5/log(5*this.BER);
        end
        
        function req_power_W = reqTXPower(this, channel_power_dBm)
            CNR_lin = 10.^(channel_power_dBm/10) / this.RXNoise;
            %divide by 1000, else we have mW
            req_power_W = (((2^this.R - 1)/this.K)*(1./CNR_lin))/1000;
        end
    end
end

