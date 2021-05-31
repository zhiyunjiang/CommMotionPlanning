classdef TXPwr
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dBm;
        dBW;
        mW;
        W;
    end
    
    methods (Access = public)
        function this = TXPwr(dBm)
            this.dBm = dBm;
            this.dBW = TXPwr.dBm2dBW(dBm); 
            this.mW = TXPwr.dBm2mW(dBm);
            this.W = TXPwr.dBm2W(dBm);
        end
        
        function view(this)
           fprintf("%f dBm\t%f dBW\t%f mW\t%f W\n",this.dBm, this.dBW, this.mW, this.W);
        end
    
    end
    
    methods (Static, Access = public)
        function dBW = dBm2dBW(dBm)
            dBW = dBm - 30;
        end
        
        function mW = dBm2mW(dBm)
           mW = 10^(dBm/10); 
        end
        
        function W = dBm2W(dBm)
            W = TXPwr.dBm2mW(dBm)/1000;
        end
    end
    
end

