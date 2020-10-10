classdef MotionParams < handle
    
    properties (SetAccess = private)
        K1;
        K2;
        VConst;
        MotionPower;
    end
    
    methods
        function this = MotionParams(k1, k2, v_const)
            this.K1 = k1;
            this.K2 = k2;
            this.VConst = v_const;
            this.MotionPower = (k1 + (k2/v_const));
        end
        
        function energy_J = motionEnergy(this, dist)
           energy_J = this.MotionPower*dist;
        end
        
    end
end

