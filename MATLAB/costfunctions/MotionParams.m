%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MotionParams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters describing robots power and energy conusmption from movement.
% Based on 
% Y. Mei, Y. Lu, Y. Hu, and C. Lee.  Deployment of mobile robots withenergy
% and timing constraints.IEEE Transactions on Robotics, 2006.


classdef MotionParams < handle
    % Properties (private):
    % K1 - first motion power constant
    % K2 - second motion power constant
    % VConst - constant velocity in m/s (we assume constant velocity throughout
    %           trajectory
    % MotionEnergyScalar - Proportionality constant for energy and distance,
    %                       i.e. energy = MotionEnergyScalar*distance
    properties (SetAccess = private)
        K1;
        K2;
        VConst;
        MotionEnergyScalar;
    end
    
    % Methods (public):
    % (Constructor) MotionParams - creates new MotionParams object
    % motionEnergy - calculates the motion energy expended over a distance
    methods
        
        function this = MotionParams(k1, k2, v_const)
            this.K1 = k1;
            this.K2 = k2;
            this.VConst = v_const;
            this.MotionEnergyScalar = (k1 + (k2/v_const));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % motionEnergy
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Calculates motion energy consumed over a distance
        % Input:
        % this - reference to the MotionParams object
        % dist - distance traveled
        %
        % Ouput:
        % energy_J - the energy in Joules used to travel the given distance
        function energy_J = motionEnergy(this, dist)
           energy_J = this.MotionEnergyScalar*dist;
        end
        
    end
end

