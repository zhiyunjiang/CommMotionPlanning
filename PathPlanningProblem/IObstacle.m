%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IObstacle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obstacle interface. Allows us to create arbitrary obstacles

classdef IObstacle < handle & matlab.mixin.Heterogeneous
    %IObstacle - Generic Obstacle interface
    
    %The only need to be able to check if an obstale obstructs the a line
    %segement with end points a and b
    methods(Abstract)
        obstructs(obj, a, b)
        plotObstacle(obj, scale, offset)
    end
    
    methods (Access = protected, Static = true)
        function defaultObject = getDefaultScalarElement()
            defaultObject = RectObs([0,0,0,0]);
        end
    end
end

