%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LICommEnergy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear approximation of the line integral which gives the communication
% energy consumed over a path. We assume a linear interpolation of power
% between the two given points. We also assume the way points' locations 
% are measured in meters and we move at a velocity of 1m/s. Otherwise, 
% this is proportional but not equal to energy.
%
% Inputs:
% path - array of path waypoints. Path(i) = [xi, yi].
% cc - comm channel object
% qp - quality of service parameters and settings
%
% Outputs:
% double energy_J - approximate comm energy consumed over path in Joules

function energy_J = LICommEnergy(path, cc, qp)
    %based on (5) in "Motion-Communication Co-optimization with
    %Cooperative Load Transferin Mobile Robotics: an Optimal Control Perspective"
    %For variable transmit power, fully observable channel
    
    energy_J = 0;
    req_power_a = qp.reqTXPower( cc.getGammaTOTdBAtPoint( path(1,:) ) );
    for i=2:length(path)
        dist = norm(path(i-1,:) - path(i,:));
        req_power_b = qp.reqTXPower( cc.getGammaTOTdBAtPoint( path(1,:) ) );
        %scale by 1/1000 to get Joules
        energy_J = energy_J + dist*(req_power_a + req_power_b)/(2*1000);
        req_power_a = req_power_b;
    end
end

