%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LITotalEnergyWithPOI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear approximation of the line integral which gives the total
% energy consumed over a path (both motion and communication). We assume
% a linear interpolation of power between two  waypoints. With multiple
% remote stations (poi), the rqeuired communication power depends on the
% robots task (scenario).
%
% Inputs:
% path - array of path waypoints. Path(i) = [xi, yi].
% bs_cc - communication channel with the base station.
% poi_cc - communication channel with the other remote station.
% qp - quality of service parameters and settings
% mp - motion parameters
% scenario - the robots tak. 1 - sensing/surveillance. 2 - broadcasting.
%                            3 - relaying
% 
% Outputs:
% double energy_J - approximate energy consumed over path in Joules
function energy_J = LITotalEnergyWithPOI(path, bs_cc, poi_cc, qos_params, mp, scenario)
    %based on (5) in "Motion-Communication Co-optimization with
    %Cooperative Load Transferin Mobile Robotics: an Optimal Control Perspective"
    %For variable transmit power, fully observale chanel
    %Motion energy as described in "Human-Robot Collaborative Site 
    %Inspection under Resource Constraints"
    
    v_const = mp.VConst; dist_scale = 1/bs_cc.res;
    
    energy_J = 0;
    bs_req_comm_power_a = qos_params.reqTXPower(bs_cc.getGammaTOTdBAtPoint(path(1,:)));
    poi_req_comm_power_a = qos_params.reqTXPower(poi_cc.getGammaTOTdBAtPoint(path(1,:)));
    for i=2:length(path)
        dist = dist_scale*norm(path(i-1,:) - path(i,:));%convert distance to meters
        %compute power requred to communicate with the two points
        bs_req_comm_power_b = qos_params.reqTXPower(bs_cc.getGammaTOTdBAtPoint(path(i,:)));
        poi_req_comm_power_b = qos_params.reqTXPower(poi_cc.getGammaTOTdBAtPoint(path(i,:)));
        
        %now compute the required comm energy, depending on which scenario
        %we're in
        bs_comm_energy = (dist/v_const)*(bs_req_comm_power_a + bs_req_comm_power_b)/2;
        pos_comm_energy = (dist/v_const)*(poi_req_comm_power_a + poi_req_comm_power_b)/2;
        if scenario == 1
            %find the min! Only have to communicate with one base station
            comm_energy = min([bs_comm_energy, pos_comm_energy]);
        elseif scenario == 2
            %find the max! Must communicate with all base stations
            comm_energy = max([bs_comm_energy, pos_comm_energy]);
        elseif scenario == 3
            %Sum - simultaneously communicate with both
            comm_energy = bs_comm_energy + pos_comm_energy;
        else
            error('Scenario must take on value 1,2, or 3');
        end
        
        energy_J = energy_J + mp.motionEnergy(dist) + comm_energy;
        
        bs_req_comm_power_a = bs_req_comm_power_b;
        poi_req_comm_power_a = poi_req_comm_power_b;
    end
end
