function total_cost = LITotalEnergyWithPOI(path, bs_cc, qos_params, mp, poi_cc, scenario)
    %based on (5) in "Motion-Communication Co-optimization with
    %Cooperative Load Transferin Mobile Robotics: an Optimal Control Perspective"
    %For variable transmit power, fully observale chanel
    %Motion energy as described in "Human-Robot Collaborative Site 
    %Inspection under Resource Constraints"
    
    v_const = mp.VConst; dist_scale = 1/bs_cc.res;
    
    total_cost = 0;
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
        
        total_cost = total_cost + mp.motionEnergy(dist) + comm_energy;
        
        bs_req_comm_power_a = bs_req_comm_power_b;
        poi_req_comm_power_a = poi_req_comm_power_b;
    end
end
