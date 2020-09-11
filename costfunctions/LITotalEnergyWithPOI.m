function total_cost = LITotalEnergyWithPOI(n1, n2, path, bs_cc, receiver_noise, R, K, poi_cc)
    %based on (5) in "Motion-Communication Co-optimization with
    %Cooperative Load Transferin Mobile Robotics: an Optimal Control Perspective"
    %For variable transmit power, partially observale chanel
    
    %also using motion energy as described in "Human-Robot Collaborative Site 
    %Inspection under Resource Constraints"
    k_1 = 7.4; k_2 = 0.29; v_const = 1*bs_cc.res;
    
    total_cost = 0;
    bs_req_comm_power_a = reqPower(bs_cc, receiver_noise, R, K, path(1,:));
    poi_req_comm_power_a = reqPower(poi_cc, receiver_noise, R, K, path(1,:));
    for i=2:length(path)
        dist = norm(path(i-1,:) - path(i,:));
        bs_req_comm_power_b = reqPower(bs_cc, receiver_noise, R, K, path(i,:));
        poi_req_comm_power_b = reqPower(poi_cc, receiver_noise, R, K, path(i,:));
        motion_energy = (k_1 + (k_2/v_const))*dist;
        %scale by 1/1000 to convert to Joules
        bs_comm_energy = (dist/v_const)*(bs_req_comm_power_a + bs_req_comm_power_b)/(2*1000);
        pos_comm_energy = (dist/v_const)*(poi_req_comm_power_a + poi_req_comm_power_b)/(2*1000);
        total_cost = total_cost + motion_energy + bs_comm_energy + pos_comm_energy;
        bs_req_comm_power_a = bs_req_comm_power_b;
        poi_req_comm_power_a = poi_req_comm_power_b;
    end
end

function req_power = reqPower(cc, receiver_noise, R, K, pos)
    channel_power_dBm = cc.getGammaTOTdBAtPoint(pos);
    CNR_lin = 10.^(channel_power_dBm/10) / receiver_noise;
    req_power = ((2^R - 1)/K)*(1/CNR_lin);
end
