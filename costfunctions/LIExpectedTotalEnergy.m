function total_cost = LIExpectedTotalEnergy(n1, n2, path, cawo, receiver_noise, R, K)
    %based on (5) in "Motion-Communication Co-optimization with
    %Cooperative Load Transferin Mobile Robotics: an Optimal Control Perspective"
    %For variable transmit power, partially observale chanel
    
    %also using motion energy as described in "Human-Robot Collaborative Site 
    %Inspection under Resource Constraints"
    k_1 = 7.4; k_2 = 0.29; v_const = 1*cawo.cc.res;
    
    total_cost = 0;
    req_comm_power_a = reqPower(cawo, receiver_noise, R, K, path(1,:));
    for i=2:length(path)
        dist = norm(path(i-1,:) - path(i,:));
        req_comm_power_b = reqPower(cawo, receiver_noise, R, K, path(i,:));
        motion_energy = (k_1 + (k_2/v_const))*dist;
        %scale by 1/1000 to convert to Joules
        comm_energy = (dist/v_const)*(req_comm_power_a + req_comm_power_b)/(2*1000);
        total_cost = total_cost + motion_energy + comm_energy;
        req_comm_power_a = req_comm_power_b;
    end
end

function req_power = reqPower(cawo, receiver_noise, R, K, pos)
    channel_power_dBm = cawo.posteriorExpecteddB(pos);
    CNR_lin = 10.^(channel_power_dBm/10) / receiver_noise;
    req_power = ((2^R - 1)/K)*(1/CNR_lin);
end

