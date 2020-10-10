function total_cost = LIExpectedTotalEnergy(path, cawo, qos, mp)
    %based on (5) in "Motion-Communication Co-optimization with
    %Cooperative Load Transferin Mobile Robotics: an Optimal Control Perspective"
    %For variable transmit power, partially observale chanel
    %Motion energy as described in "Human-Robot Collaborative Site 
    %Inspection under Resource Constraints"
    
    v_const = mp.VConst; dist_scale = 1/bs_cc.res;
    
    total_cost = 0;
    req_comm_power_a = reqPower(cawo, receiver_noise_lin, R, K, path(1,:));
    for i=2:length(path)
        dist = dist_scale*norm(path(i-1,:) - path(i,:));%convert distance to meters
        req_comm_power_b = qos.reqTXPower(cawo.posteriorExpecteddB(path(i,:)));

        comm_energy = (dist/v_const)*(req_comm_power_a + req_comm_power_b)/2;
        total_cost = total_cost + mp.motionEnergy(dist) + comm_energy;
        req_comm_power_a = req_comm_power_b;
    end
end
