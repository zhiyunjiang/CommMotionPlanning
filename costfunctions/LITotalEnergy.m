function total_cost = LITotalEnergy(path, cc, qos, mp)
    %based on (5) in "Motion-Communication Co-optimization with
    %Cooperative Load Transferin Mobile Robotics: an Optimal Control Perspective"
    %For variable transmit power, fully observable channel
    %also using motion energy as described in "Human-Robot Collaborative Site 
    %Inspection under Resource Constraints"
    
    v_const = mp.VConst; dist_scale = 1/cc.res;

    total_cost = 0;
    req_comm_power_a = qos.reqTXPower(cc.getGammaTOTdBAtPoint(path(1,:)));
    for i=2:length(path)
        dist = dist_scale*norm(path(i-1,:) - path(i,:));
        
        req_comm_power_b = qos.reqTXPower(cc.getGammaTOTdBAtPoint(path(i,:)));
        comm_energy = (dist/v_const)*(req_comm_power_a + req_comm_power_b)/2;
        
        total_cost = total_cost + mp.motionEnergy(dist) + comm_energy;
        req_comm_power_a = req_comm_power_b;
    end
end


