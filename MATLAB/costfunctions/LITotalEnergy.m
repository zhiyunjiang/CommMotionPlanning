function total_cost = LITotalEnergy(path, cc, qos, mp, delta)
    %based on (5) in "Motion-Communication Co-optimization with
    %Cooperative Load Transferin Mobile Robotics: an Optimal Control Perspective"
    %For variable transmit power, fully observable channel
    %also using motion energy as described in "Human-Robot Collaborative Site 
    %Inspection under Resource Constraints"
    
    if nargin < 5
        delta = 1;
    end
    
    v_const = mp.VConst;

    total_cost = 0;
    req_comm_power_a = qos.reqTXPower(cc.getGammaTOTdBAtPt(path(1,:)));
    for i=2:length(path)
        dist = norm(path(i-1,:) - path(i,:));
        
        req_comm_power_b = qos.reqTXPower(cc.getGammaTOTdBAtPt(path(i,:)));
        comm_energy = (dist/v_const)*(req_comm_power_a + req_comm_power_b)/2;
        
        total_cost = total_cost + delta*mp.motionEnergy(dist) + comm_energy;
        req_comm_power_a = req_comm_power_b;
    end
end


