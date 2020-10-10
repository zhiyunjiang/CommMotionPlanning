function total_cost = LIExpectedCommEnergy(path, cawo, qp)
    %based on (5) in "Motion-Communication Co-optimization with
    %Cooperative Load Transferin Mobile Robotics: an Optimal Control Perspective"
    %For variable transmit power, partially observale chanell
    
    total_cost = 0;
    req_power_a = reqPower(cawo, receiver_noise, R, K, path(1,:));
    for i=2:length(path)
        dist = norm(path(i-1,:) - path(i,:));
        req_power_b = qp.reqTXPower(cawo.posteriorExpecteddB(path(i,:)));
        total_cost = total_cost + dist*(req_power_a + req_power_b)/2;
        req_power_a = req_power_b;
    end 
end

