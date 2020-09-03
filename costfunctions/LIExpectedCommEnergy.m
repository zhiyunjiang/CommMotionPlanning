function total_cost = LIExpectedCommEnergy(n1, n2, path, cawo, receiver_noise, R, K)
    %based on (5) in "Motion-Communication Co-optimization with
    %Cooperative Load Transferin Mobile Robotics: an Optimal Control Perspective"
    %For variable transmit power, partially observale chanell
    
    total_cost = 0;
    req_power_a = reqPower(cawo, receiver_noise, R, K, path(1,:));
    for i=2:length(path)
        dist = norm(path(i-1,:) - path(i,:));
        req_power_b = reqPower(cawo, receiver_noise, R, K, path(i,:));
        %scale by 1/1000 to get Joules
        total_cost = total_cost + dist*(req_power_a + req_power_b)/(2*1000);
        req_power_a = req_power_b;
    end
end

function req_power = reqPower(cawo, receiver_noise, R, K, pos)
    channel_power_dBm = cawo.posteriorExpecteddB(pos);
    CNR_lin = 10.^(channel_power_dBm/10) / receiver_noise;
    req_power = ((2^R - 1)/K)*(1/CNR_lin);
end

