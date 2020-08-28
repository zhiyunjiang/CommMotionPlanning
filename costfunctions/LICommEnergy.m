function total_cost = LICommEnergy(n1, n2, path, cawo, receiver_noise, R, K)
    %based on (5) in "Motion-Communication Co-optimization with
    %Cooperative Load Transferin Mobile Robotics: an Optimal Control Perspective"
    %For variable transmit power
    total_cost = 0;

    req_power_a = reqPower(cawo, receiver_noise, R, K, path(1,:));
    for i=2:length(path)
        dist = norm(path(i-1,:) - path(i,:));
        req_power_b = reqPower(cawo, receiver_noise, R, K, path(i,:));
        total_cost = total_cost + dist*(req_power_a + req_power_b)/2;
        req_power_a = req_power_b;
    end
end

function req_power = reqPower(cawo, receiver_noise, R, K, pos)
    %K = -1.5/log(5*BER);
    req_power = ((2^R - 1)/K)*(1/(cawo.posteriorExpecteddB(pos) - receiver_noise));
end

