function total_cost = MinPNoConn(n1, n2, path, cawo, pth, theta)
    %for fixed transmit power, partially observable channel
    total_cost = 0;
    for i=2:length(path)
        no_conn = ( cawo.posteriorPConn(path(i-1,:)) < pth ) ;
        dist = norm(path(i-1,:) - path(i,:));
        total_cost = total_cost + dist*(theta + (1-theta)*no_conn);
    end
end

