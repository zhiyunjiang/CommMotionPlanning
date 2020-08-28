function total_cost = LIPNoConn(n1, n2, path, cawo)
    %cost is line integral of probability of no connection (linearly interpolated) 
    %for variable transmit power
    total_cost = 0;
    p_no_conn_a = 1 - cawo.posteriorPConn(path(1,:));
    for i=2:length(path)
        dist = norm(path(i-1,:) - path(i,:));
        p_no_conn_b = 1 - cawo.posteriorPConn(path(i,:));
        total_cost = total_cost + dist*(p_no_conn_a + p_no_conn_b)/2;
        p_no_conn_a = p_no_conn_b;
    end
end

