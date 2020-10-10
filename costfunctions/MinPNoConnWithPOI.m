function total_cost = MinPNoConnWithPOI(n1, n2, path, bs_cawo, poi_cawo, pth, theta, use_or)
    %for fixed transmit power, partially observable channel
    if nargin < 8
        use_or = 0;
    end
    total_cost = 0;
    if use_or
        no_conn = @(conn1, conn2) ~(conn1 || conn2);
    else
        no_conn = @(conn1, conn2) ~(conn1 && conn2);
    end
    for i=2:length(path)
        bs_conn = ( bs_cawo.posteriorPConn(path(i-1,:)) >= pth );
        poi_conn = ( poi_cawo.posteriorPConn(path(i-1,:)) >= pth );
        dist = norm(path(i-1,:) - path(i,:));
        total_cost = total_cost + dist*(theta + (1-theta)*no_conn(bs_conn, poi_conn));
    end
end

