function total_cost = MinNoConnWithPOI(n1, n2, path, bs_cc, poi_cc, gamma_th, theta, use_or)
    %for fixed transmit power, fully observable channel
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
        bs_conn = ( bs_cc.getGammaTOTdBAtPoint(path(i-1,:)) >= gamma_th) ;
        poi_conn = ( poi_cc.getGammaTOTdBAtPoint(path(i-1,:)) >= gamma_th) ;
        dist = norm(path(i-1,:) - path(i,:));
        total_cost = total_cost + dist*(theta + (1-theta)*no_conn(bs_conn, poi_conn));
    end
end

