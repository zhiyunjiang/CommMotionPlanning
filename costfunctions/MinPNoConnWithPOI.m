function total_cost = MinPNoConnWithPOI(n1, n2, path, bs_cawo, poi_cawo, pth, theta)
    %for fixed transmit power, partially observable channel
    total_cost = 0;
    for i=2:length(path)
        bs_no_conn = ( bs_cawo.posteriorPConn(path(i-1,:)) < pth ) ;
        poi_no_conn = ( poi_cawo.posteriorPConn(path(i-1,:)) < pth ) ;
        no_conn = bs_no_conn * poi_no_conn;
        dist = norm(path(i-1,:) - path(i,:));
        total_cost = total_cost + dist*(theta + (1-theta)*no_conn);
    end
end

