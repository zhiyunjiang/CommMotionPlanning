function total_cost = MinNoConnWithPOI(n1, n2, path, bs_cc, poi_cc, gamma_th, theta)
    %for fixed transmit power, fully observable channel
    total_cost = 0;
    for i=2:length(path)
        bs_no_conn = ( bs_cc.getGammaTOTdBAtPoint(path(i-1,:)) < gamma_th) ;
        poi_no_conn = ( poi_cc.getGammaTOTdBAtPoint(path(i-1,:)) < gamma_th) ;
        no_conn = bs_no_conn * poi_no_conn;
        dist = norm(path(i-1,:) - path(i,:));
        total_cost = total_cost + dist*(theta + (1-theta)*no_conn);
    end
end

