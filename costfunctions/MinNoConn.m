function total_cost = MinNoConn(n1, n2, path, cc, gamma_th, theta)
    %for fixed transmit power, fully observable channel
    total_cost = 0;
    for i=2:length(path)
        no_conn = ( cc.getGammaTOTdBAtPoint(path(i-1,:)) < gamma_th) ;
        dist = norm(path(i-1,:) - path(i,:));
        total_cost = total_cost + dist*(theta + (1-theta)*no_conn);
    end
end

