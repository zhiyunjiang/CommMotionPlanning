%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MinPNoConn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the disconnected distance along path as well as a small
% penalty for total distance. Used for fixed transmit power situation, i.e.
% the robot either can or cannot communicate. Working with partially observed
% channel, so we check that the probability of connectivity is above a
% given threshold.
%
% Inputs:
% path - array of path waypoints. Path(i) = [xi, yi].
% cawo - channel based on observations. Required channel power encoded here
% pth - threshold probability. If probbaility of conenction is above this,
%       count as connected.
% eps - small weight for total distance
% 
% Outputs:
% double total_cost - disconnected distance plus small weight for total
%                       distance

function total_cost = MinPNoConn(path, cawo, pth, eps)
    %for fixed transmit power, partially observable channel
    total_cost = 0;
    for i=2:length(path)
        no_conn_a = ( cawo.posteriorPConn(path(i-1,:)) < pth ) ;
        no_conn_b = ( cawo.posteriorPConn(path(i,:)) < pth ) ;
        dist = norm(path(i-1,:) - path(i,:));
        total_cost = total_cost + dist*(eps + 0.5*(no_conn_a + no_conn_b) );
    end
end

