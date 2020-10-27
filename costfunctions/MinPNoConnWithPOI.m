%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MinPNoConnWithPOI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the disconnected distance along path as well as a small
% penalty for total distance. Used for fixed transmit power situation, i.e.
% the robot either can or cannot communicate, with two remote stations.
% Working with partially observed channel, so we check that the probability
% of connectivity is above a given threshold.
%
% Inputs:
% path - array of path waypoints. Path(i) = [xi, yi].
% bs_cawo - communication channel based on observations.
% poi_cawo - communication channel based on observations of second remote
%           station
% pth - threshold probability. If probbaility of conenction is above this,
%       count as connected.
% gamma_th - threshold channel power
% eps - small weight for total distance
% scenario - the robots tak. 1 - sensing/surveillance. 2 - broadcasting.
%                            3 - relaying
% 
% Outputs:
% double total_cost - disconnected distance plus small weight for total
%                       distance

function total_cost = MinPNoConnWithPOI(path, bs_cawo, poi_cawo, pth, gamma_th, eps, scenario)
    %for fixed transmit power, partially observable channel
    if nargin < 7
        scenario = 1;
    end
    total_cost = 0;
    if scenario == 1
        no_conn = @(conn1, conn2) ~(conn1 || conn2);
    else
        no_conn = @(conn1, conn2) ~(conn1 && conn2);
    end
    
    bs_conn_a = ( bs_cawo.posteriorPConn(path(1,:), gamma_th) >= pth );
    poi_conn_a = ( poi_cawo.posteriorPConn(path(1,:), gamma_th) >= pth );
    for i=2:length(path)
        bs_conn_b = ( bs_cawo.posteriorPConn(path(i,:), gamma_th) >= pth );
        poi_conn_b = ( poi_cawo.posteriorPConn(path(i,:), gamma_th) >= pth );
        dist = norm(path(i-1,:) - path(i,:));
        total_cost = total_cost + dist*...
            (eps + 0.5*( no_conn(bs_conn_a, poi_conn_a) + no_conn(bs_conn_b, poi_conn_b) ));
        bs_conn_a = bs_conn_b;
        poi_conn_a = poi_conn_b;
    end
end

