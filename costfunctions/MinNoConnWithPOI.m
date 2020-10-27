%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MinNoConnWithPOI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the disconnected distance along path as well as a small
% penalty for total distance. Used for fixed transmit power situation, i.e.
% the robot either can or cannot communicate, with two remote stations.
% Working with fully observed channel.
%
% Inputs:
% path - array of path waypoints. Path(i) = [xi, yi].
% bs_cc - communication channel with the base station.
% poi_cc - communication channel for the second remote station
% gamma_th_dBm - threshold channel power in dBm. If channel is less than this,
% no connection.
% eps - small weight for total distance
% scenario - the robots tak. 1 - sensing/surveillance. 2 - broadcasting.
%                            3 - relaying
% 
% Outputs:
% double total_cost - disconnected distance plus small weight for total
%                       distance

function total_cost = MinNoConnWithPOI(path, bs_cc, poi_cc, gamma_th_dBm, eps, scenario)
    %for fixed transmit power, fully observable channel
    if nargin < 6
        scenario = 1;
    end
    total_cost = 0;
    if scenario == 1
        no_conn = @(conn1, conn2) ~(conn1 || conn2);
    else
        no_conn = @(conn1, conn2) ~(conn1 && conn2);
    end
    bs_conn_a = ( bs_cc.getGammaTOTdBAtPoint(path(1,:)) >= gamma_th_dBm);
    poi_conn_a = ( poi_cc.getGammaTOTdBAtPoint(path(1,:)) >= gamma_th_dBm);
    for i=2:length(path)
        bs_conn_b = ( bs_cc.getGammaTOTdBAtPoint(path(i,:)) >= gamma_th_dBm) ;
        poi_conn_b = ( poi_cc.getGammaTOTdBAtPoint(path(i,:)) >= gamma_th_dBm) ;
        dist = norm(path(i-1,:) - path(i,:));
        total_cost = total_cost + dist*(eps + 0.5*(no_conn(bs_conn_a, poi_conn_a) + no_conn(bs_conn_b, poi_conn_b)));
        bs_conn_a = bs_conn_b;
        poi_conn_a = poi_conn_b;
    end
end

