%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MinNoConn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the disconnected distance along path as well as a small
% penalty for total distance. Used for fixed transmit power situation, i.e.
% the robot either can or cannot communicate. Working with fully observed
% channel.
%
% Inputs:
% path - array of path waypoints. Path(i) = [xi, yi].
% cc - communication channel with the base station.
% gamma_th_dBm - threshold channel power in dBm. If channel is less than this,
% no connection.
% eps - small weight for total distance
% scenario - the robots tak. 1 - sensing/surveillance. 2 - broadcasting.
%                            3 - relaying
% 
% Outputs:
% double total_cost - disconnected distance plus small weight for total
%                       distance
function total_cost = MinNoConn(path, cc, gamma_th_dBm, eps)
    %for fixed transmit power, fully observable channel
    total_cost = 0;
    for i=2:length(path)
        no_conn = ( cc.getGammaTOTdBAtPt(path(i-1,:)) < gamma_th_dBm) ;
        dist = norm(path(i-1,:) - path(i,:));
        total_cost = total_cost + dist*(eps + no_conn);
    end
end

