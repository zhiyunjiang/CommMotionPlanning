%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIExpectedTotalEnergyWithPOI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear approximation of the line integral which gives the total
% energy consumed over a path (both motion and communication). We assume
% a linear interpolation of power between two  waypoints. With multiple
% remote stations (poi), the rqeuired communication power depends on the
% robots task (scenario).
%
% Inputs:
% path - array of path waypoints. Path(i) = [xi, yi].
% cawo - channle analyzer with observations object
% cawo_poi - channle analyzer with observations object for the other remote
%             station
% qp - quality of service parameters and settings
% mp - motion parameters
% scenario - the robots tak. 1 - sensing/surveillance. 2 - broadcasting.
%                            3 - relaying
% delta - weighting parameter for motion energy. If not provided, defaults
%           to 1 (equal weight with communication energy)
% 
% Outputs:
% double energy_J - approximate energy consumed over path in Joules

function energy_J = LIExpTotalEnergyWithPOI(path, cawo, cawo_poi, qos, mp, scenario, delta)
    v_const = mp.VConst; 
    
    energy_J = 0;
    bs_req_comm_power_a = qos.reqTXPower(cawo.getMeanAtGridPoint(path(1,:)));
    poi_req_comm_power_a = qos.reqTXPower(cawo_poi.getMeanAtGridPoint(path(1,:)));
    
    for i=2:length(path)
        dist = norm(path(i-1,:) - path(i,:));
        
        bs_req_comm_power_b = qos.reqTXPower(cawo.getMeanAtGridPoint(path(i,:)));
        poi_req_comm_power_b = qos.reqTXPower(cawo_poi.getMeanAtGridPoint(path(i,:)));

        if scenario == 1
            req_comm_power_a = min([bs_req_comm_power_a, poi_req_comm_power_a]);
            req_comm_power_b = min([bs_req_comm_power_b, poi_req_comm_power_b]);
        elseif scenario == 2
            req_comm_power_a = max([bs_req_comm_power_a, poi_req_comm_power_a]);
            req_comm_power_b = max([bs_req_comm_power_b, poi_req_comm_power_b]);
        else
            req_comm_power_a = bs_req_comm_power_a + poi_req_comm_power_a;
            req_comm_power_b = bs_req_comm_power_b + poi_req_comm_power_b;
        end

        comm_energy = (dist/v_const)*(req_comm_power_a + req_comm_power_b);
        energy_J = energy_J + delta*mp.motionEnergy(dist) + comm_energy;
        
        bs_req_comm_power_a = bs_req_comm_power_b;
        poi_req_comm_power_a = poi_req_comm_power_b;
    end
end