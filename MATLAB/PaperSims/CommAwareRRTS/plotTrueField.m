%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotTrueField
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots the true communication field on which the path planning
% problem solution is later evaluated.
%
% Inputs:
% objective - 1 indicates minimum disconnected distance (fixed TX power)
%             2 indicates minimum energy (variable TX power)
% scenario - 1 upload, 2 broadcast, 3 relay
% num_rs - number of remote stations. Can only be 1 or 2
% cc1 - Comm channel object for first remote station
% cc2 - Comm channel object for second remote station
% qos - quality of service parameters
% gamma_TH - channel gain threshold
function plotTrueField(objective, scenario, num_rs, cc1, cc2, qos, gamma_th)
    if objective == 1%min disconnected distance
        if num_rs == 1
            %single base station
            connection_field = cc1.getConnectionField(gamma_th);
            cc1.plotField(connection_field, 'Connected Regions')
        elseif num_rs == 2 %
            %multi-base station, min-disconnected distance
            conn_field = cc1.getConnectionField(gamma_th) + 2*cc2.getConnectionField(gamma_th);
            cc1.plotDoubleConnectField(conn_field, 'Connected Regions')
        else
            error('Simulations only work with up to 2 base stations');
        end
    elseif objective == 2%min energy
        if num_rs == 1
            req_power = cc1.getReqTXPowerW(qos, 1);
            
        elseif num_rs ==2
            if scenario == 1
                req_power = min(cc1.getReqTXPowerW(qos, 1), cc2.getReqTXPowerW(qos, 1));
            elseif scenario == 2
                req_power = max(cc1.getReqTXPowerW(qos, 1), cc2.getReqTXPowerW(qos, 1));
            elseif scenario == 3
                req_power = cc1.getReqTXPowerW(qos, 1) + cc2.getReqTXPowerW(qos, 1);
            end
        else
            error('Simulations only work with up to 2 base stations');
        end
        cc1.plotField(10*log10(req_power), 'Required TX Power for Communication', 'Required Transmit Power (dB)');
    else
      error('Invalid objective ID');  
    end
end

