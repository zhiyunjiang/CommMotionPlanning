%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotPredictedField
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots the predicted communication field on which the path planning
% problem is overlayed.
%
% Inputs:
% objective - 1 indicates minimum disconnected distance (fixed TX power)
%             2 indicates minimum energy (variable TX power)
% scenario - 1 upload, 2 broadcast, 3 relay
% num_rs - number of remote stations. Can only be 1 or 2
% cawo1 - Channel analyzer with observations object, representing the
%          predicted channel from the first remote station
% cawo2 - Channel analyzer with observations object, representing the
%          predicted channel from the second remote station
% qos - quality of service parameters
% p_th - connection probability threshold
% gamma_TH - channel gain threshold

function plotPredictedField(objective, scenario, num_rs, cawo1, cawo2, qos, p_th, gamma_th)
    if objective == 1%min disconnected distance
        if num_rs == 1
            %single base station
            connection_field = cawo1.getConnectionField(p_th, gamma_th);
            cawo1.cc.plotField(connection_field, 'Connected Regions')
        elseif num_rs == 2 %
            %multi-base station, min-disconnected distance
            conn_field = cawo1.getConnectionField(p_th, gamma_th) + 2*cawo2.getConnectionField(p_th, gamma_th);
            cawo1.cc.plotDoubleConnectField(conn_field, 'Connected Regions')
        else
            error('Simulations only work with up to 2 base stations');
        end
    elseif objective == 2%min energy
        if num_rs == 1
            req_power = cawo1.getEReqTXPowerW(qos);
            
        elseif num_rs ==2
            if scenario == 1
                req_power = min(cawo1.getEReqTXPowerW(qos), cawo2.getEReqTXPowerW(qos));
            elseif scenario == 2
                req_power = max(cawo1.getEReqTXPowerW(qos), cawo2.getEReqTXPowerW(qos));
            elseif scenario == 3
                req_power = cawo1.getEReqTXPowerW(qos) + cawo2.getEReqTXPowerW(qos);
            end
        else
            error('Simulations only work with up to 2 base stations');
        end
        cawo1.cc.plotField(10*log10(req_power'), 'Required TX Power for Communication', 'Required Transmit Power (dB)');
    else
      error('Invalid objective ID');  
    end
end
