%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getCostFxnPredicted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For given objective, scenario, and number of base stations, return the
% appropriate cost function associated with the predicted channel.
% Inputs:
% objective - 0 indicates minimum distance
%             1 indicates minimum disconnected distance (fixed TX power),
%             2 indicates minimum energy (variable TX power)
% scenario - 1 upload, 2 broadcast, 3 relay
% num_bs - number of remote stations. Can only be 1 or 2
% cawo1 - Channel Analyzer with Observations object for first remote station
% cawo2 - Channel Analyzer with Observations object for second remote station
% mp - MotionParameters object
% qos - Quality of Service parameters object
% delta - positive weight on motion path (generaly <1)
% p_th - connection probability threshold
% gamma_th - channel power threshold
%
% Output:
% cost_fxn - the cost function to be used in an instance of RRT*

function cost_fxn = getCostFxnPredicted(objective, scenario, num_bs, cawo1, cawo2, mp, qos, delta, p_th, gamma_th)

    if objective == 0%min distance cost function
        cost_fxn = @(n1, n2, path, mode) GridDist(path);
    elseif objective == 1%min disconnected distance
        if num_bs == 1
            cost_fxn = @(n1, n2, path, mode) MinPNoConn(path, cawo1, p_th, gamma_th, delta);
        elseif num_bs == 2
            cost_fxn = @(n1, n2, path, mode) MinPNoConnWithPOI(path, cawo1, cawo2, p_th, gamma_th, delta, scenario);
        else
            error('Current implementation of disconnected distance cost functions only handles up to 2 remote stations');
        end
    elseif objective == 2%min energy
        if num_bs == 1
            cost_fxn = @(n1, n2, path, mode) LIExpectedTotalEnergy(path, cawo1, qos, mp, delta);
        elseif num_bs == 2
            cost_fxn = @(n1, n2, path, mode) LIExpTotalEnergyWithPOI(path, cawo1, cawo2, qos, mp, scenario, delta);
        else
            error('Current implementation of energy cost functions only handles up to 2 remote stations');
        end
    else
        error('Invalid objective id');
    end
end
