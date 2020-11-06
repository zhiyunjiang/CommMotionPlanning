%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% toRawFromGrid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Translates point in discretized, normalized space to a point in
% continuous space
%
% Inputs:
% region - [x_max, x_min, y_max, y_min]
% res - resolution of the discretized region, i.e. 1/res = step size
% grid_pt - point in discretized, normalized space
%
% Output:
% raw_pt - x,y in continuous, non-shifted space
function raw_pt = toRawFromGrid(region, res, grid_pt)
    if isinf(res)
       res = 1;%we just need to shift 
    end
    
    raw_pt = (grid_pt-1)/res + [region(2), region(4)];
    
    if any( raw_pt(:,1) > region(1) ) || any( raw_pt(:,1) < region(2) )...
        || any( raw_pt(:,2) > region(3) ) || any( raw_pt(:,2) < region(4) )
        warning('Some points are outside of the region');
    end
end

