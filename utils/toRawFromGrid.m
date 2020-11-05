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
    gMin = 1;
    xGMax = ceil( (region(1) - region(2))*res);
    yGMax = ceil( (region(3) - region(4))*res);
    
    if any( grid_pt(:,1) > xGMax ) || any( grid_pt(:,1) < gMin )...
        || any( grid_pt(:,2) > yGMax ) || any( grid_pt(:,2) <gMin )
        warning('Some points are outside of the region');
    end
    
    raw_pt = (grid_pt-1)/res + [region(2), region(4)];
end

