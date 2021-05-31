%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% toGridFromRaw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Translates arbitrary point to a point on discretized version of the
% space. The discretized space is 1-indexed, for compatability with matlab
% indexing.
%
% Inputs:
% region - [x_max, x_min, y_max, y_min]
% res - resolution of the discretized region, i.e. 1/res = step size
% raw_pt - arbitrary point x,y
%
% Output:
% grid_pt - x,y in grid (index) coordiantes
function grid_pt = toGridFromRaw(region, res, raw_pt)
    if any( raw_pt(:,1) > region(1) ) || any( raw_pt(:,1) < region(2) )...
        || any( raw_pt(:,2) > region(3) ) || any( raw_pt(:,2) < region(4) )
        warning('Some points are outside of the region');
    end
    
    if res == Inf% just offset so that min is 1, 1
        grid_pt = raw_pt - [region(2), region(4)] + 1;
    else
        xGMax = ceil( (region(1) - region(2))*res)+1;
        yGMax = ceil( (region(3) - region(4))*res)+1;
        xG = min(round((raw_pt(:,1) - region(2)).*res) + 1, xGMax);
        yG = min(round((raw_pt(:,2) - region(4)).*res) + 1, yGMax);
        grid_pt = [xG,yG];
    end
end

