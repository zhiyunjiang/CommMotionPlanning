%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GridDist
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the Euclidean distance of a path lineraly interpolated between a
% series of waypoints. This is the standard RRT* cost function.
%
% Inputs:
% path - array of path waypoints. Path(i) = [xi, yi].
%
% Outputs:
% double dist - path's Euclidean distance

function dist = GridDist(path)
    diffs = path(1:end-1,:) - path(2:end,:);
    dist = sum(sqrt(sum(diffs.^2,2)));
end

