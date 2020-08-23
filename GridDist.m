function dist = GridDist(n1,n2, path, mode)
    diffs = path(1:end-1,:) - path(2:end,:);
    dist = sum(sqrt(sum(diffs.^2,2)));
end

